#pragma once
#include <cuda_runtime.h>
#include <atomic>
#include <algorithm>
#include <utility>
#include "blst/common/defines.h"
#include <cstdio>

// -------- nodiscard that works across compilers --------
#if defined(__has_cpp_attribute)
# if __has_cpp_attribute(nodiscard)
#   define NODISCARD [[nodiscard]]
# else
#   define NODISCARD __attribute__((warn_unused_result))
# endif
#elif defined(_MSC_VER)
# define NODISCARD _Check_return_
#else
# define NODISCARD __attribute__((warn_unused_result))
#endif

#ifdef __INTELLISENSE__
#define BLST_LAUNCH(...) 
#else
#define BLST_LAUNCH(...) <<< __VA_ARGS__ >>>
#endif

namespace blst {

    inline void   clear_last_error() { (void)cudaGetLastError(); }

    // A small future-like handle for one kernel launch.
    // Not thread-safe. Move-only. Event is destroyed on destruction.
    class cuda_future {
    public:
        cuda_future() = default;

        cuda_future(cudaEvent_t evt, cudaStream_t stream, cudaError_t launch_err) noexcept
			: evt_(evt), stream_(stream), launch_err_(launch_err) {
		}

        NODISCARD cudaError_t error() const noexcept { return launch_err_; }

        // Waits until the recorded event completes and returns any async error from that work.
        // If launch failed, returns the launch error without waiting.
        NODISCARD cudaError_t get() noexcept {
            if (launch_err_ != cudaSuccess) return launch_err_;
            
            // Handle zero-size operations (no event created)
            if (!evt_) return cudaSuccess;
            
            cudaError_t e = cudaEventSynchronize(evt_);
            // Destroy event once we're done; idempotent if called again.
            cudaEventDestroy(evt_); evt_ = nullptr;
            return e;
        }

        // Non-blocking: returns true if the op has finished (and no error so far)
        bool ready() const noexcept {
            if (launch_err_ != cudaSuccess) return false;
            
            // Handle zero-size operations (no event created) - they're always "ready"
            if (!evt_) return true;
            
            return cudaEventQuery(evt_) == cudaSuccess;
        }

        // Stream used for the launch (not owned)
        cudaStream_t stream() const noexcept { return stream_; }

        // --- move support ---
        cuda_future(cuda_future&& o) noexcept { move_from(std::move(o)); }
        cuda_future& operator=(cuda_future&& o) noexcept {
            if (this != &o) { cleanup(); move_from(std::move(o)); }
            return *this;
        }

        // No copying
        cuda_future(const cuda_future&) = delete;
        cuda_future& operator=(const cuda_future&) = delete;

        ~cuda_future() { cleanup(); }

    private:
        void cleanup() noexcept {
            if (evt_) cudaEventDestroy(evt_);
            evt_ = nullptr;
        }
        void move_from(cuda_future&& o) noexcept {
            evt_ = o.evt_; o.evt_ = nullptr;
            stream_ = o.stream_;
            launch_err_ = o.launch_err_;
        }

        cudaEvent_t  evt_ = nullptr;
        cudaStream_t stream_ = 0;
        cudaError_t  launch_err_ = cudaSuccess;

        friend struct launcher_access;
    };

    // ---- runtime-tuned launcher: works with any non-templated __global__ kernel ----
    template <typename Kernel, typename... Args>
    NODISCARD cuda_future launch_auto(Kernel kernel,
        size_t n,
        Args&&... args)
    {
		int block_limit = 0;
		size_t dynamic_smem = 0;
		cudaStream_t stream = 0;

        // Handle zero-size operations - return successful no-op
        if (n == 0) {
            return cuda_future(nullptr, stream, cudaSuccess);
        }

        // Cache a good TB once per (Kernel) instantiation
        static std::atomic<int> cached_tb{ 0 };
        int tb = cached_tb.load(std::memory_order_acquire);

        if (tb <= 0) {
            int minGrid = 0, sugg = 0;
            cudaError_t occ = cudaOccupancyMaxPotentialBlockSize(
                &minGrid, &sugg, kernel, dynamic_smem, block_limit);
            if (occ != cudaSuccess || sugg <= 0) sugg = 128;
            sugg = (sugg + 31) & ~31; // warp multiple
            cached_tb.store(sugg, std::memory_order_release);
            tb = sugg;
        }

        // Respect kernel's cap
        cudaFuncAttributes fa{};
        if (cudaFuncGetAttributes(&fa, (const void*)kernel) == cudaSuccess)
            tb = std::min(tb, fa.maxThreadsPerBlock);

        dim3 block((unsigned)tb, 1, 1);
        dim3 grid((unsigned)ceil_div(n, (size_t)block.x), 1, 1);

        clear_last_error(); // clear any sticky error from earlier

        // Create an event we record right after the launch (lightweight: no timing)
		cudaEvent_t evt = nullptr;
		cudaError_t launch_err = cudaEventCreateWithFlags(&evt, cudaEventDisableTiming);
        if (launch_err != cudaSuccess) {
			return cuda_future(nullptr, stream, launch_err);
        }

        // Launch
#ifndef __INTELLISENSE__ // suppress VS IntelliSense error
        kernel << <grid, block, dynamic_smem, stream >> > (std::forward<Args>(args)..., n);
#endif

        // Capture immediate launch/config error; if OK, record the event after this kernel
        launch_err = cudaGetLastError();
        if (launch_err == cudaSuccess) {
            cudaError_t rec = cudaEventRecord(evt, stream);
            if (rec != cudaSuccess) { launch_err = rec; }
        }
        // If launch failed, cleanup event now so get() is trivial
        if (launch_err != cudaSuccess) { cudaEventDestroy(evt); evt = nullptr; }
		return cuda_future(evt, stream, launch_err);
    }

    // ---- optional: compile-time TB specialization path ----------------------------
    // Write your kernel as template<int TB> with __launch_bounds__(TB):
    //   template<int TB>
    //   __launch_bounds__(TB)
    //   __global__ void my_kernel_tb(..., size_t n) { /* grid-stride */ }
    //
    // Then dispatch by TB:
    template<int TB, typename KernelT, typename... Args>
    NODISCARD cuda_future launch_tb(KernelT kernel, size_t n, Args&&... args,
        size_t dynamic_smem = 0, cudaStream_t stream = 0)
    {
        // Handle zero-size operations - return successful no-op
        if (n == 0) {
            return cuda_future(nullptr, stream, cudaSuccess);
        }

        dim3 block((unsigned)TB, 1, 1);
        dim3 grid((unsigned)ceil_div(n, (size_t)TB), 1, 1);

        clear_last_error();

        cudaEvent_t evt = nullptr;
        auto launch_err = cudaEventCreateWithFlags(&evt, cudaEventDisableTiming);
        if (launch_err != cudaSuccess) {
            return cuda_future(nullptr, stream, launch_err);
		}

#ifndef __INTELLISENSE__ // suppress VS IntelliSense error
        kernel << <grid, block, dynamic_smem, stream >> > (std::forward<Args>(args)..., n);
#endif

		// Capture immediate launch/config error; if OK, record the event after this kernel
        launch_err = cudaGetLastError();
        if (launch_err == cudaSuccess) {
            cudaError_t rec = cudaEventRecord(evt, stream);
            if (rec != cudaSuccess) launch_err = rec;
        }
        if (launch_err != cudaSuccess) { cudaEventDestroy(evt); evt = nullptr; }
        
		return cuda_future(evt, stream, launch_err);
    }











    inline const char* mtype_str(int t)
    {
        switch (t)
        {
        case 1: return "host";
        case 2: return "device";
        case 3: return "managed";
        case 4: return "unified";
        default: return "unknown";
        }
    }

    inline void print_dev_limits(int dev)
    {
        int maxThreadsPerBlock = 0;
        int maxGridX = 0;
        int maxGridY = 0;
        int maxGridZ = 0;
        int maxSmemLegacy = 0;
        int maxSmemOptin = 0;
        cudaDeviceGetAttribute(&maxThreadsPerBlock, cudaDevAttrMaxThreadsPerBlock, dev);
        cudaDeviceGetAttribute(&maxGridX, cudaDevAttrMaxGridDimX, dev);
        cudaDeviceGetAttribute(&maxGridY, cudaDevAttrMaxGridDimY, dev);
        cudaDeviceGetAttribute(&maxGridZ, cudaDevAttrMaxGridDimZ, dev);
        cudaDeviceGetAttribute(&maxSmemLegacy, cudaDevAttrMaxSharedMemoryPerBlock, dev);
        cudaDeviceGetAttribute(&maxSmemOptin, cudaDevAttrMaxSharedMemoryPerBlockOptin, dev);
        std::fprintf(stderr,
            "[device %d] maxTPB=%d gridMax=(%d,%d,%d) maxSmemLegacy=%d maxSmemOptin=%d\n",
            dev, maxThreadsPerBlock, maxGridX, maxGridY, maxGridZ, maxSmemLegacy, maxSmemOptin);
    }

    template <typename T>
    inline void check_ptr_arg(const T&) {}

    template <typename T>
    inline void check_ptr_arg(T* p)
    {
        cudaPointerAttributes a;
        cudaError_t s = cudaPointerGetAttributes(&a, (const void*)p);
        if (s == cudaSuccess)
        {
            std::fprintf(stderr, "  arg %p -> %s\n", (const void*)p, mtype_str((int)a.type));
        }
        else
        {
            std::fprintf(stderr, "  arg %p -> cudaPointerGetAttributes error %d\n", (const void*)p, (int)s);
        }
    }

    template <typename... Args>
    inline void check_all_ptr_args(const Args&... args)
    {
        int dummy[] = { (check_ptr_arg(args), 0)... };
        (void)dummy;
    }

    template <typename KernelT>
    inline void print_kernel_attrs(KernelT k, const char* name)
    {
        cudaFuncAttributes fa{};
        if (cudaFuncGetAttributes(&fa, (const void*)k) == cudaSuccess)
        {
            std::fprintf(stderr,
                "[ker %s] regs=%d staticSmem=%zu local=%zu maxTPB=%d bin=%d ptx=%d\n",
                name, fa.numRegs, (size_t)fa.sharedSizeBytes, (size_t)fa.localSizeBytes,
                fa.maxThreadsPerBlock, fa.binaryVersion, fa.ptxVersion);
        }
        else
        {
            std::fprintf(stderr, "[ker %s] cudaFuncGetAttributes failed\n", name);
        }
    }

    template <typename KernelT>
    inline int suggest_block_size(KernelT k, size_t dyn_smem)
    {
        int minGrid = 0;
        int block = 0;
        // If dyn_smem varies with block, you can pass a lambda via the VariableSMem API.
        if (cudaOccupancyMaxPotentialBlockSize(&minGrid, &block, k, (size_t)dyn_smem, 0) == cudaSuccess && block > 0)
        {
            return block;
        }
        return 128; // fallback guess
    }

    // NEW: clamp is templated on the kernel so we can set attributes on the right symbol.
    template <typename KernelT>
    inline cudaError_t clamp_dims_to_limits(KernelT k,
        dim3& grid,
        dim3& block,
        size_t& dyn_smem,
        const cudaFuncAttributes& fa,
        int dev)
    {
        int maxTPBDev = 0, maxGridX = 0, maxGridY = 0, maxGridZ = 0;
        int maxSmemLegacy = 0, maxSmemOptin = 0;

        cudaDeviceGetAttribute(&maxTPBDev, cudaDevAttrMaxThreadsPerBlock, dev);
        cudaDeviceGetAttribute(&maxGridX, cudaDevAttrMaxGridDimX, dev);
        cudaDeviceGetAttribute(&maxGridY, cudaDevAttrMaxGridDimY, dev);
        cudaDeviceGetAttribute(&maxGridZ, cudaDevAttrMaxGridDimZ, dev);
        cudaDeviceGetAttribute(&maxSmemLegacy, cudaDevAttrMaxSharedMemoryPerBlock, dev);
        cudaDeviceGetAttribute(&maxSmemOptin, cudaDevAttrMaxSharedMemoryPerBlockOptin, dev);

        // Clamp block to both device and kernel limits
        const int tpb = (int)block.x * (int)block.y * (int)block.z;
        const int maxTPBKer = (fa.maxThreadsPerBlock > 0 ? fa.maxThreadsPerBlock : maxTPBDev);
        if (tpb <= 0 || tpb > maxTPBKer || tpb > maxTPBDev) {
            int clamped = maxTPBKer < maxTPBDev ? maxTPBKer : maxTPBDev;
            if (clamped < 1) clamped = 1;
            block = dim3((unsigned)clamped, 1, 1);
            std::fprintf(stderr, "[clamp] block -> (%u,%u,%u)\n", block.x, block.y, block.z);
        }

        // Clamp grid
        if (grid.x == 0) grid.x = 1; if (grid.y == 0) grid.y = 1; if (grid.z == 0) grid.z = 1;
        if (grid.x > (unsigned)maxGridX) grid.x = maxGridX;
        if (grid.y > (unsigned)maxGridY) grid.y = maxGridY;
        if (grid.z > (unsigned)maxGridZ) grid.z = maxGridZ;

        // Dynamic shared memory opt-in (only if requested)
        if (dyn_smem > 0) {
            int cap = (maxSmemOptin > 0 ? maxSmemOptin : maxSmemLegacy);
            int want = (int)dyn_smem;
            if (want > cap) want = cap;

            // If we want > legacy, try to opt-in on this kernel symbol
            if (want > maxSmemLegacy) {
                cudaError_t s = cudaFuncSetAttribute((const void*)k,
                    cudaFuncAttributeMaxDynamicSharedMemorySize,
                    want);
                if (s != cudaSuccess) {
                    std::fprintf(stderr, "[clamp] opt-in failed (%s); using legacy %d\n",
                        cudaGetErrorString(s), maxSmemLegacy);
                    dyn_smem = (size_t)maxSmemLegacy;
                }
                else {
                    dyn_smem = (size_t)want;
                }
            }
            else {
                // Stay within legacy cap
                if ((int)dyn_smem > maxSmemLegacy) dyn_smem = (size_t)maxSmemLegacy;
            }
        }

        return cudaSuccess;
    }

    template <typename KernelT, typename... Args>
    inline cudaError_t launch_checked(KernelT k,
        dim3 grid,
        dim3 block,
        size_t dyn_smem,
        cudaStream_t stream,
        const char* name,
        const Args&... args)
    {
        int dev = 0; cudaGetDevice(&dev);
        print_dev_limits(dev);
        print_kernel_attrs(k, name);
        cudaDeviceSetLimit(cudaLimitStackSize, 64 * 1024);

        std::fprintf(stderr, "[args ptr types]\n");
        check_all_ptr_args(args...);

        cudaFuncAttributes fa{};
        cudaFuncGetAttributes(&fa, (const void*)k);

        // NOTE: pass the kernel symbol to clamp so it can set func attributes on the right thing
        clamp_dims_to_limits(k, grid, block, dyn_smem, fa, dev);

        std::fprintf(stderr, "[launch %s] grid=(%u,%u,%u) block=(%u,%u,%u) dyn_smem=%zu\n",
            name, grid.x, grid.y, grid.z, block.x, block.y, block.z, dyn_smem);

        k << <grid, block, dyn_smem, stream >> > (args...);
        cudaError_t e = cudaGetLastError();
        if (e == cudaSuccess) return e;

        std::fprintf(stderr, "[launch err %s] %s (%d)\n", name, cudaGetErrorString(e), (int)e);

        // Fallback: occupancy-suggested block and tiny grid
        int sugg = suggest_block_size(k, dyn_smem);
        if (sugg < 32) sugg = 32;
        std::fprintf(stderr, "[fallback] suggest block=%d, retry grid=1 block=%d dyn_smem=0\n", sugg, sugg);
        k << <dim3(1, 1, 1), dim3((unsigned)sugg, 1, 1), 0, stream >> > (args...);
        cudaError_t e2 = cudaGetLastError();
        std::fprintf(stderr, "[fallback result] %s (%d)\n", cudaGetErrorString(e2), (int)e2);
        return e;
    }




} // namespace blst
