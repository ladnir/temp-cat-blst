#pragma once
#include <cuda_runtime.h>
#include <atomic>
#include <algorithm>
#include <utility>
#include "blst/common/defines.h"

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

namespace blst {

    inline size_t ceil_div(size_t x, size_t y) { return (x + y - 1) / y; }
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








} // namespace blst
