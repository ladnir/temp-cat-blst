#pragma once

#include "blst/config.h"

#ifdef BLST_CUDA

#include <array>
#include <cstddef>
#include <vector>
#include <cuda_runtime.h>
#include "blst/fields/Fr.h"
#include "blst/fields/Fp12.h"
#include "blst/curve/point.h"
#include "blst/curve/GT.h"
#include "blst/fft/context.h"

namespace blst {

    // Error checking macro for CUDA calls
#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        printf("CUDA error at %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        throw std::runtime_error(cudaGetErrorString(err)); \
    } \
} while(0)

    template<typename T>
    class CudaOwnedBuffer
    {
    public:
        CudaOwnedBuffer() = default;
        ~CudaOwnedBuffer()
        {
            release();
        }

        CudaOwnedBuffer(const CudaOwnedBuffer&) = delete;
        CudaOwnedBuffer& operator=(const CudaOwnedBuffer&) = delete;

        void ensure(std::size_t count)
        {
            if (count <= capacity_) {
                return;
            }

            release();
            CUDA_CHECK(cudaMalloc(&data_, count * sizeof(T)));
            capacity_ = count;
        }

        void upload(const T* src, std::size_t count)
        {
            if (count == 0) {
                return;
            }

            ensure(count);
            CUDA_CHECK(cudaMemcpy(data_, src, count * sizeof(T), cudaMemcpyHostToDevice));
        }

        auto data() -> T*
        {
            return data_;
        }

        auto data() const -> const T*
        {
            return data_;
        }

    private:
        void release()
        {
            if (data_ != nullptr) {
                (void)cudaFree(data_);
                data_ = nullptr;
                capacity_ = 0;
            }
        }

        T* data_ = nullptr;
        std::size_t capacity_ = 0;
    };

    inline int fft_cuda_log_n(int n)
    {
        int log_n = 0;
        for (int temp = n; temp > 1; temp >>= 1) {
            log_n++;
        }
        return log_n;
    }

    inline Fr fft_cuda_inverse_scale(std::size_t n)
    {
        Fr n_fr = {};
        Fr::from_uint64(n_fr, static_cast<u64>(n));

        Fr inv_n = {};
        Fr::reciprocal(inv_n, n_fr);
        return inv_n;
    }

    template<typename G>
    struct CudaFftWorkspace
    {
        void ensure(std::size_t n)
        {
            values.ensure(n);
        }

        CudaOwnedBuffer<G> values;
    };

    template<typename G>
    struct CudaFftPreparedRoots;

    template<>
    struct CudaFftPreparedRoots<Fr>
    {
        CudaFftPreparedRoots() = default;

        explicit CudaFftPreparedRoots(const std::vector<Fr>& roots)
        {
            reset(roots);
        }

        void reset(const std::vector<Fr>& roots)
        {
            host_roots = roots;
            if (host_roots.empty()) {
                return;
            }
            device_roots.upload(host_roots.data(), host_roots.size());

            inv_n = fft_cuda_inverse_scale(host_roots.size());
            device_inv_n.upload(&inv_n, 1);
        }

        auto size() const -> std::size_t
        {
            return host_roots.size();
        }

        std::vector<Fr> host_roots;
        Fr inv_n = {};
        CudaOwnedBuffer<Fr> device_roots;
        CudaOwnedBuffer<Fr> device_inv_n;
    };

    struct CudaScalarMulPreparedRoots
    {
        CudaScalarMulPreparedRoots() = default;

        explicit CudaScalarMulPreparedRoots(const std::vector<Fr>& roots)
        {
            reset(roots);
        }

        void reset(const std::vector<Fr>& roots)
        {
            roots_non_mont.resize(roots.size());
            if (roots.empty()) {
                return;
            }
            for (std::size_t i = 0; i < roots.size(); ++i) {
                Fr::from_mont(roots_non_mont[i], roots[i]);
            }
            device_roots.upload(roots_non_mont.data(), roots_non_mont.size());

            const Fr inv_n = fft_cuda_inverse_scale(roots.size());
            Fr::from_mont(inv_n_non_mont, inv_n);
            device_inv_n.upload(&inv_n_non_mont, 1);
        }

        auto size() const -> std::size_t
        {
            return roots_non_mont.size();
        }

        std::vector<Fr> roots_non_mont;
        Fr inv_n_non_mont = {};
        CudaOwnedBuffer<Fr> device_roots;
        CudaOwnedBuffer<Fr> device_inv_n;
    };

    template<>
    struct CudaFftPreparedRoots<G1Point> : CudaScalarMulPreparedRoots
    {
        using CudaScalarMulPreparedRoots::CudaScalarMulPreparedRoots;
    };

    template<>
    struct CudaFftPreparedRoots<G2Point> : CudaScalarMulPreparedRoots
    {
        using CudaScalarMulPreparedRoots::CudaScalarMulPreparedRoots;
    };

    struct CudaGtPreparedRoots
    {
        CudaGtPreparedRoots() = default;

        explicit CudaGtPreparedRoots(const std::vector<Fr>& roots)
        {
            reset(roots);
        }

        void reset(const std::vector<Fr>& roots)
        {
            compiled_roots.resize(roots.size());
            if (roots.empty()) {
                return;
            }
            for (std::size_t i = 0; i < roots.size(); ++i) {
                Fr canonical = {};
                Fr::from_mont(canonical, roots[i]);
                compiled_roots[i] = Fp12::glv4_compile(canonical.data);
            }
            device_roots.upload(compiled_roots.data(), compiled_roots.size());

            const Fr inv_n = fft_cuda_inverse_scale(roots.size());
            Fr canonical_inv_n = {};
            Fr::from_mont(canonical_inv_n, inv_n);
            scale_compiled = Fp12::glv4_compile(canonical_inv_n.data);
            device_inv_n.upload(&scale_compiled, 1);
        }

        auto size() const -> std::size_t
        {
            return compiled_roots.size();
        }

        std::vector<GLV4_Exp> compiled_roots;
        GLV4_Exp scale_compiled = {};
        CudaOwnedBuffer<GLV4_Exp> device_roots;
        CudaOwnedBuffer<GLV4_Exp> device_inv_n;
    };

    template<>
    struct CudaFftPreparedRoots<GT> : CudaGtPreparedRoots
    {
        using CudaGtPreparedRoots::CudaGtPreparedRoots;
    };

/**
 * CUDA kernel for bit-reversal permutation.
 * Each thread handles one swap if i < j.
 */
    template<typename F>
    __global__ void bit_reverse_kernel(F* a, int n, int log_n)
    {
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) return;

        // Compute bit-reversed index
        int j = 0;
        int temp = i;
        for (int k = 0; k < log_n; ++k) {
            j = (j << 1) | (temp & 1);
            temp >>= 1;
        }

        // Only swap once (when i < j)
        if (i < j) {
            F tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
        }
    }

    /**
     * CUDA kernel for one FFT stage using Cooley-Tukey butterfly.
     * Each thread handles one butterfly operation.
     *
     * @param a         Data array (in-place)
     * @param roots     Precomputed roots of unity
     * @param n         Total size of the array
     * @param len       Current butterfly span (2, 4, 8, ..., n)
     * @param half      len / 2
     */
    template<typename G, typename F>
    __global__ void fft_butterfly_kernel(G* a, const F* roots, int n, int len, int half)
    {
        // Total number of butterflies in this stage
        int num_butterflies = n / 2;

        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        if (tid >= num_butterflies) return;

        // Determine which butterfly group this thread belongs to
        int group = tid / half;          // Which len-sized group
        int j = tid % half;              // Position within the group

        // Calculate indices into the data array
        int i = group * len;             // Start of this group
        int idx0 = i + j;
        int idx1 = i + j + half;

        // Root of unity for this butterfly
        F w = roots[j * (n / len)];

        // Cooley-Tukey butterfly: 
        //   u = a[idx0]
        //   t = a[idx1] * w
        //   a[idx0] = u + t
        //   a[idx1] = u - t
        G u = a[idx0];
        G t;

        // For Fr: use static multiply
        Fr::mul(reinterpret_cast<Fr&>(t), reinterpret_cast<const Fr&>(a[idx1]), reinterpret_cast<const Fr&>(w));

        Fr::add(reinterpret_cast<Fr&>(a[idx0]), reinterpret_cast<const Fr&>(u), reinterpret_cast<const Fr&>(t));
        Fr::sub(reinterpret_cast<Fr&>(a[idx1]), reinterpret_cast<const Fr&>(u), reinterpret_cast<const Fr&>(t));
    }

    /**
     * Specialized CUDA FFT butterfly kernel for Fr elements.
     * Uses Fr::ct_bfly which is already __host__ __device__.
     */
    __global__ void fft_butterfly_kernel_fr(Fr* a, const Fr* roots, int n, int len, int half)
    {
        int num_butterflies = n / 2;

        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        if (tid >= num_butterflies) return;

        int group = tid / half;
        int j = tid % half;

        int i = group * len;
        int idx0 = i + j;
        int idx1 = i + j + half;

        Fr w = roots[j * (n / len)];

        // Use the existing ct_bfly which is BLST_HD
        Fr::ct_bfly(a[idx0], a[idx1], w);
    }

    /**
     * CUDA FFT butterfly kernel for G1 points.
     * Performs: t = a[idx1] * w, a[idx0] = u + t, a[idx1] = u - t
     */
    __global__ void fft_butterfly_kernel_g1(G1Point* a, const Fr* roots_non_mont, int n, int len, int half)
    {
        int num_butterflies = n / 2;

        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        if (tid >= num_butterflies) return;

        int group = tid / half;
        int j = tid % half;

        int i = group * len;
        int idx0 = i + j;
        int idx1 = i + j + half;

        // Cooley-Tukey butterfly for G1:
        // t = a[idx1] * w (scalar multiplication)
        // a[idx0] = u + t
        // a[idx1] = u - t
        G1Point u = a[idx0];
        G1Point t;

        // Scalar multiplication: t = a[idx1] * w
        const Fr& w_non_mont = roots_non_mont[j * (n / len)];
        G1Point::mult(t, a[idx1], reinterpret_cast<const byte*>(&w_non_mont), 256);

        // a[idx0] = u + t
        G1Point::dadd(a[idx0], u, t);

        // a[idx1] = u - t
        G1Point neg_t;
        G1Point::neg(neg_t, t);
        G1Point::dadd(a[idx1], u, neg_t);
    }

    /**
     * CUDA FFT butterfly kernel for G2 points.
     */
    __global__ void fft_butterfly_kernel_g2(G2Point* a, const Fr* roots_non_mont, int n, int len, int half)
    {
        int num_butterflies = n / 2;

        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        if (tid >= num_butterflies) return;

        int group = tid / half;
        int j = tid % half;

        int i = group * len;
        int idx0 = i + j;
        int idx1 = i + j + half;

        G2Point u = a[idx0];
        G2Point t;

        // Scalar multiplication: t = a[idx1] * w
        const Fr& w_non_mont = roots_non_mont[j * (n / len)];
        G2Point::mult(t, a[idx1], reinterpret_cast<const byte*>(&w_non_mont), 256);

        // a[idx0] = u + t
        G2Point::dadd(a[idx0], u, t);

        // a[idx1] = u - t
        G2Point neg_t;
        G2Point::neg(neg_t, t);
        G2Point::dadd(a[idx1], u, neg_t);
    }

    /**
     * CUDA FFT butterfly kernel for GT elements.
     * GT addition is Fp12 multiplication, scalar mult is exponentiation.
     */
    __global__ void fft_butterfly_kernel_gt(GT* a, const GLV4_Exp* compiled_roots, int n, int len, int half)
    {
        int num_butterflies = n / 2;

        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        if (tid >= num_butterflies) return;

        int group = tid / half;
        int j = tid % half;

        int i = group * len;
        int idx0 = i + j;
        int idx1 = i + j + half;

        GT u = a[idx0];
        GT t = a[idx1];

        // Match the CPU FFT fast path exactly when w == 1.
        if (j != 0) {
            Fp12::cyclotomic_exp_glv4(
                t.data,
                a[idx1].data,
                compiled_roots[j * (n / len)]);
        }

        // Mirror the CPU GTCoeffCtx path exactly: subtraction is negation
        // followed by the GT group operation, rather than the GT::sub helper.
        GT::dadd(a[idx0], u, t);

        GT neg_t;
        GT::neg(neg_t, t);
        GT::dadd(a[idx1], u, neg_t);
    }

    /**
     * CUDA kernel for inverse FFT scaling (multiply each element by 1/n).
     */
    __global__ void fft_scale_kernel_fr(Fr* a, const Fr* inv_n, int n)
    {
        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        if (tid >= n) return;

        Fr temp;
        Fr::mul(temp, a[tid], *inv_n);
        a[tid] = temp;
    }

    /**
     * CUDA kernel for inverse FFT scaling for G1 points.
     */
    __global__ void fft_scale_kernel_g1(G1Point* a, const Fr* inv_n_non_mont, int n)
    {
        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        if (tid >= n) return;

        G1Point temp;
        G1Point::mult(temp, a[tid], reinterpret_cast<const byte*>(inv_n_non_mont), 256);
        a[tid] = temp;
    }

    /**
     * CUDA kernel for inverse FFT scaling for G2 points.
     */
    __global__ void fft_scale_kernel_g2(G2Point* a, const Fr* inv_n_non_mont, int n)
    {
        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        if (tid >= n) return;

        G2Point temp;
        G2Point::mult(temp, a[tid], reinterpret_cast<const byte*>(inv_n_non_mont), 256);
        a[tid] = temp;
    }

    /**
     * CUDA kernel for inverse FFT scaling for GT elements.
     */
    __global__ void fft_scale_kernel_gt(GT* a, const GLV4_Exp* inv_n, int n)
    {
        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        if (tid >= n) return;

        GT temp;
        Fp12::cyclotomic_exp_glv4(temp.data, a[tid].data, *inv_n);
        a[tid] = temp;
    }

    /**
     * Host function to perform FFT on GPU for Fr elements.
     *
     * @param a         Input/output vector (modified in place)
     * @param roots     Precomputed roots of unity (or inverse roots for iFFT)
     * @param invert    If true, perform inverse FFT and scale by 1/n
     */
    inline void fft_cuda(
        std::vector<Fr>& a,
        const CudaFftPreparedRoots<Fr>& prepared_roots,
        CudaFftWorkspace<Fr>& workspace,
        bool invert)
    {
        int n = static_cast<int>(a.size());
        if (n <= 1) return;
        BLST_ASSERT(prepared_roots.size() == a.size());

        const int log_n = fft_cuda_log_n(n);
        workspace.ensure(a.size());
        Fr* d_a = workspace.values.data();

        // Copy data to device
        CUDA_CHECK(cudaMemcpy(d_a, a.data(), n * sizeof(Fr), cudaMemcpyHostToDevice));

        // Determine block/grid sizes
        const int blockSize = 256;

        // Step 1: Bit-reversal permutation
        {
            int gridSize = (n + blockSize - 1) / blockSize;
            bit_reverse_kernel<Fr> << <gridSize, blockSize >> > (d_a, n, log_n);
            CUDA_CHECK(cudaGetLastError());
        }

        // Step 2: FFT butterfly stages
        int num_butterflies = n / 2;
        int gridSize = (num_butterflies + blockSize - 1) / blockSize;

        for (int len = 2; len <= n; len <<= 1) {
            int half = len / 2;
            fft_butterfly_kernel_fr << <gridSize, blockSize >> > (
                d_a,
                prepared_roots.device_roots.data(),
                n,
                len,
                half);
            CUDA_CHECK(cudaGetLastError());
        }

        // Step 3: Scale for inverse FFT
        if (invert) {
            int gridSizeScale = (n + blockSize - 1) / blockSize;
            fft_scale_kernel_fr << <gridSizeScale, blockSize >> > (
                d_a,
                prepared_roots.device_inv_n.data(),
                n);
            CUDA_CHECK(cudaGetLastError());
        }

        // Copy result back to host
        CUDA_CHECK(cudaMemcpy(a.data(), d_a, n * sizeof(Fr), cudaMemcpyDeviceToHost));
    }

    inline void fft_cuda(
        std::vector<Fr>& a,
        const CudaFftPreparedRoots<Fr>& prepared_roots,
        bool invert)
    {
        CudaFftWorkspace<Fr> workspace;
        fft_cuda(a, prepared_roots, workspace, invert);
    }

    inline void fft_cuda(std::vector<Fr>& a, const std::vector<Fr>& roots, bool invert)
    {
        CudaFftPreparedRoots<Fr> prepared_roots(roots);
        CudaFftWorkspace<Fr> workspace;
        fft_cuda(a, prepared_roots, workspace, invert);
    }

    /**
     * Host function to perform FFT on GPU for G1 points.
     *
     * @param a         Input/output vector of G1 points (modified in place)
     * @param roots     Precomputed roots of unity (Fr elements)
     * @param invert    If true, perform inverse FFT and scale by 1/n
     */
    inline void fft_cuda(
        std::vector<G1Point>& a,
        const CudaFftPreparedRoots<G1Point>& prepared_roots,
        CudaFftWorkspace<G1Point>& workspace,
        bool invert)
    {
        int n = static_cast<int>(a.size());
        if (n <= 1) return;
        BLST_ASSERT(prepared_roots.size() == a.size());

        const int log_n = fft_cuda_log_n(n);
        workspace.ensure(a.size());
        G1Point* d_a = workspace.values.data();

        // Copy data to device
        CUDA_CHECK(cudaMemcpy(d_a, a.data(), n * sizeof(G1Point), cudaMemcpyHostToDevice));

        // Determine block/grid sizes
        const int blockSize = 128; // Smaller for larger elements

        // Step 1: Bit-reversal permutation
        {
            int gridSize = (n + blockSize - 1) / blockSize;
            bit_reverse_kernel<G1Point> << <gridSize, blockSize >> > (d_a, n, log_n);
            CUDA_CHECK(cudaGetLastError());
        }

        // Step 2: FFT butterfly stages
        int num_butterflies = n / 2;
        int gridSize = (num_butterflies + blockSize - 1) / blockSize;

        for (int len = 2; len <= n; len <<= 1) {
            int half = len / 2;
            fft_butterfly_kernel_g1 << <gridSize, blockSize >> > (
                d_a,
                prepared_roots.device_roots.data(),
                n,
                len,
                half);
            CUDA_CHECK(cudaGetLastError());
        }

        // Step 3: Scale for inverse FFT
        if (invert) {
            int gridSizeScale = (n + blockSize - 1) / blockSize;
            fft_scale_kernel_g1 << <gridSizeScale, blockSize >> > (
                d_a,
                prepared_roots.device_inv_n.data(),
                n);
            CUDA_CHECK(cudaGetLastError());
        }

        // Copy result back to host
        CUDA_CHECK(cudaMemcpy(a.data(), d_a, n * sizeof(G1Point), cudaMemcpyDeviceToHost));
    }

    inline void fft_cuda(
        std::vector<G1Point>& a,
        const CudaFftPreparedRoots<G1Point>& prepared_roots,
        bool invert)
    {
        CudaFftWorkspace<G1Point> workspace;
        fft_cuda(a, prepared_roots, workspace, invert);
    }

    inline void fft_cuda(std::vector<G1Point>& a, const std::vector<Fr>& roots, bool invert)
    {
        CudaFftPreparedRoots<G1Point> prepared_roots(roots);
        CudaFftWorkspace<G1Point> workspace;
        fft_cuda(a, prepared_roots, workspace, invert);
    }

    /**
     * Host function to perform FFT on GPU for G2 points.
     */
    inline void fft_cuda(
        std::vector<G2Point>& a,
        const CudaFftPreparedRoots<G2Point>& prepared_roots,
        CudaFftWorkspace<G2Point>& workspace,
        bool invert)
    {
        int n = static_cast<int>(a.size());
        if (n <= 1) return;
        BLST_ASSERT(prepared_roots.size() == a.size());

        const int log_n = fft_cuda_log_n(n);
        workspace.ensure(a.size());
        G2Point* d_a = workspace.values.data();

        CUDA_CHECK(cudaMemcpy(d_a, a.data(), n * sizeof(G2Point), cudaMemcpyHostToDevice));

        const int blockSize = 64; // Even smaller for G2 (larger elements)

        // Bit-reversal
        {
            int gridSize = (n + blockSize - 1) / blockSize;
            bit_reverse_kernel<G2Point> << <gridSize, blockSize >> > (d_a, n, log_n);
            CUDA_CHECK(cudaGetLastError());
        }

        // Butterfly stages
        int num_butterflies = n / 2;
        int gridSize = (num_butterflies + blockSize - 1) / blockSize;

        for (int len = 2; len <= n; len <<= 1) {
            int half = len / 2;
            fft_butterfly_kernel_g2 << <gridSize, blockSize >> > (
                d_a,
                prepared_roots.device_roots.data(),
                n,
                len,
                half);
            CUDA_CHECK(cudaGetLastError());
        }

        // Scale for inverse FFT
        if (invert) {
            int gridSizeScale = (n + blockSize - 1) / blockSize;
            fft_scale_kernel_g2 << <gridSizeScale, blockSize >> > (
                d_a,
                prepared_roots.device_inv_n.data(),
                n);
            CUDA_CHECK(cudaGetLastError());
        }

        CUDA_CHECK(cudaMemcpy(a.data(), d_a, n * sizeof(G2Point), cudaMemcpyDeviceToHost));
    }

    inline void fft_cuda(
        std::vector<G2Point>& a,
        const CudaFftPreparedRoots<G2Point>& prepared_roots,
        bool invert)
    {
        CudaFftWorkspace<G2Point> workspace;
        fft_cuda(a, prepared_roots, workspace, invert);
    }

    inline void fft_cuda(std::vector<G2Point>& a, const std::vector<Fr>& roots, bool invert)
    {
        CudaFftPreparedRoots<G2Point> prepared_roots(roots);
        CudaFftWorkspace<G2Point> workspace;
        fft_cuda(a, prepared_roots, workspace, invert);
    }

    /**
     * Host function to perform FFT on GPU for GT elements.
     */
    inline void fft_cuda(
        std::vector<GT>& a,
        const CudaFftPreparedRoots<GT>& prepared_roots,
        CudaFftWorkspace<GT>& workspace,
        bool invert)
    {
        int n = static_cast<int>(a.size());
        if (n <= 1) return;
        BLST_ASSERT(prepared_roots.size() == a.size());

        const int log_n = fft_cuda_log_n(n);
        workspace.ensure(a.size());
        GT* d_a = workspace.values.data();

        CUDA_CHECK(cudaMemcpy(d_a, a.data(), n * sizeof(GT), cudaMemcpyHostToDevice));

        const int blockSize = 32; // Small for GT (very large elements - Fp12)

        // Bit-reversal
        {
            int gridSize = (n + blockSize - 1) / blockSize;
            bit_reverse_kernel<GT> << <gridSize, blockSize >> > (d_a, n, log_n);
            CUDA_CHECK(cudaGetLastError());
        }

        // Butterfly stages
        int num_butterflies = n / 2;
        int gridSize = (num_butterflies + blockSize - 1) / blockSize;

        for (int len = 2; len <= n; len <<= 1) {
            int half = len / 2;
            fft_butterfly_kernel_gt << <gridSize, blockSize >> > (
                d_a,
                prepared_roots.device_roots.data(),
                n,
                len,
                half);
            CUDA_CHECK(cudaGetLastError());
        }

        // Scale for inverse FFT
        if (invert) {
            int gridSizeScale = (n + blockSize - 1) / blockSize;
            fft_scale_kernel_gt << <gridSizeScale, blockSize >> > (
                d_a,
                prepared_roots.device_inv_n.data(),
                n);
            CUDA_CHECK(cudaGetLastError());
        }

        CUDA_CHECK(cudaMemcpy(a.data(), d_a, n * sizeof(GT), cudaMemcpyDeviceToHost));
    }

    inline void fft_cuda(
        std::vector<GT>& a,
        const CudaFftPreparedRoots<GT>& prepared_roots,
        bool invert)
    {
        CudaFftWorkspace<GT> workspace;
        fft_cuda(a, prepared_roots, workspace, invert);
    }

    inline void fft_cuda(std::vector<GT>& a, const std::vector<Fr>& roots, bool invert)
    {
        CudaFftPreparedRoots<GT> prepared_roots(roots);
        CudaFftWorkspace<GT> workspace;
        fft_cuda(a, prepared_roots, workspace, invert);
    }

} // namespace blst

#endif // BLST_CUDA
