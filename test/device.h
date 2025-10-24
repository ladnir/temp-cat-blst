
#include "blst/curve/multi_scalar_par.h"

namespace blst
{

	// CUDA kernel wrapper macros for device execution
#define CHECK(call) do { \
    DeviceError err = call; \
    if (err != DeviceError::cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, getErrorString<B>(err)); \
        exit(1); \
    } \
} while(0)

		// Helper function to create test points on host
	template<typename Group>
	void create_test_group(Group* points, size_t num, InsecureRng& rng) {
		Group generator = Group::generator();
		unsigned char scalar[32] = { 0 };

		for (size_t i = 0; i < num; i++) {
			scalar[0] = rng.uniform_u64();
			Group::mult(points[i], generator, scalar, 8);
		}
	}

	template<Backend B, typename T>
	auto deviceMalloc(T** ptr, size_t size) -> DeviceError {
		if constexpr (B == Backend::CPU)
		{

			*ptr = (T*)malloc(size);
			return DeviceError::cudaSuccess;
		}
		else
		{

#ifdef BLST_CUDA
			static_assert(B == Backend::CUDA);
			DeviceError err = cudaMalloc(ptr, size);
			return static_cast<DeviceError>(err);
#endif // BLST_CUDA
		}
	}

	template<Backend B>
	DeviceError deviceMemcpy(void* dest, const void* src, size_t size, MemcpyKind kind)
	{
		if constexpr (B == Backend::CPU)
		{
			memcpy(dest, src, size);
			return DeviceError::cudaSuccess;
		}
		else
		{
#ifdef BLST_CUDA
			return cudaMemcpy(dest, src, size, static_cast<cudaMemcpyKind>(kind));
#endif // BLST_CUDA

		}
	}

	template<Backend B>
	auto deviceFree(void* ptr) -> DeviceError
	{
		if constexpr (B == Backend::CPU)
		{
			free(ptr);
			return DeviceError::cudaSuccess;
		}
		else
		{
#ifdef BLST_CUDA
			DeviceError err = cudaFree(ptr);
			return static_cast<DeviceError>(err);
#endif // BLST_CUDA
		}
	}

}