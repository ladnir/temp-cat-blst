#include "multi_scalar_par_impl.h"

namespace blst
{

	template
		DeviceError msm_booth_batched_pipeline<Backend::CPU, G1Point, G1AffinePoint, Pointxyzz<Fp>>(
			void* d_workspace,                          // nullptr => size query
			size_t& workspace_bytes,                    // in/out
			d2Span<const G1AffinePoint> d_points,  // [M][N]
			d3Span<const uint8_t> d_scalars,            // [M][N][nbytes],
			int nbytes,
			int nbits,
			int N,                                      // points per MSM
			int M,                                      // number of MSMs
			int wbits,                                  // window bits (Booth)
			/* final outputs: */
			span<G1Point> d_outputs,                    // [M]
			Stream stream,
			bool debug,
			bool print);

#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)

	template
		DeviceError msm_booth_batched_pipeline<Backend::CUDA, G1Point, G1AffinePoint, Pointxyzz<Fp>>(
			void* d_workspace,                          // nullptr => size query
			size_t& workspace_bytes,                    // in/out
			d2Span<const G1AffinePoint> d_points,  // [M][N]
			d3Span<const uint8_t> d_scalars,            // [M][N][nbytes],
			int nbytes,
			int nbits,
			int N,                                      // points per MSM
			int M,                                      // number of MSMs
			int wbits,                                  // window bits (Booth)
			/* final outputs: */
			span<G1Point> d_outputs,                    // [M]
			Stream stream,
			bool debug,
			bool print);
#endif
//
//	template DeviceError extract_digits_flat<Backend::CPU>(
//		d2Span<const uint8_t> scalars, // [M*N][nbytes]
//		int num_scalars,                // number of scalars
//		span<WindowExtract> d_windows, // [D]
//		int nbytes,           // bytes per scalar
//		int nbits,            // significant bits per scalar
//		int wbits,            // window size (e.g., 5). W=2^wbits
//		int D,                // number of windows = ceil(nbits / wbits)
//		d2Span<u8> digit,    // [D][M*N]
//		d2Span<IdxSign> idx_sign, // [D][M*N]
//		Stream stream);
//	
//#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)
//
//	template DeviceError extract_digits_flat<Backend::CUDA>(
//		d2Span<const uint8_t> scalars, // [M*N][nbytes]
//		int num_scalars,                // number of scalars
//		span<WindowExtract> d_windows, // [D]
//		int nbytes,           // bytes per scalar
//		int nbits,            // significant bits per scalar
//		int wbits,            // window size (e.g., 5). W=2^wbits
//		int D,                // number of windows = ceil(nbits / wbits)
//		d2Span<u8> digit,    // [D][M*N]
//		d2Span<IdxSign> idx_sign, // [D][M*N]
//		Stream stream);
//#endif // BLST_CUDA

}