#include "multi_scalar_par_impl.h"

namespace blst
{

	template
	DeviceError msm_booth_batched_pipeline<Backend::CPU, G2Point, G2AffinePoint, Pointxyzz<Fp2>> (
		void* d_workspace,                          // nullptr => size query
		size_t& workspace_bytes,                    // in/out
		d2Span<const G2AffinePoint> d_points,  // [M][N]
		d3Span<const uint8_t> d_scalars,            // [M][N][nbytes],
		int nbytes,
		int nbits,
		int N,                                      // points per MSM
		int M,                                      // number of MSMs
		int wbits,                                  // window bits (Booth)
		/* final outputs: */
		span<G2Point> d_outputs,                    // [M]
		Stream stream,
		bool debug,
		bool print);

#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)

	template
		DeviceError msm_booth_batched_pipeline<Backend::CUDA, G2Point, G2AffinePoint, Pointxyzz<Fp2>>(
			void* d_workspace,                          // nullptr => size query
			size_t& workspace_bytes,                    // in/out
			d2Span<const G2AffinePoint> d_points,  // [M][N]
			d3Span<const uint8_t> d_scalars,            // [M][N][nbytes],
			int nbytes,
			int nbits,
			int N,                                      // points per MSM
			int M,                                      // number of MSMs
			int wbits,                                  // window bits (Booth)
			/* final outputs: */
			span<G2Point> d_outputs,                    // [M]
			Stream stream,
			bool debug,
			bool print);
#endif
}