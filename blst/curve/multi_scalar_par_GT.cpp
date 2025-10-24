#include "multi_scalar_par_impl.h"

namespace blst
{

	template
	DeviceError msm_booth_batched_pipeline<Backend::CPU, G1Point, G1AffinePoint, Pointxyzz<Fp>> (
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

#ifdef BLST_CUDA

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
}