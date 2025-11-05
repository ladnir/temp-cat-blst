#include <benchmark/benchmark.h>
#include <random>
#include <cstring>
#include <stdint.h>

// Include the actual BLST field headers
#include "blst/fields/Fp.h"
#include "blst/fields/Fr.h"
#include "blst/fields/Fp_cuda.h"
#include "blst/backend/backend.h"
#include "blst/backend/cuda/mul_mont_cuda.h"


using namespace blst;

#define LIMBS NLIMBS(384)

// Helper function to set Fp from uint64_t value
void set_Fp(Fp& v, uint64_t x) {
	for (int i = 0; i < LIMBS; ++i)
		v[i] = 0;
	memcpy(&v, &x, sizeof(x));
}

// Helper function to generate random Fp elements
void generate_random_fp(Fp& fp) {
	static std::random_device rd;
	static std::mt19937_64 gen(rd());
	static std::uniform_int_distribution<uint64_t> dis;

	// Initialize with random values
	for (int i = 0; i < LIMBS; ++i) {
		fp[i] = dis(gen);
	}

	while (Fp::is_redc(fp) == 0)
	{
		Fp::rshift(fp, fp, 1);
	}
}

// ========== ADDITION BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp addition (ASM version)
static void BM_FpAddition_asm(benchmark::State& state) {
	Fp a, b, c;
	generate_random_fp(a);
	generate_random_fp(b);
	auto p = get_BLST_P();

	for (auto _ : state) {
		add_mod_384(&c, &a, &b, &p);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpAddition_asm);
#endif

// Benchmark for Fp addition (C version)
static void BM_FpAddition_cc(benchmark::State& state) {
	Fp a, b, c;
	generate_random_fp(a);
	generate_random_fp(b);
	auto p = get_BLST_P();

	for (auto _ : state) {
		//Fp::add(c, a, b);
		cc_add_mod_n<NLIMBS(384)>(c.data, a.data, b.data, get_BLST_P().data);

		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}


// Benchmark for Fp addition (High-level API)
static void BM_FpAddition_api(benchmark::State& state) {
	Fp a, b, c;
	generate_random_fp(a);
	generate_random_fp(b);

	for (auto _ : state) {
		Fp::add(c, a, b);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpAddition_api);

// ========== SUBTRACTION BENCHMARKS ==========

#ifdef BLST_ASM

// Benchmark for Fp subtraction (ASM version)
static void BM_FpSubtraction_asm(benchmark::State& state) {
	Fp a, b, c;
	generate_random_fp(a);
	generate_random_fp(b);
	auto p = get_BLST_P();

	for (auto _ : state) {
		sub_mod_384(&c, &a, &b, &p);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSubtraction_asm);

#endif

// Benchmark for Fp subtraction (C version)
static void BM_FpSubtraction_cc(benchmark::State& state) {
	Fp a, b, c;
	generate_random_fp(a);
	generate_random_fp(b);
	auto p = get_BLST_P();

	for (auto _ : state) {
		cc_sub_mod_n<NLIMBS(384)>(c.data, a.data, b.data, p.data);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSubtraction_cc);

// Benchmark for Fp subtraction (High-level API)
static void BM_FpSubtraction_api(benchmark::State& state) {
	Fp a, b, c;
	generate_random_fp(a);
	generate_random_fp(b);

	for (auto _ : state) {
		Fp::sub(c, a, b);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSubtraction_api);

// ========== MULTIPLICATION BENCHMARKS ==========
#ifdef BLST_ASM

// Benchmark for Fp multiplication (ASM version)
static void BM_FpMultiplication_asm(benchmark::State& state) {
	Fp a, b, c;
	generate_random_fp(a);
	generate_random_fp(b);

	// Convert to Montgomery form for multiplication
	Fp::to_mont(a, a);
	Fp::to_mont(b, b);

	auto p = get_BLST_P();

	for (auto _ : state) {
		mulx_mont_384(&c, &a, &b, &p, BLST_p0);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMultiplication_asm);
#endif

// Benchmark for Fp multiplication (C version)
static void BM_FpMultiplication_cc(benchmark::State& state) {
	Fp a, b, c;
	generate_random_fp(a);
	generate_random_fp(b);

	// Convert to Montgomery form for multiplication
	Fp::to_mont(a, a);
	Fp::to_mont(b, b);

	//auto p = get_BLST_P();

	for (auto _ : state) {
		//Fp::mul(c, a, b);
		cc_mul_mont_n<NLIMBS(384)>(c.data, a.data, b.data, get_BLST_P().data, BLST_p0);

		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMultiplication_cc);

// Benchmark for Fp multiplication (High-level API)
static void BM_FpMultiplication_api(benchmark::State& state) {
	Fp a, b, c;
	generate_random_fp(a);
	generate_random_fp(b);

	// Convert to Montgomery form for multiplication
	Fp::to_mont(a, a);
	Fp::to_mont(b, b);

	for (auto _ : state) {
		Fp::mul(c, a, b);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMultiplication_api);

// ========== SQUARING BENCHMARKS ==========
#ifdef BLST_ASM

// Benchmark for Fp squaring (ASM version)
static void BM_FpSquaring_asm(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);

	// Convert to Montgomery form for squaring
	Fp::to_mont(a, a);

	auto p = get_BLST_P();

	for (auto _ : state) {
		sqrx_mont_384(&c, &a, &p, BLST_p0);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSquaring_asm);
#endif

// Benchmark for Fp squaring (C version)
static void BM_FpSquaring_cc(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);

	// Convert to Montgomery form for squaring
	Fp::to_mont(a, a);

	auto p = get_BLST_P();

	for (auto _ : state) {
		//Fp::sqr(c, a);
		cc_mul_mont_n<NLIMBS(384)>(c.data, a.data, a.data, get_BLST_P().data, BLST_p0);

		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSquaring_cc);

// Benchmark for Fp squaring (High-level API)
static void BM_FpSquaring_api(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);

	// Convert to Montgomery form for squaring
	Fp::to_mont(a, a);

	for (auto _ : state) {
		Fp::sqr(c, a);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSquaring_api);

// ========== MULTIPLY BY 3 BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp multiplication by 3 (ASM version)
static void BM_FpMulBy3_asm(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);
	auto p = get_BLST_P();

	for (auto _ : state) {
		mul_by_3_mod_384(&c, &a, &p);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy3_asm);
#endif

// Benchmark for Fp multiplication by 3 (C version)
static void BM_FpMulBy3_cc(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);
	auto p = get_BLST_P();

	for (auto _ : state) {
		//Fp::mul_by_3(c, a);
		cc_mul_by_3_mod_n<NLIMBS(384)>(c.data, a.data, get_BLST_P().data);

		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy3_cc);

// Benchmark for Fp multiplication by 3 (High-level API)
static void BM_FpMulBy3_api(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);

	for (auto _ : state) {
		Fp::mul_by_3(c, a);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy3_api);

// ========== MULTIPLY BY 8 BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp multiplication by 8 (ASM version)
static void BM_FpMulBy8_asm(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);
	auto p = get_BLST_P();

	for (auto _ : state) {
		mul_by_8_mod_384(&c, &a, &p);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy8_asm);
#endif

// Benchmark for Fp multiplication by 8 (C version)
static void BM_FpMulBy8_cc(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);

	for (auto _ : state) {
		Fp::mul_by_8(c, a);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy8_cc);

// Benchmark for Fp multiplication by 8 (High-level API)
static void BM_FpMulBy8_api(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);

	for (auto _ : state) {
		Fp::mul_by_8(c, a);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy8_api);

// ========== RIGHT SHIFT BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp right shift (ASM version)
static void BM_FpRshift_asm(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);
	auto p = get_BLST_P();

	for (auto _ : state) {
		rshift_mod_384(&c, &a, 2, &p);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpRshift_asm);
#endif

// Benchmark for Fp right shift (C version)
static void BM_FpRshift_cc(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);
	auto p = get_BLST_P();

	for (auto _ : state) {
		//Fp::rshift(c, a, 2);
		cc_rshift_mod_n<NLIMBS(384)>(c.data, a.data, 2, get_BLST_P().data);

		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpRshift_cc);

// Benchmark for Fp right shift (High-level API)
static void BM_FpRshift_api(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);

	for (auto _ : state) {
		Fp::rshift(c, a, 2);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpRshift_api);



// ========== LEFT  SHIFT BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp right shift (ASM version)
static void BM_FpLshift_asm(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);
	auto p = get_BLST_P();

	for (auto _ : state) {
		lshift_mod_384(&c, &a, 2, &p);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpLshift_asm);
#endif

// Benchmark for Fp left  shift (C version)
static void BM_FpLshift_cc(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);
	auto p = get_BLST_P();

	for (auto _ : state) {
		//Fp::rshift(c, a, 2);
		cc_lshift_mod_n<NLIMBS(384)>(c.data, a.data, 2, get_BLST_P().data);

		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpLshift_cc);

// Benchmark for Fp left  shift (High-level API)
static void BM_FpLshift_api(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);

	for (auto _ : state) {
		Fp::lshift(c, a, 2);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpLshift_api);


// ========== DIVISION BY 2 BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp division by 2 (ASM version)
static void BM_FpDivBy2_asm(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);
	auto p = get_BLST_P();

	for (auto _ : state) {
		div_by_2_mod_384(&c, &a, &p);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpDivBy2_asm);
#endif

// Benchmark for Fp division by 2 (C version)
static void BM_FpDivBy2_cc(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);
	auto p = get_BLST_P();

	for (auto _ : state) {
		//cc_div_by_2_mod_384(c, a, p);
		cc_rshift_mod_n<NLIMBS(384)>(c.data, a.data, 1, get_BLST_P().data);

		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpDivBy2_cc);

// Benchmark for Fp division by 2 (High-level API)
static void BM_FpDivBy2_api(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);

	for (auto _ : state) {
		Fp::div_by_2(c, a);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpDivBy2_api);

// ========== MONTGOMERY CONVERSION BENCHMARKS ==========


#ifdef BLST_ASM
// Benchmark for Montgomery conversion from (ASM version)
static void BM_FpFromMontgomery_asm(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);
	Fp::to_mont(a, a);  // Convert to Montgomery form first
	auto p = get_BLST_P();

	for (auto _ : state) {
		fromx_mont_384(&c, &a, &p, BLST_p0);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpFromMontgomery_asm);
#endif

// Benchmark for Montgomery conversion from (C version)
static void BM_FpFromMontgomery_cc(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);
	Fp::to_mont(a, a);  // Convert to Montgomery form first
	auto p = get_BLST_P();

	for (auto _ : state) {
		//cc_from_mont_384(c, a, p, BLST_p0);
		cc_from_mont_n<NLIMBS(384)>(c.data, a.data, get_BLST_P().data, BLST_p0);

		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpFromMontgomery_cc);

// Benchmark for Montgomery conversion from (High-level API)
static void BM_FpFromMontgomery_api(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);
	Fp::to_mont(a, a);  // Convert to Montgomery form first

	for (auto _ : state) {
		Fp::from_mont(c, a);
		benchmark::DoNotOptimize(c);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpFromMontgomery_api);

// ========== CONDITIONAL NEGATION BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp conditional negation (ASM version)
static void BM_FpCneg_asm(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);
	bool flag = true;
	auto p = get_BLST_P();

	for (auto _ : state) {
		cneg_mod_384(&c, &a, flag, &p);
		benchmark::DoNotOptimize(c);
		flag = !flag;
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpCneg_asm);
#endif

// Benchmark for Fp conditional negation (C version)
static void BM_FpCneg_cc(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);
	bool flag = true;
	auto p = get_BLST_P();

	for (auto _ : state) {
		//cc_cneg_mod_384(c, a, flag, p);
		cc_cneg_mod_n<NLIMBS(384)>(c.data, a.data, flag, p.data);

		benchmark::DoNotOptimize(c);
		flag = !flag;
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpCneg_cc);

// Benchmark for Fp conditional negation (High-level API)
static void BM_FpCneg_api(benchmark::State& state) {
	Fp a, c;
	generate_random_fp(a);
	bool flag = true;

	for (auto _ : state) {
		Fp::cneg(c, a, flag);
		benchmark::DoNotOptimize(c);
		flag = !flag;
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpCneg_api);

// ========== EXPONENTIATION BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for exponentiation using sqr_n_mul (ASM version)
static void BM_FpExponentiation_asm(benchmark::State& state) {
	Fp a, result, one;
	set_Fp(a, 3);  // Base value
	Fp::to_mont(a, a);       // Convert to Montgomery form
	set_Fp(one, 1);
	Fp::to_mont(one, one);   // One in Montgomery form

	auto p = get_BLST_P();

	for (auto _ : state) {
		// Use sqr_n_mul_mont_383 for efficient exponentiation
		// This computes a^(2^count) * b
		sqrx_n_mul_mont_383(&result, &a, 3, &p, BLST_p0, &one); // a^8 * 1 = a^8
		benchmark::DoNotOptimize(result);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpExponentiation_asm);
#endif

// Benchmark for exponentiation using sqr_n_mul (C version)
static void BM_FpExponentiation_cc(benchmark::State& state) {
	Fp a, result, one;
	set_Fp(a, 3);  // Base value
	Fp::to_mont(a, a);       // Convert to Montgomery form
	set_Fp(one, 1);
	Fp::to_mont(one, one);   // One in Montgomery form

	auto p = get_BLST_P();

	for (auto _ : state) {
		// Use C version of sqr_n_mul_mont
		cc_sqr_n_mul_mont_384(result.data, a.data, 3, p.data, BLST_p0, one.data); // a^8 * 1 = a^8
		benchmark::DoNotOptimize(result);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpExponentiation_cc);

// Benchmark for exponentiation using sqr_n_mul (High-level API)
static void BM_FpExponentiation_api(benchmark::State& state) {
	Fp a, result, one;
	set_Fp(a, 3);  // Base value
	Fp::to_mont(a, a);       // Convert to Montgomery form
	set_Fp(one, 1);
	Fp::to_mont(one, one);   // One in Montgomery form

	for (auto _ : state) {
		// Use high-level API
		Fp::sqr_n_mul(result, a, 3, one); // a^8 * 1 = a^8
		benchmark::DoNotOptimize(result);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpExponentiation_api);

// ========== INVERSION BENCHMARKS ==========

// Benchmark for Fp inversion (fast reciprocal)
static void BM_FpInversion_fast(benchmark::State& state) {
	Fp a, result;
	generate_random_fp(a);
	Fp::to_mont(a, a);  // Convert to Montgomery form

	for (auto _ : state) {
		Fp::flt_reciprocal(result, a);
		benchmark::DoNotOptimize(result);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpInversion_fast);

// Benchmark for Fp inversion (constant-time reciprocal)
static void BM_FpInversion_ct(benchmark::State& state) {
	Fp a, result;
	generate_random_fp(a);
	Fp::to_mont(a, a);  // Convert to Montgomery form

	for (auto _ : state) {
		Fp::reciprocal(result, a);
		benchmark::DoNotOptimize(result);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpInversion_ct);

// ========== COMPLEX OPERATIONS BENCHMARKS ==========

// Benchmark for complex operation combining multiple operations
static void BM_FpComplexOperation(benchmark::State& state) {
	Fp a, b, c, d, temp1, temp2;
	generate_random_fp(a);
	generate_random_fp(b);

	// Convert to Montgomery form
	Fp::to_mont(a, a);
	Fp::to_mont(b, b);

	for (auto _ : state) {
		// Perform: d = (a * b) + (a^2) - b
		Fp::mul(temp1, a, b);      // temp1 = a * b
		Fp::sqr(c, a);             // c = a^2  
		Fp::add(temp2, temp1, c);  // temp2 = temp1 + c
		Fp::sub(d, temp2, b);      // d = temp2 - b
		benchmark::DoNotOptimize(d);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpComplexOperation);

// Benchmark for Montgomery round-trip (Fp::to_mont + from_fp)
static void BM_FpMontgomeryRoundTrip(benchmark::State& state) {
	Fp a, mont, result;
	generate_random_fp(a);

	for (auto _ : state) {
		Fp::to_mont(mont, a);         // Convert to Montgomery form
		Fp::from_mont(result, mont);  // Convert back from Montgomery form
		benchmark::DoNotOptimize(result);
	}

	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMontgomeryRoundTrip);









// ==== small cuda helpers ====
#define CUDA_OK(x) do { cudaError_t _e = (x); if (_e != cudaSuccess) { \
  fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(_e), __FILE__, __LINE__); std::abort(); } } while (0)

// ==== __global__ kernels (1 element per thread, strided) ====
BLST_GLOBAL void kernel_add_mod_384_u64
(
	uint64_t* BLST_restrict out,
	const uint64_t* BLST_restrict a,
	const uint64_t* BLST_restrict b,
	size_t n_elems
)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	size_t stride = gridDim.x * blockDim.x;
	for (size_t i = idx; i < n_elems; i += stride)
	{
		const uint64_t* pa = a + 6 * i;
		const uint64_t* pb = b + 6 * i;
		uint64_t* pr = out + 6 * i;
		ptx_add_mod_u64<Fp>(pr, pa, pb);
	}
}

#ifdef BLST_ADD_MOD_384_U32_PTX
BLST_GLOBAL void kernel_add_mod_384_u32
(
	uint32_t* BLST_restrict out,
	const uint32_t* BLST_restrict a,
	const uint32_t* BLST_restrict b,
	size_t n_elems
)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	size_t stride = gridDim.x * blockDim.x;
	for (size_t i = idx; i < n_elems; i += stride)
	{
		const uint32_t* pa = a + 12 * i;
		const uint32_t* pb = b + 12 * i;
		uint32_t* pr = out + 12 * i;
		add_mod_384_u32_ptx(pr, pa, pb, (u32*)Fp::modulus().data);
	}
}
#endif

// ==== Benchmark fixtures (Google Benchmark) ====
// We use manual timing with CUDA events so we measure device time only.
// Memory traffic per element ~= 144 bytes (read a,b = 96B; write r = 48B). We expose that as a counter.

struct GpuBuf64 {
	uint64_t* d_a = nullptr;
	uint64_t* d_b = nullptr;
	uint64_t* d_r = nullptr;
	size_t n = 0;
};

struct GpuBuf32 {
	uint32_t* d_a = nullptr;
	uint32_t* d_b = nullptr;
	uint32_t* d_r = nullptr;
	size_t n = 0;
};

struct GpuBufFp {
	Fp* d_a = nullptr;
	Fp* d_b = nullptr;
	Fp* d_r = nullptr;
	size_t n = 0;
};

static void fill_host_u64(std::vector<uint64_t>& v)
{
	std::mt19937_64 rng(12345);
	for (auto& x : v) x = rng();
}

static void fill_host_u32(std::vector<uint32_t>& v)
{
	std::mt19937 rng(12345);
	for (auto& x : v) x = (uint32_t)rng();
}

static void fill_host_Fp(std::vector<Fp>& v)
{
	std::mt19937 rng(12345);
	for (auto& X : v)
		for (u64 i = 0; i < NLIMBS(384); ++i)
			X[i] = (uint32_t)rng();
}


static void choose_launch(size_t n, dim3& grid, dim3& block)
{
	int sm_count = 0;
	cudaDeviceProp prop{};
	CUDA_OK(cudaGetDeviceProperties(&prop, 0));
	sm_count = prop.multiProcessorCount;

	const int blk = 256;
	size_t blocks = (n + blk - 1) / blk;
	blocks = std::max(blocks, (size_t)sm_count * 4);
	grid = dim3((unsigned)blocks);
	block = dim3(blk);
}

// ---- 64-bit benchmark ----
static void BM_add_mod_384_u64(benchmark::State& state)
{
	const size_t n = static_cast<size_t>(state.range(0));
	GpuBuf64 buf;

	// host data
	std::vector<uint64_t> h_a(6 * n);
	std::vector<uint64_t> h_b(6 * n);

	fill_host_u64(h_a);
	fill_host_u64(h_b);

	// device alloc + copy
	buf.n = n;
	CUDA_OK(cudaMalloc(&buf.d_a, 6 * n * sizeof(uint64_t)));
	CUDA_OK(cudaMalloc(&buf.d_b, 6 * n * sizeof(uint64_t)));
	CUDA_OK(cudaMalloc(&buf.d_r, 6 * n * sizeof(uint64_t)));

	CUDA_OK(cudaMemcpy(buf.d_a, h_a.data(), 6 * n * sizeof(uint64_t), cudaMemcpyHostToDevice));
	CUDA_OK(cudaMemcpy(buf.d_b, h_b.data(), 6 * n * sizeof(uint64_t), cudaMemcpyHostToDevice));

	dim3 grid, block;
	choose_launch(n, grid, block);

	// warmup
	kernel_add_mod_384_u64 << <grid, block >> > (buf.d_r, buf.d_a, buf.d_b, n);
	CUDA_OK(cudaDeviceSynchronize());

	cudaEvent_t evStart, evStop;
	CUDA_OK(cudaEventCreate(&evStart));
	CUDA_OK(cudaEventCreate(&evStop));

	double total_ms = 0.0;

	for (auto _ : state)
	{
		CUDA_OK(cudaEventRecord(evStart));
		kernel_add_mod_384_u64 << <grid, block >> > (buf.d_r, buf.d_a, buf.d_b, n);
		CUDA_OK(cudaEventRecord(evStop));
		CUDA_OK(cudaEventSynchronize(evStop));
		float ms = 0.0f;
		CUDA_OK(cudaEventElapsedTime(&ms, evStart, evStop));
		total_ms += ms;
		state.SetIterationTime(ms * 1e-3);
	}

	CUDA_OK(cudaEventDestroy(evStart));
	CUDA_OK(cudaEventDestroy(evStop));

	// bytes/second (device time)
	const double bytes_per_op = 144.0; // per element
	const double iters = static_cast<double>(state.iterations());
	const double total_s = total_ms * 1e-3;
	const double total_bytes = bytes_per_op * n * iters;

	const double gib_per_s = (total_bytes / total_s) / (1024.0 * 1024.0 * 1024.0);
	state.counters["GiB/s"] = gib_per_s;  // plain double works across GB versions

	state.SetItemsProcessed(static_cast<int64_t>(n) * state.iterations());

	CUDA_OK(cudaFree(buf.d_a));
	CUDA_OK(cudaFree(buf.d_b));
	CUDA_OK(cudaFree(buf.d_r));
}
BENCHMARK(BM_add_mod_384_u64)
->UseManualTime()
->Unit(benchmark::kNanosecond)
->Arg(1 << 10)   // 1K elems
->Arg(1 << 14)   // 16K
->Arg(1 << 18)   // 256K
->Arg(1 << 20);  // 1M


#ifdef BLST_ADD_MOD_384_U32_PTX
// ---- 32-bit benchmark ----
static void BM_add_mod_384_u32(benchmark::State& state)
{
	const size_t n = static_cast<size_t>(state.range(0));
	GpuBuf32 buf;

	// host data
	std::vector<uint32_t> h_a(12 * n);
	std::vector<uint32_t> h_b(12 * n);

	fill_host_u32(h_a);
	fill_host_u32(h_b);

	// device alloc + copy
	buf.n = n;
	CUDA_OK(cudaMalloc(&buf.d_a, 12 * n * sizeof(uint32_t)));
	CUDA_OK(cudaMalloc(&buf.d_b, 12 * n * sizeof(uint32_t)));
	CUDA_OK(cudaMalloc(&buf.d_r, 12 * n * sizeof(uint32_t)));

	CUDA_OK(cudaMemcpy(buf.d_a, h_a.data(), 12 * n * sizeof(uint32_t), cudaMemcpyHostToDevice));
	CUDA_OK(cudaMemcpy(buf.d_b, h_b.data(), 12 * n * sizeof(uint32_t), cudaMemcpyHostToDevice));

	dim3 grid, block;
	choose_launch(n, grid, block);

	// warmup
	kernel_add_mod_384_u32 << <grid, block >> > (buf.d_r, buf.d_a, buf.d_b, n);
	CUDA_OK(cudaDeviceSynchronize());

	cudaEvent_t evStart, evStop;
	CUDA_OK(cudaEventCreate(&evStart));
	CUDA_OK(cudaEventCreate(&evStop));

	double total_ms = 0.0;

	for (auto _ : state)
	{
		CUDA_OK(cudaEventRecord(evStart));
		kernel_add_mod_384_u32 << <grid, block >> > (buf.d_r, buf.d_a, buf.d_b, n);
		CUDA_OK(cudaEventRecord(evStop));
		CUDA_OK(cudaEventSynchronize(evStop));
		float ms = 0.0f;
		CUDA_OK(cudaEventElapsedTime(&ms, evStart, evStop));
		total_ms += ms;
		state.SetIterationTime(ms * 1e-3);
	}

	CUDA_OK(cudaEventDestroy(evStart));
	CUDA_OK(cudaEventDestroy(evStop));

	const double bytes_per_op = 144.0; // same as 64-bit layout
	const double iters = static_cast<double>(state.iterations());
	const double total_s = total_ms * 1e-3;
	const double total_bytes = bytes_per_op * n * iters;

	const double gib_per_s = (total_bytes / total_s) / (1024.0 * 1024.0 * 1024.0);
	state.counters["GiB/s"] = gib_per_s;  // plain double works across GB versions
	state.SetItemsProcessed(static_cast<int64_t>(n) * state.iterations());

	CUDA_OK(cudaFree(buf.d_a));
	CUDA_OK(cudaFree(buf.d_b));
	CUDA_OK(cudaFree(buf.d_r));
}
BENCHMARK(BM_add_mod_384_u32)
->UseManualTime()
->Unit(benchmark::kNanosecond)
->Arg(1 << 10)
->Arg(1 << 14)
->Arg(1 << 18)
->Arg(1 << 20);

#endif


// ---- c impl 384-bit benchmark ----
static void BM_add_mod_384_basic(benchmark::State& state)
{
	const size_t n = static_cast<size_t>(state.range(0));
	GpuBufFp buf;

	// host data
	std::vector<Fp> h_a(n);
	std::vector<Fp> h_b(n);

	fill_host_Fp(h_a);
	fill_host_Fp(h_b);


	// device alloc + copy
	buf.n = n;
	CUDA_OK(cudaMalloc(&buf.d_a, n * sizeof(Fp)));
	CUDA_OK(cudaMalloc(&buf.d_b, n * sizeof(Fp)));
	CUDA_OK(cudaMalloc(&buf.d_r, n * sizeof(Fp)));

	CUDA_OK(cudaMemcpy(buf.d_a, h_a.data(), n * sizeof(Fp), cudaMemcpyHostToDevice));
	CUDA_OK(cudaMemcpy(buf.d_b, h_b.data(), n * sizeof(Fp), cudaMemcpyHostToDevice));

	dim3 grid, block;
	choose_launch(n, grid, block);

	// warmup
	//kernel_add_mod_384_u32 << <grid, block >> > (buf.d_r, buf.d_a, buf.d_b, buf.d_p, n);
	add_fp_vec(buf.d_r, buf.d_a, buf.d_b, n);
	CUDA_OK(cudaDeviceSynchronize());

	cudaEvent_t evStart, evStop;
	CUDA_OK(cudaEventCreate(&evStart));
	CUDA_OK(cudaEventCreate(&evStop));

	double total_ms = 0.0;

	for (auto _ : state)
	{
		CUDA_OK(cudaEventRecord(evStart));
		//kernel_add_mod_384_u32 << <grid, block >> > (buf.d_r, buf.d_a, buf.d_b, buf.d_p, n);
		add_fp_vec(buf.d_r, buf.d_a, buf.d_b, n);
		CUDA_OK(cudaEventRecord(evStop));
		CUDA_OK(cudaEventSynchronize(evStop));
		float ms = 0.0f;
		CUDA_OK(cudaEventElapsedTime(&ms, evStart, evStop));
		total_ms += ms;
		state.SetIterationTime(ms * 1e-3);
	}

	CUDA_OK(cudaEventDestroy(evStart));
	CUDA_OK(cudaEventDestroy(evStop));

	const double bytes_per_op = 144.0; // same as 64-bit layout
	const double iters = static_cast<double>(state.iterations());
	const double total_s = total_ms * 1e-3;
	const double total_bytes = bytes_per_op * n * iters;

	const double gib_per_s = (total_bytes / total_s) / (1024.0 * 1024.0 * 1024.0);
	state.counters["GiB/s"] = gib_per_s;  // plain double works across GB versions
	state.SetItemsProcessed(static_cast<int64_t>(n) * state.iterations());

	CUDA_OK(cudaFree(buf.d_a));
	CUDA_OK(cudaFree(buf.d_b));
	CUDA_OK(cudaFree(buf.d_r));
}
BENCHMARK(BM_add_mod_384_basic)
->UseManualTime()
->Unit(benchmark::kNanosecond)
->Arg(1 << 10)
->Arg(1 << 14)
->Arg(1 << 18)
->Arg(1 << 20);












namespace temp {


	// 384-bit Montgomery multiply, base = 2^64, n = 6
	// n0 = -p^{-1} mod 2^64
	// No __int128. CUDA device code only (uses __umul64hi).
	struct wide64
	{
		uint64_t lo;
		uint64_t hi;
	};


	static BLST_DEVICE BLST_FORCEINLINE
		void mul_add_wide64(
			wide64& acc,
			const uint64_t a,
			const uint64_t b
		)
	{
		uint64_t prod_lo;
		uint64_t prod_hi;
		uint64_t lo_old;
		uint64_t sum_lo;
		uint64_t c1;

		prod_lo = a * b;
		prod_hi = __umul64hi(a, b);

		lo_old = acc.lo;
		sum_lo = lo_old + prod_lo;
		c1 = (sum_lo < lo_old) ? 1ull : 0ull;

		acc.lo = sum_lo;
		acc.hi = acc.hi + prod_hi + c1;
	}

	static BLST_DEVICE BLST_FORCEINLINE
		void add_limb_into_wide(
			wide64& acc,
			const uint64_t x
		)
	{
		uint64_t lo_old;
		uint64_t sum_lo;
		uint64_t c;

		lo_old = acc.lo;
		sum_lo = lo_old + x;
		c = (sum_lo < lo_old) ? 1ull : 0ull;

		acc.lo = sum_lo;
		acc.hi = acc.hi + c;
	}

	// If you want to use Peter's helper instead of add_limb_into_wide:
	// BLST_HD BLST_FORCEINLINE
	// void cc_add_wide_wide_limb(wide64* ret, const wide64& a, const uint64_t& b)
	// { wide64 aa = a; ret->lo = aa.lo + b; ret->hi = aa.hi + ((ret->lo < aa.lo) ? 1ull : 0ull); }
	__constant__ uint64_t p_[] = { BLST_P_LIMBS };

	//    BLST_DEVICE BLST_FORCEINLINE
	//        void ptx_mul_mont_384_u64(
	//            uint64_t* BLST_restrict r,
	//            const uint64_t* BLST_restrict a,
	//            const uint64_t* BLST_restrict b
	//            //const uint64_t* BLST_restrict p,
	//            //const uint64_t n0
	//        )
	//    {
	//		auto p = p_;
	//        const uint64_t n0 = BLST_p0;
	//        const int n = 6;
	//
	//        uint64_t t[7];
	//        int i;
	//        for (i = 0; i < n + 1; i++)
	//            t[i] = 0;
	//
	//        int j;
	//        for (j = 0; j < n; j++)
	//        {
	//            uint64_t bj;
	//            uint64_t carry;
	//            int k;
	//
	//            bj = b[j];
	//            carry = 0;
	//
	//            // t += a * b_j
	//#pragma unroll 6
	//            for (k = 0; k < n; k++)
	//            {
	//                wide64 acc;
	//                acc.lo = t[k];
	//                acc.hi = carry;
	//                mul_add_wide64(acc, a[k], bj);
	//                t[k] = acc.lo;
	//                carry = acc.hi;
	//            }
	//
	//            // propagate into t[n]
	//            {
	//                wide64 acc;
	//                acc.lo = t[n];
	//                acc.hi = 0;
	//                add_limb_into_wide(acc, carry);
	//                t[n] = acc.lo;
	//                // acc.hi is at most 1 here; keep it in t[n] by design (carry into next stage)
	//                // In CIOS we keep the extra word in t[n].
	//            }
	//
	//            // m = (t[0] * n0) mod 2^64
	//            uint64_t m;
	//            m = t[0] * n0;
	//
	//            // t += m * p, and shift (drop the least limb) as we store
	//            carry = 0;
	//
	//#pragma unroll 6
	//            for (k = 0; k < n; k++)
	//            {
	//                wide64 acc;
	//                acc.lo = t[k];
	//                acc.hi = carry;
	//                mul_add_wide64(acc, p[k], m);
	//
	//                if (k != 0)
	//                    t[k - 1] = acc.lo;
	//
	//                carry = acc.hi;
	//            }
	//
	//            {
	//                wide64 acc;
	//                acc.lo = t[n];
	//                acc.hi = 0;
	//                add_limb_into_wide(acc, carry);
	//                t[n - 1] = acc.lo;
	//                t[n] = acc.hi;
	//            }
	//        }
	//
	//        // Conditional subtract: diff = t - p
	//        uint64_t diff[6];
	//        uint64_t borrow;
	//        borrow = 0;
	//
	//        for (i = 0; i < n; i++)
	//        {
	//            uint64_t ti;
	//            uint64_t pi;
	//            uint64_t ri;
	//            uint64_t bi;
	//
	//            ti = t[i];
	//            pi = p[i];
	//
	//            ri = ti - pi;
	//            bi = (ti < pi) ? 1ull : 0ull;
	//
	//            if (borrow != 0)
	//            {
	//                uint64_t old = ri;
	//                ri = ri - 1ull;
	//                bi = bi | (old == 0ull ? 1ull : 0ull);
	//            }
	//
	//            diff[i] = ri;
	//            borrow = bi;
	//        }
	//
	//        // If borrow == 0 => t >= p => choose diff
	//        // If borrow == 1 => t < p  => choose t
	//        // mask = 0xFFFF.. if borrow == 0, else 0
	//        uint64_t mask;
	//        mask = 0ull - (uint64_t)(borrow ^ 1ull);
	//
	//        for (i = 0; i < n; i++)
	//        {
	//            uint64_t di;
	//            uint64_t ti;
	//            uint64_t ri;
	//
	//            di = diff[i];
	//            ti = t[i];
	//            ri = (di & mask) | (ti & ~mask);
	//            r[i] = ri;
	//        }
	//    }


	// ===== 384-bit Montgomery multiply (6x64 CIOS, no __int128) =====
	// r = a * b * R^{-1} mod p, base = 2^64, n = 6
	BLST_DEVICE BLST_FORCEINLINE
		void ptx_mul_mont_384_u64(
			uint64_t* BLST_restrict r,
			const uint64_t* BLST_restrict a,
			const uint64_t* BLST_restrict b
		)
	{
		const uint64_t* BLST_restrict p = p_;
		const uint64_t n0 = BLST_p0;
		const int n = 6;

		uint64_t t[7];
		int i;
		for (i = 0; i < n + 1; i++)
			t[i] = 0ull;

		int j;
		for (j = 0; j < n; j++)
		{
			uint64_t bj;
			uint64_t carry;
			bj = b[j];

			// ---- t += a * bj  (schoolbook column) ----
			carry = 0ull;
			int k;
#pragma unroll 6
			for (k = 0; k < n; k++)
			{
				uint64_t lo;
				uint64_t plo;
				uint64_t phi;
				uint64_t sum;
				uint64_t c1;

				uint64_t ak;
				ak = a[k];

				lo = t[k];
				plo = ak * bj;
				phi = __umul64hi(ak, bj);

				sum = lo + plo;
				c1 = (sum < lo) ? 1ull : 0ull;

				t[k] = sum;
				carry = carry + phi + c1;
			}

			// propagate carry into t[n]
			{
				uint64_t old;
				uint64_t sum;
				uint64_t c;

				old = t[n];
				sum = old + carry;
				c = (sum < old) ? 1ull : 0ull;

				t[n] = sum;
				// c would be an extra limb; CIOS keeps t as 7 limbs, so we drop c.
			}

			// ---- Montgomery step: m = t[0] * n0 (mod 2^64) ----
			uint64_t m;
			m = t[0] * n0;

			// ---- t += m * p, with the implicit right-shift (drop t[0]) ----
			carry = 0ull;
#pragma unroll 6
			for (k = 0; k < n; k++)
			{
				uint64_t lo;
				uint64_t plo;
				uint64_t phi;
				uint64_t sum;
				uint64_t c1;

				lo = t[k];
				plo = p[k] * m;
				phi = __umul64hi(p[k], m);

				sum = lo + plo;
				c1 = (sum < lo) ? 1ull : 0ull;

				if (k != 0)
					t[k - 1] = sum;

				carry = carry + phi + c1;
			}

			// write the most significant kept limb and update t[n]
			{
				uint64_t old;
				uint64_t sum;
				uint64_t c;

				old = t[n];
				sum = old + carry;
				c = (sum < old) ? 1ull : 0ull;

				t[n - 1] = sum;
				t[n] = c;
			}
		}

		// ---- Final conditional subtract: diff = t - p (branchless per-limb) ----
		uint64_t diff[6];
		uint64_t borrow;

		borrow = 0ull;
#pragma unroll 6
		for (i = 0; i < n; i++)
		{
			uint64_t ti;
			uint64_t pi;
			uint64_t spi;   // pi + borrow
			uint64_t ri;    // ti - (pi + borrow)

			ti = t[i];
			pi = p[i];
			spi = pi + borrow;
			ri = ti - spi;

			borrow = (ti < spi) ? 1ull : 0ull;
			diff[i] = ri;
		}

		// If borrow == 0 => t >= p => choose diff; else choose t
		// mask = 0xFFFF..FFFF when borrow == 0, else 0
		uint64_t mask;
		mask = 0ull - (borrow ^ 1ull);

#pragma unroll 6
		for (i = 0; i < n; i++)
		{
			uint64_t di;
			uint64_t ti;
			uint64_t ri;

			di = diff[i];
			ti = t[i];
			ri = (di & mask) | (ti & ~mask);
			r[i] = ri;
		}
	}

	// 384-bit Montgomery multiply, base = 2^32, n = 12
	// n0 must be -p^{-1} mod 2^32
	static BLST_DEVICE BLST_FORCEINLINE
		uint64_t ptx_madwide_u32(uint32_t a, uint32_t b, uint64_t add64)
	{
		uint64_t out;
		asm volatile(
			"{\n\t"
			" .reg .u64 t64;\n\t"
			" mad.wide.u32 t64, %1, %2, %3;\n\t"
			" mov.b64 %0, t64;\n\t"
			"}\n\t"
			: "=l"(out)
			: "r"(a), "r"(b), "l"(add64)
			);
		return out;
	}

	BLST_DEVICE BLST_FORCEINLINE
		void ptx_mul_mont_384_u32(
			uint32_t* BLST_restrict r,
			const uint32_t* BLST_restrict a,
			const uint32_t* BLST_restrict b
			//const uint32_t* BLST_restrict p,
			//uint32_t n0
		)
	{
		auto p = (uint32_t*)p_;
		const uint64_t n0 = BLST_p0;

		const int n = 12;
		uint32_t t[13];
		int i;
		int j;

		for (i = 0; i < n + 1; i++)
			t[i] = 0;

		for (j = 0; j < n; j++)
		{
			uint32_t bj;
			uint32_t carry32;
			uint64_t acc64;
			int k;

			bj = b[j];
			carry32 = 0;

#pragma unroll 12
			for (k = 0; k < n; k++)
			{
				acc64 = ptx_madwide_u32(a[k], bj, (uint64_t)t[k] + (uint64_t)carry32);
				t[k] = (uint32_t)acc64;
				carry32 = (uint32_t)(acc64 >> 32);
			}

			t[n] = t[n] + carry32;

			uint32_t m;
			m = (uint32_t)((uint64_t)t[0] * (uint64_t)n0);

			carry32 = 0;

#pragma unroll 12
			for (k = 0; k < n; k++)
			{
				acc64 = ptx_madwide_u32(p[k], m, (uint64_t)t[k] + (uint64_t)carry32);
				if (k != 0)
					t[k - 1] = (uint32_t)acc64;
				carry32 = (uint32_t)(acc64 >> 32);
			}

			acc64 = (uint64_t)t[n] + (uint64_t)carry32;
			t[n - 1] = (uint32_t)acc64;
			t[n] = (uint32_t)(acc64 >> 32);
		}

		uint32_t diff[12];
		uint32_t borrow;
		uint64_t sub64;

		borrow = 0;

		for (i = 0; i < n; i++)
		{
			sub64 = (uint64_t)t[i] - (uint64_t)p[i] - (uint64_t)borrow;
			diff[i] = (uint32_t)sub64;
			borrow = (uint32_t)((sub64 >> 63) & 1);
		}

		// Decide between t and diff:
		// In CIOS, choose t if final carry (t[n] != 0) > borrow; otherwise choose diff
		// Here, in base 2^32, t[n] is 0 or 1
		uint32_t tn;
		uint32_t gt;
		uint32_t mask;

		tn = t[n];
		gt = (tn > borrow) ? 1u : 0u;
		mask = 0u - gt;

		for (i = 0; i < n; i++)
			r[i] = (t[i] & mask) | (diff[i] & ~mask);
	}

}





// ==== __global__ kernels (1 element per thread, strided) ====
BLST_GLOBAL //__launch_bounds__(128, 2)
void kernel_mul_mont_384_u64
(
	uint64_t* BLST_restrict out,
	const uint64_t* BLST_restrict a,
	const uint64_t* BLST_restrict b,
	size_t n_elems
)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	size_t stride = gridDim.x * blockDim.x;
	for (size_t i = idx; i < n_elems; i += stride)
	{
		const uint64_t* pa = a + 6 * i;
		const uint64_t* pb = b + 6 * i;
		uint64_t* pr = out + 6 * i;
		Fp tt;
		//ptx_mul_mont_u64<Fp>(pr, pa, pb/*, (u64*)Fp::modulus().data, BLST_p0*/);
		ptx_mul_mont_u64<Fp>(tt.data, pa, pb);
		ptx_mul_mont_u64<Fp>(tt.data, tt.data, pb);
		ptx_mul_mont_u64<Fp>(tt.data, tt.data, pb);
		ptx_mul_mont_u64<Fp>(tt.data, tt.data, pb);
		ptx_mul_mont_u64<Fp>(tt.data, tt.data, pb);
		ptx_mul_mont_u64<Fp>(tt.data, tt.data, pb);
		ptx_mul_mont_u64<Fp>(pr, tt.data, pb);
	}
}

BLST_GLOBAL void kernel_mul_mont_384_cc
(
	uint32_t* BLST_restrict out,
	const uint32_t* BLST_restrict a,
	const uint32_t* BLST_restrict b,
	size_t n_elems
)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	size_t stride = gridDim.x * blockDim.x;
	for (size_t i = idx; i < n_elems; i += stride)
	{
		const Fp& pa = *(Fp*)(a + 12 * i);
		const Fp& pb = *(Fp*)(b + 12 * i);
		Fp& pr = *(Fp*)(out + 12 * i);
		Fp tt;
		
		//cc_mul_mont_n<6>(pr.data, pa.data, pb.data, Fp::modulus().data, BLST_p0);
		cc_mul_mont_n<6>(tt.data, pa.data, pb.data, Fp::modulus().data, BLST_p0);
		cc_mul_mont_n<6>(tt.data, tt.data, pb.data, Fp::modulus().data, BLST_p0);
		cc_mul_mont_n<6>(tt.data, tt.data, pb.data, Fp::modulus().data, BLST_p0);
		cc_mul_mont_n<6>(tt.data, tt.data, pb.data, Fp::modulus().data, BLST_p0);
		cc_mul_mont_n<6>(tt.data, tt.data, pb.data, Fp::modulus().data, BLST_p0);
		cc_mul_mont_n<6>(tt.data, tt.data, pb.data, Fp::modulus().data, BLST_p0);
		cc_mul_mont_n<6>(pr.data, tt.data, pb.data, Fp::modulus().data, BLST_p0);
	}
}

template<
	typename KernelT
>
static void choose_launch_occupancy(
	KernelT kernel,
	size_t n,
	int dynamic_smem_bytes,
	dim3& grid,
	dim3& block
)
{
	int minGridSize;
	int blockSize;
	cudaError_t st;

	minGridSize = 0;
	blockSize = 0;

	st = cudaOccupancyMaxPotentialBlockSize(
		&minGridSize,
		&blockSize,
		kernel,
		dynamic_smem_bytes,
		0
	);
	if (st != cudaSuccess)
	{
		// Fallback
		blockSize = 128;
	}

	size_t blocks64;
	unsigned blocks32;

	blocks64 = (n + blockSize - 1) / blockSize;
	if (blocks64 > std::numeric_limits<unsigned>::max())
		blocks32 = std::numeric_limits<unsigned>::max();
	else
		blocks32 = static_cast<unsigned>(blocks64);

	block = dim3(static_cast<unsigned>(blockSize));
	grid = dim3(blocks32);
}


// ---- 64-bit benchmark ----
static void BM_mul_mont_384_u64(benchmark::State& state)
{
	const size_t n = static_cast<size_t>(state.range(0));
	GpuBuf64 buf;

	// host data
	std::vector<uint64_t> h_a(6 * n);
	std::vector<uint64_t> h_b(6 * n);

	fill_host_u64(h_a);
	fill_host_u64(h_b);

	// device alloc + copy
	buf.n = n;
	CUDA_OK(cudaMalloc(&buf.d_a, 6 * n * sizeof(uint64_t)));
	CUDA_OK(cudaMalloc(&buf.d_b, 6 * n * sizeof(uint64_t)));
	CUDA_OK(cudaMalloc(&buf.d_r, 6 * n * sizeof(uint64_t)));

	CUDA_OK(cudaMemcpy(buf.d_a, h_a.data(), 6 * n * sizeof(uint64_t), cudaMemcpyHostToDevice));
	CUDA_OK(cudaMemcpy(buf.d_b, h_b.data(), 6 * n * sizeof(uint64_t), cudaMemcpyHostToDevice));

	dim3 grid, block;
	choose_launch(n, grid, block);

	// warmup
	kernel_mul_mont_384_u64 << <grid, block >> > (buf.d_r, buf.d_a, buf.d_b, n);
	CUDA_OK(cudaDeviceSynchronize());

	cudaEvent_t evStart, evStop;
	CUDA_OK(cudaEventCreate(&evStart));
	CUDA_OK(cudaEventCreate(&evStop));

	double total_ms = 0.0;

	for (auto _ : state)
	{
		CUDA_OK(cudaEventRecord(evStart));
		kernel_mul_mont_384_u64 << <grid, block >> > (buf.d_r, buf.d_a, buf.d_b, n);
		CUDA_OK(cudaEventRecord(evStop));
		CUDA_OK(cudaEventSynchronize(evStop));
		float ms = 0.0f;
		CUDA_OK(cudaEventElapsedTime(&ms, evStart, evStop));
		total_ms += ms;
		state.SetIterationTime(ms * 1e-3);
	}

	CUDA_OK(cudaEventDestroy(evStart));
	CUDA_OK(cudaEventDestroy(evStop));

	// bytes/second (device time)
	const double bytes_per_op = 144.0; // per element
	const double iters = static_cast<double>(state.iterations());
	const double total_s = total_ms * 1e-3;
	const double total_bytes = bytes_per_op * n * iters;

	const double gib_per_s = (total_bytes / total_s) / (1024.0 * 1024.0 * 1024.0);
	state.counters["GiB/s"] = gib_per_s;  // plain double works across GB versions

	state.SetItemsProcessed(static_cast<int64_t>(n) * state.iterations());

	CUDA_OK(cudaFree(buf.d_a));
	CUDA_OK(cudaFree(buf.d_b));
	CUDA_OK(cudaFree(buf.d_r));
}
BENCHMARK(BM_mul_mont_384_u64)
->UseManualTime()
->Unit(benchmark::kNanosecond)
->Arg(1 << 10)   // 1K elems
->Arg(1 << 14)   // 16K
->Arg(1 << 18)   // 256K
->Arg(1 << 20);  // 1M


// ---- 32-bit benchmark ----
static void BM_mul_mont_384_cc(benchmark::State& state)
{
	const size_t n = static_cast<size_t>(state.range(0));
	GpuBuf32 buf;

	// host data
	std::vector<uint32_t> h_a(12 * n);
	std::vector<uint32_t> h_b(12 * n);

	fill_host_u32(h_a);
	fill_host_u32(h_b);

	// device alloc + copy
	buf.n = n;
	CUDA_OK(cudaMalloc(&buf.d_a, 12 * n * sizeof(uint32_t)));
	CUDA_OK(cudaMalloc(&buf.d_b, 12 * n * sizeof(uint32_t)));
	CUDA_OK(cudaMalloc(&buf.d_r, 12 * n * sizeof(uint32_t)));

	CUDA_OK(cudaMemcpy(buf.d_a, h_a.data(), 12 * n * sizeof(uint32_t), cudaMemcpyHostToDevice));
	CUDA_OK(cudaMemcpy(buf.d_b, h_b.data(), 12 * n * sizeof(uint32_t), cudaMemcpyHostToDevice));

	dim3 grid, block;
	choose_launch(n, grid, block);

	// warmup
	kernel_mul_mont_384_cc << <grid, block >> > (buf.d_r, buf.d_a, buf.d_b, n);
	CUDA_OK(cudaDeviceSynchronize());

	cudaEvent_t evStart, evStop;
	CUDA_OK(cudaEventCreate(&evStart));
	CUDA_OK(cudaEventCreate(&evStop));

	double total_ms = 0.0;

	for (auto _ : state)
	{
		CUDA_OK(cudaEventRecord(evStart));
		kernel_mul_mont_384_cc << <grid, block >> > (buf.d_r, buf.d_a, buf.d_b, n);
		CUDA_OK(cudaEventRecord(evStop));
		CUDA_OK(cudaEventSynchronize(evStop));
		float ms = 0.0f;
		CUDA_OK(cudaEventElapsedTime(&ms, evStart, evStop));
		total_ms += ms;
		state.SetIterationTime(ms * 1e-3);
	}

	CUDA_OK(cudaEventDestroy(evStart));
	CUDA_OK(cudaEventDestroy(evStop));

	const double bytes_per_op = 144.0; // same as 64-bit layout
	const double iters = static_cast<double>(state.iterations());
	const double total_s = total_ms * 1e-3;
	const double total_bytes = bytes_per_op * n * iters;

	const double gib_per_s = (total_bytes / total_s) / (1024.0 * 1024.0 * 1024.0);
	state.counters["GiB/s"] = gib_per_s;  // plain double works across GB versions
	state.SetItemsProcessed(static_cast<int64_t>(n) * state.iterations());

	CUDA_OK(cudaFree(buf.d_a));
	CUDA_OK(cudaFree(buf.d_b));
	CUDA_OK(cudaFree(buf.d_r));
}
BENCHMARK(BM_mul_mont_384_cc)
->UseManualTime()
->Unit(benchmark::kNanosecond)
->Arg(1 << 10)
->Arg(1 << 14)
->Arg(1 << 18)
->Arg(1 << 20);





// ==== __global__ kernels (1 element per thread, strided) ====
BLST_GLOBAL //__launch_bounds__(128, 2)
void kernel_sqr_mont_384_u64
(
	Fp* BLST_restrict out,
	const Fp* BLST_restrict a,
	size_t n_elems
)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	size_t stride = gridDim.x * blockDim.x;
	for (size_t i = idx; i < n_elems; i += stride)
	{
		const Fp& pa = a[i];
		Fp& pr = out[i];
		Fp tt;
		//ptx_sqr_mont_u64<Fp>(pr, pa, pb/*, (u64*)Fp::modulus().data, BLST_p0*/);
		ptx_sqr_mont_u64<Fp>(tt.data, pa.data);
		ptx_sqr_mont_u64<Fp>(tt.data, tt.data);
		ptx_sqr_mont_u64<Fp>(tt.data, tt.data);
		ptx_sqr_mont_u64<Fp>(tt.data, tt.data);
		ptx_sqr_mont_u64<Fp>(tt.data, tt.data);
		ptx_sqr_mont_u64<Fp>(tt.data, tt.data);
		ptx_sqr_mont_u64<Fp>(pr.data, tt.data);
	}
}


// ---- c impl 384-bit benchmark ----
static void BM_sqr_mont_384_u64(benchmark::State& state)
{
	const size_t n = static_cast<size_t>(state.range(0));
	GpuBufFp buf;

	// host data
	std::vector<Fp> h_a(n);

	fill_host_Fp(h_a);


	// device alloc + copy
	buf.n = n;
	CUDA_OK(cudaMalloc(&buf.d_a, n * sizeof(Fp)));
	CUDA_OK(cudaMalloc(&buf.d_r, n * sizeof(Fp)));

	CUDA_OK(cudaMemcpy(buf.d_a, h_a.data(), n * sizeof(Fp), cudaMemcpyHostToDevice));

	dim3 grid, block;
	choose_launch(n, grid, block);

	// warmup
	kernel_sqr_mont_384_u64 << <grid, block >> > (buf.d_r, buf.d_a, n);
	CUDA_OK(cudaDeviceSynchronize());

	cudaEvent_t evStart, evStop;
	CUDA_OK(cudaEventCreate(&evStart));
	CUDA_OK(cudaEventCreate(&evStop));

	double total_ms = 0.0;

	for (auto _ : state)
	{
		CUDA_OK(cudaEventRecord(evStart));
		kernel_sqr_mont_384_u64 << <grid, block >> > (buf.d_r, buf.d_a, n);
		CUDA_OK(cudaEventRecord(evStop));
		CUDA_OK(cudaEventSynchronize(evStop));
		float ms = 0.0f;
		CUDA_OK(cudaEventElapsedTime(&ms, evStart, evStop));
		total_ms += ms;
		state.SetIterationTime(ms * 1e-3);
	}

	CUDA_OK(cudaEventDestroy(evStart));
	CUDA_OK(cudaEventDestroy(evStop));

	const double bytes_per_op = 144.0; // same as 64-bit layout
	const double iters = static_cast<double>(state.iterations());
	const double total_s = total_ms * 1e-3;
	const double total_bytes = bytes_per_op * n * iters;

	const double gib_per_s = (total_bytes / total_s) / (1024.0 * 1024.0 * 1024.0);
	state.counters["GiB/s"] = gib_per_s;  // plain double works across GB versions
	state.SetItemsProcessed(static_cast<int64_t>(n) * state.iterations());

	CUDA_OK(cudaFree(buf.d_a));
	CUDA_OK(cudaFree(buf.d_r));
}
BENCHMARK(BM_sqr_mont_384_u64)
->UseManualTime()
->Unit(benchmark::kNanosecond)
->Arg(1 << 10)
->Arg(1 << 14)
->Arg(1 << 18)
->Arg(1 << 20);



//void FpMain()
//{
//	const size_t n = 1 << 18;
//	GpuBuf64 buf;
//
//	// host data
//	std::vector<uint64_t> h_a(6 * n);
//	std::vector<uint64_t> h_b(6 * n);
//
//	fill_host_u64(h_a);
//	fill_host_u64(h_b);
//
//	// device alloc + copy
//	buf.n = n;
//	CUDA_OK(cudaMalloc(&buf.d_a, 6 * n * sizeof(uint64_t)));
//	CUDA_OK(cudaMalloc(&buf.d_b, 6 * n * sizeof(uint64_t)));
//	CUDA_OK(cudaMalloc(&buf.d_r, 6 * n * sizeof(uint64_t)));
//
//	CUDA_OK(cudaMemcpy(buf.d_a, h_a.data(), 6 * n * sizeof(uint64_t), cudaMemcpyHostToDevice));
//	CUDA_OK(cudaMemcpy(buf.d_b, h_b.data(), 6 * n * sizeof(uint64_t), cudaMemcpyHostToDevice));
//
//	dim3 grid, block;
//	choose_launch(n, grid, block);
//
//	for (u64 i = 0; i < 100; ++i)
//	{
//		kernel_mul_mont_384_u64 << <grid, block >> > (buf.d_r, buf.d_a, buf.d_b, n);
//	}
//
//}








































// ==== __global__ kernels (1 element per thread, strided) ====
BLST_GLOBAL //__launch_bounds__(128, 2)
void kernel_mul_mont_256_u64
(
	uint64_t* BLST_restrict out,
	const uint64_t* BLST_restrict a,
	const uint64_t* BLST_restrict b,
	size_t n_elems
)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	size_t stride = gridDim.x * blockDim.x;
	for (size_t i = idx; i < n_elems; i += stride)
	{
		const uint64_t* pa = a + 4 * i;
		const uint64_t* pb = b + 4 * i;
		uint64_t* pr = out + 4 * i;
		ptx_mul_mont_u64<Fr>(pr, pa, pb/*, (u64*)Fp::modulus().data, BLST_p0*/);
	}
}


// ---- 64-bit benchmark ----
static void BM_mul_mont_256_u64(benchmark::State& state)
{
	const size_t n = static_cast<size_t>(state.range(0));
	GpuBuf64 buf;

	static_assert(sizeof(Fr) / sizeof(uint64_t) == 4, "Fr size mismatch");
	// host data
	std::vector<uint64_t> h_a(4 * n);
	std::vector<uint64_t> h_b(4 * n);

	fill_host_u64(h_a);
	fill_host_u64(h_b);

	// device alloc + copy
	buf.n = n;
	CUDA_OK(cudaMalloc(&buf.d_a, 4 * n * sizeof(uint64_t)));
	CUDA_OK(cudaMalloc(&buf.d_b, 4 * n * sizeof(uint64_t)));
	CUDA_OK(cudaMalloc(&buf.d_r, 4 * n * sizeof(uint64_t)));

	CUDA_OK(cudaMemcpy(buf.d_a, h_a.data(), 4 * n * sizeof(uint64_t), cudaMemcpyHostToDevice));
	CUDA_OK(cudaMemcpy(buf.d_b, h_b.data(), 4 * n * sizeof(uint64_t), cudaMemcpyHostToDevice));

	dim3 grid, block;
	choose_launch(n, grid, block);

	// warmup
	kernel_mul_mont_256_u64 << <grid, block >> > (buf.d_r, buf.d_a, buf.d_b, n);
	CUDA_OK(cudaDeviceSynchronize());

	cudaEvent_t evStart, evStop;
	CUDA_OK(cudaEventCreate(&evStart));
	CUDA_OK(cudaEventCreate(&evStop));

	double total_ms = 0.0;

	for (auto _ : state)
	{
		CUDA_OK(cudaEventRecord(evStart));
		kernel_mul_mont_256_u64 << <grid, block >> > (buf.d_r, buf.d_a, buf.d_b, n);
		CUDA_OK(cudaEventRecord(evStop));
		CUDA_OK(cudaEventSynchronize(evStop));
		float ms = 0.0f;
		CUDA_OK(cudaEventElapsedTime(&ms, evStart, evStop));
		total_ms += ms;
		state.SetIterationTime(ms * 1e-3);
	}

	CUDA_OK(cudaEventDestroy(evStart));
	CUDA_OK(cudaEventDestroy(evStop));

	// bytes/second (device time)
	const double bytes_per_op = 144.0; // per element
	const double iters = static_cast<double>(state.iterations());
	const double total_s = total_ms * 1e-3;
	const double total_bytes = bytes_per_op * n * iters;

	const double gib_per_s = (total_bytes / total_s) / (1024.0 * 1024.0 * 1024.0);
	state.counters["GiB/s"] = gib_per_s;  // plain double works across GB versions

	state.SetItemsProcessed(static_cast<int64_t>(n) * state.iterations());

	CUDA_OK(cudaFree(buf.d_a));
	CUDA_OK(cudaFree(buf.d_b));
	CUDA_OK(cudaFree(buf.d_r));
}
BENCHMARK(BM_mul_mont_256_u64)
->UseManualTime()
->Unit(benchmark::kNanosecond)
->Arg(1 << 10)   // 1K elems
->Arg(1 << 14)   // 16K
->Arg(1 << 18)   // 256K
->Arg(1 << 20);  // 1M