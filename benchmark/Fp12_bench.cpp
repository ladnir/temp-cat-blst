#include <benchmark/benchmark.h>
#include <random>
#include <cstring>
#include <stdint.h>
#include <thrust/device_vector.h>

// Include the actual BLST field headers
#include "blst/fields/Fp12.h"
#include "blst/curve/point.h"
#include "blst/curve/GT.h"

using namespace blst;

// Helper function to generate random Fp12 elements
void generate_random_fp12(Fp12& fp12) {
	static std::random_device rd;
	static std::mt19937_64 gen(rd());
	static std::uniform_int_distribution<uint64_t> dis;

	limb_t* p = (limb_t*)&fp12;
	for (size_t i = 0; i < sizeof(Fp12) / sizeof(limb_t); ++i) {
		p[i] = dis(gen);
	}

	// Ensure the generated element is reduced modulo the field characteristic
	for (auto& comp : fp12.data) {
		for (auto& part : comp.data) {
			for (auto& fp : part.data) {
				while (Fp::is_redc(fp) == 0) {
					Fp::rshift(fp, fp, 1);
				}
			}
		}
	}
}

// ========== MULTIPLICATION BENCHMARKS ==========

// Benchmark for Fp12 multiplication
static void BM_Fp12Multiplication(benchmark::State& state) {
	Fp12 a, b, c;
	generate_random_fp12(a);
	generate_random_fp12(b);

	for (auto _ : state) {
		Fp12::mul(c, a, b);
		benchmark::DoNotOptimize(c);
	}
	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12Multiplication);

// ========== SQUARING BENCHMARKS ==========

// Benchmark for Fp12 squaring
static void BM_Fp12Squaring(benchmark::State& state) {
	Fp12 a, c;
	generate_random_fp12(a);

	for (auto _ : state) {
		Fp12::sqr(c, a);
		benchmark::DoNotOptimize(c);
	}
	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12Squaring);

// Benchmark for Fp12 cyclotomic squaring
static void BM_Fp12CyclotomicSquaring(benchmark::State& state) {
	Fp12 a, c;
	generate_random_fp12(a);

	for (auto _ : state) {
		Fp12::cyclotomic_sqr(c, a);
		benchmark::DoNotOptimize(c);
	}
	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12CyclotomicSquaring);

// ========== INVERSION BENCHMARK ==========

// Benchmark for Fp12 inversion
static void BM_Fp12Inversion(benchmark::State& state) {
	Fp12 a, c;
	generate_random_fp12(a);

	for (auto _ : state) {
		Fp12::inverse(c, a);
		benchmark::DoNotOptimize(c);
	}
	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12Inversion);

// ========== CONJUGATION BENCHMARK ==========

// Benchmark for Fp12 conjugation
static void BM_Fp12Conjugate(benchmark::State& state) {
	Fp12 a;
	generate_random_fp12(a);

	for (auto _ : state) {
		// Conjugation is in-place
		Fp12::conjugate(a);
		benchmark::DoNotOptimize(a);
	}
	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12Conjugate);

// ========== FROBENIUS MAP BENCHMARK ==========

// Benchmark for Fp12 Frobenius map
static void BM_Fp12Frobenius(benchmark::State& state) {
	Fp12 a, c;
	generate_random_fp12(a);

	for (auto _ : state) {
		Fp12::frobenius_map(c, a, 1);
		benchmark::DoNotOptimize(c);
	}
	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12Frobenius);

// ========== EXPONENTIATION BENCHMARK ==========

// Benchmark for Fp12 exponentiation
static void BM_Fp12Exponentiation(benchmark::State& state) {
	Fp12 a, c;
	generate_random_fp12(a);

	// A small, simple exponent for benchmarking purposes
	byte exponent[32];

	std::memset(exponent, 128 | 0xcc, 32); // Example exponent

	for (auto _ : state) {
		Fp12::exp(c, a, exponent, 256);
		benchmark::DoNotOptimize(c);
	}
	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12Exponentiation);


// Helper: Generate a random scalar (256-bit)
void generate_random_scalar(byte scalar[33]) {
	static std::random_device rd;
	static std::mt19937_64 gen(rd());
	static std::uniform_int_distribution<uint64_t> dis;

	uint64_t* p = (uint64_t*)scalar;
	for (int i = 0; i < 4; ++i) {
		p[i] = dis(gen);
	}
}

// ... existing benchmarks ...

static void BM_Fp12CyclotomicExponentiation(benchmark::State& state) {
	Fp12 a, base, c;
	generate_random_fp12(a);
	Fp12::final_exp(base, a); // Ensure base is in cyclotomic subgroup

	byte scalar[33] = { 0 };
	generate_random_scalar(scalar);

	for (auto _ : state) {
		Fp12::cyclotomic_exp(c, base, scalar, 257);
		benchmark::DoNotOptimize(c);
	}
	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12CyclotomicExponentiation);

static void BM_Fp12ExponentiationGLV(benchmark::State& state) {
	Fp12 a, base, c;
	generate_random_fp12(a);
	Fp12::final_exp(base, a); // Ensure base is in cyclotomic subgroup

	byte scalar[32];
	generate_random_scalar(scalar);

	for (auto _ : state) {
		Fp12::cyclotomic_exp_glv(c, base, scalar);
		benchmark::DoNotOptimize(c);
	}
	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12ExponentiationGLV);

// We need a helper to split scalars for the pre-split benchmark.
// Since the library function might be internal, we duplicate the call logic
// or use the Fr helper if available. 
// Assuming Fr.h is available and has the constants/helpers we need.
// If gt_split_scalar_phi2_u256_barrett is not public, we can't call it directly.
// HOWEVER, we can just benchmark exp_glv with "random" 128-bit k1/k2 
// which is effectively what the split produces performance-wise.

static void BM_Fp12ExponentiationGLV_PreSplit(benchmark::State& state) {
	Fp12 a, base, c;
	generate_random_fp12(a);
	Fp12::final_exp(base, a);

	// Generate random 128-bit splits (simulation of pre-computation)
	// k1, k2 are approx 128 bits for 256-bit scalars.
	std::vector<std::pair<std::array<limb_t, 4>, std::array<limb_t, 4>>> splits(128);

	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_int_distribution<uint64_t> dis;

	for (auto& s : splits) {
		// Fill 128 bits (2 limbs) randomly, zero upper limbs
		s.first[0] = dis(gen); s.first[1] = dis(gen); s.first[2] = 0; s.first[3] = 0;
		s.second[0] = dis(gen); s.second[1] = dis(gen); s.second[2] = 0; s.second[3] = 0;
	}

	size_t idx = 0;
	const auto& s = splits[idx & 127];
	for (auto _ : state) {
		Fp12::cyclotomic_exp_glv(c, base, s.first.data(), s.second.data());
		benchmark::DoNotOptimize(c);
		idx++;
	}
	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12ExponentiationGLV_PreSplit);


//// We need a helper to split scalars for the pre-split benchmark.
//// Since the library function might be internal, we duplicate the call logic
//// or use the Fr helper if available. 
//// Assuming Fr.h is available and has the constants/helpers we need.
//// If gt_split_scalar_phi2_u256_barrett is not public, we can't call it directly.
//// HOWEVER, we can just benchmark exp_glv with "random" 128-bit k1/k2 
//// which is effectively what the split produces performance-wise.
//
//static void BM_Fp12ExponentiationGLV_opt(benchmark::State& state) {
//    Fp12 a, base, c;
//    generate_random_fp12(a);
//    Fp12::final_exp(base, a);
//
//    // Generate random 128-bit splits (simulation of pre-computation)
//    // k1, k2 are approx 128 bits for 256-bit scalars.
//    std::vector<std::pair<std::array<limb_t, 4>, std::array<limb_t, 4>>> splits(128);
//
//    std::random_device rd;
//    std::mt19937_64 gen(rd());
//    std::uniform_int_distribution<uint64_t> dis;
//
//    for (auto& s : splits) {
//        // Fill 128 bits (2 limbs) randomly, zero upper limbs
//        s.first[0] = dis(gen); s.first[1] = dis(gen); s.first[2] = 0; s.first[3] = 0;
//        s.second[0] = dis(gen); s.second[1] = dis(gen); s.second[2] = 0; s.second[3] = 0;
//    }
//
//    size_t idx = 0;
//    const auto& s = splits[idx & 127];
//    for (auto _ : state) {
//        Fp12::cyclotomic_exp_glv_optimized(c, base, s.first.data(), s.second.data());
//        benchmark::DoNotOptimize(c);
//        idx++;
//    }
//    state.SetItemsProcessed(state.iterations());
//}
//BENCHMARK(BM_Fp12ExponentiationGLV_opt);
//namespace blst {
//    BLST_HD int recode_glv_jsf(int8_t jsf[][2], const limb_t* k1, const limb_t* k2);
//
//    BLST_HD void cyclotomic_exp_glv_opt(Fp12& out, const Fp12& base, int8_t jsf[132][2], int len);
//}

// We need a helper to split scalars for the pre-split benchmark.
// Since the library function might be internal, we duplicate the call logic
// or use the Fr helper if available. 
// Assuming Fr.h is available and has the constants/helpers we need.
// If gt_split_scalar_phi2_u256_barrett is not public, we can't call it directly.
// HOWEVER, we can just benchmark exp_glv with "random" 128-bit k1/k2 
// which is effectively what the split produces performance-wise.

static void BM_Fp12ExponentiationGLV_opt(benchmark::State& state) {
	Fp12 a, base, c;
	generate_random_fp12(a);
	Fp12::final_exp(base, a);

	// Generate random 128-bit splits (simulation of pre-computation)
	// k1, k2 are approx 128 bits for 256-bit scalars.
	//std::vector<std::pair<std::array<limb_t, 4>, std::array<limb_t, 4>>> splits(128);

	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_int_distribution<uint64_t> dis;
	std::pair<std::array<limb_t, 4>, std::array<limb_t, 4>> s;

	//for (auto& s : splits) 
	{
		// Fill 128 bits (2 limbs) randomly, zero upper limbs
		s.first[0] = dis(gen); s.first[1] = dis(gen); s.first[2] = 0; s.first[3] = 0;
		s.second[0] = dis(gen); s.second[1] = dis(gen); s.second[2] = 0; s.second[3] = 0;
	}


	auto jsf = Fp12::glv_recode_jsf(s.first.data(), s.second.data());

	size_t idx = 0;
	//const auto& s = splits[idx & 127];
	for (auto _ : state) {
		//Fp12::cyclotomic_exp_glv_optimized(c, base, s.first.data(), s.second.data());
		Fp12::cyclotomic_exp_glv_jsf(c, base, jsf);
		benchmark::DoNotOptimize(c);
		idx++;
	}
	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12ExponentiationGLV_opt);



// 1. Basic GLV-4 Benchmark (Runtime Split + Exp)
static void BM_Fp12ExponentiationGLV4(benchmark::State& state) {
	Fp12 a, base, c;
	generate_random_fp12(a);
	Fp12::final_exp(base, a); // Ensure base is in cyclotomic subgroup

	byte scalar[32];
	generate_random_scalar(scalar);

	// Convert scalar bytes to limbs for the API
	limb_t k[4];
	std::memcpy(k, scalar, 32);

	for (auto _ : state) {
		// Calls the full GLV-4 logic including decomposition
		Fp12::cyclotomic_exp_glv4(c, base, k);
		benchmark::DoNotOptimize(c);
	}
	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12ExponentiationGLV4);

// 2. Constant Scalar GLV-4 Benchmark (Offline Split, Runtime Exp)
static void BM_Fp12ExponentiationGLV4_Constant(benchmark::State& state) {
	Fp12 a, base, c;
	generate_random_fp12(a);
	Fp12::final_exp(base, a);

	byte scalar[32];
	generate_random_scalar(scalar);
	limb_t k[4];
	std::memcpy(k, scalar, 32);

	// --- Offline Compilation Phase (Not Timed) ---
	// Compile the scalar into properties (Schedule + Sign Mask)
	GLV4_Exp props = Fp12::glv4_compile(k);

	// --- Runtime Benchmark Phase ---
	for (auto _ : state) {
		// Only run the execution phase, simulating fixed scalar usage
		Fp12::cyclotomic_exp_glv4(c, base, props);
		benchmark::DoNotOptimize(c);
	}
	state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12ExponentiationGLV4_Constant);




BLST_HD void print(const Fp12& a) {
	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 2; j++) {
			printf("Fp12[%d][%d]: ", i, j);
			for (int k = 0; k < 6; k++) {
				printf("%016llx ", (unsigned long long)a.data[i / 6].data[i % 6].data[j].data[k]);
			}
			printf("\n");
		}
	}
}

// Number of parallel instances (Threads * Blocks)
const int BATCH_SIZE = 128 * 2;

#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA Error: %s in %s at line %d\n", \
                   cudaGetErrorString(err), __FILE__, __LINE__); \
            exit(1); \
        } \
    } while (0)

// --- Kernels ---

// Wrapper for Basic Exp
__global__ void kernel_cyclotomic_exp_basic(Fp12* out, const Fp12* base, const byte* scalars) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= BATCH_SIZE) return;

	// Each thread reads a 32-byte scalar
	const byte* thread_scalar = scalars + (tid * 32);

	// Note: This relies on the implementation not using stack-heavy arrays internally
	// if generic stack is small. 'basic' is usually safe.
	Fp12::cyclotomic_exp_basic(out[tid], base[tid], thread_scalar, 256);
}
//// Struct for pair return to avoid std::pair dependency on device
//struct OpSelect { const Fp12* ptr; bool inv; };
//
//BLST_HD auto get_op(int u, int v, auto& base, auto& P2, auto& P1_P2, auto& P1_nP2) {
//	const Fp12* op = nullptr;
//	bool inv = false;
//	if (u == 0) { op = &P2; if (v < 0) inv = true; }
//	else if (v == 0) { op = &base; if (u < 0) inv = true; }
//	else if (u == v) { op = &P1_P2; if (u < 0) inv = true; }
//	else { op = &P1_nP2; if (u < 0) inv = true; }
//	return OpSelect{ op, inv };
//};
//
///// @brief Precomputed schedule for constant-time GLV-2 JSF exponentiation.
///// Contains a compact bytecode sequence of operations (squares/multiplies).
///// The ecnoding has two signed digits per index. Then we pack two such
///// pairs into a byte for storage efficiency.
//struct GLV_JSF_Exp_ {
//	uint8_t sched_len;     ///< Length of the valid schedule.
//	uint8_t schedule[71];  ///< Bytecode for the exponentiation loop.
//
//	// sets the idx-th digit pair to (d1, d2)
//	BLST_HD void push(int8_t d1, int8_t d2) {
//
//		auto idx = sched_len++;
//		BLST_ASSERT(idx < sizeof(schedule) * 2);
//
//		auto enc = ((d1 + 1) | ((d2 + 1) << 2));
//		enc <<= ((idx % 2) * 4);
//		enc |= schedule[idx / 2] & (~(0xF << ((idx % 2) * 4)));
//		schedule[idx / 2] = enc;
//	}
//
//	struct Pair {
//		int8_t d1;
//		int8_t d2;
//	};
//
//	// gets the idx-th digit pair (d1, d2)
//	BLST_HD Pair get(size_t idx) const {
//		if (idx >= sched_len)
//			printf("out of bounds 1\n");
//		if (idx / 2 >= sizeof(schedule))
//			printf("out of bounds 2\n");
//		auto enc = (schedule[idx / 2] >> ((idx % 2) * 4)) & 0xF;
//		auto d1 = (int8_t)((enc & 0x3) - 1);
//		auto d2 = (int8_t)(((enc >> 2) & 0x3) - 1);
//		return { d1, d2 };
//	}
//};
//
//
//// Helper for safe bit access
//BLST_HD inline unsigned get_scalar_bit(const limb_t* k, int i) {
//	if (i < 0) return 0;
//	// Assumes 64-bit limbs.
//	return (unsigned)((k[i / 64] >> (i % 64)) & 1);
//}
//
//BLST_HD GLV_JSF_Exp_ glv_recode_jsf(const limb_t* k1, const limb_t* k2)
//{
//	GLV_JSF_Exp_ result = { 0 }; // Zero-initialize for safe bitwise OR in set()
//	unsigned d1 = 0, d2 = 0; // Carries for k1, k2
//	unsigned u, v;
//	int i = 0;
//	result.sched_len = 0;
//
//	// Map scalars to signed digits {-1, 0, 1}
//	// Loop limit 132 covers 129-bit input + worst-case carry propagation
//	while (i < 133) {
//		unsigned k1_bit = get_scalar_bit(k1, i);
//		unsigned k2_bit = get_scalar_bit(k2, i);
//
//		// Lookahead bits for JSF decision
//		unsigned k1_next = get_scalar_bit(k1, i + 1);
//		unsigned k2_next = get_scalar_bit(k2, i + 1);
//
//		u = k1_bit + d1;
//		v = k2_bit + d2;
//
//		int j1, j2;
//
//		// JSF Logic for row 1
//		if ((u & 1) == 0) {
//			j1 = 0;
//		}
//		else {
//			// Cohen's JSF Logic:
//			// If u == 3 (mod 4) OR (u == 1 (mod 4) AND v == 2 (mod 4)), choose -1
//			// Note: v % 4 checks lookahead of v.
//			unsigned u_mod_4 = (u + (k1_next << 1)) & 3;
//			unsigned v_mod_4 = (v + (k2_next << 1)) & 3;
//
//			if (u_mod_4 == 3 || (u_mod_4 == 1 && v_mod_4 == 2)) j1 = -1;
//			else j1 = 1;
//		}
//
//		// JSF Logic for row 2
//		if ((v & 1) == 0) {
//			j2 = 0;
//		}
//		else {
//			unsigned u_mod_4 = (u + (k1_next << 1)) & 3;
//			unsigned v_mod_4 = (v + (k2_next << 1)) & 3;
//
//			if (v_mod_4 == 3 || (v_mod_4 == 1 && u_mod_4 == 2)) j2 = -1;
//			else j2 = 1;
//		}
//
//		// Bounds check for safety before writing
//		if (result.sched_len < 132) {
//			result.push(j1, j2);
//		}
//
//		// Calculate next carries: d = (u - j) / 2
//		d1 = (u - j1) >> 1;
//		d2 = (v - j2) >> 1;
//
//		i++;
//
//		// Quit early if we are past bit 130 and have no carries left
//		if (i >= 130 && d1 == 0 && d2 == 0) break;
//	}
//
//	// Trim leading zeros to find real MSB
//	int len = result.sched_len;
//	while (len > 0 && result.get(len - 1).d1 == 0 && result.get(len - 1).d2 == 0) {
//		len--;
//	}
//	result.sched_len = (uint8_t)len;
//
//	return result;
//}
//
//
//BLST_HD void cyclotomic_exp_glv_jsf(Fp12& out, const Fp12& base, const GLV_JSF_Exp_& jsf)
//{
//
//
//#if defined(__CUDA_ARCH__)
//	if (threadIdx.x == 0 && blockIdx.x == 0) {
//		printf("GPU Trace enter: base[0]=%llx\n", (unsigned long long)base[0][0][0][0]);
//	}
//#endif
//	int len = jsf.sched_len;
//	if (len == 0) {
//		vec_copy(&out, &Fp12::one(), sizeof(Fp12));
//#if defined(__CUDA_ARCH__)
//		if (threadIdx.x == 0 && blockIdx.x == 0) {
//			printf("GPU Trace early exit\n");
//		}
//#endif
//		return;
//	}
//
//	// 2. Precompute Lookup Points (Stack)
//	// We compute mixed points to reduce the loop to 1 add per step max.
//	Fp12 P2, P1_P2, P1_nP2;
//
//	// P2 = phi^2(P1)
//	Fp12::frobenius_map(P2, base, 2);
//
//#if defined(__CUDA_ARCH__)
//	if (threadIdx.x == 0 && blockIdx.x == 0) {
//		printf("GPU Trace step 1\n");
//	}
//#endif
//	// P1_P2 = P1 * P2
//	Fp12::mul(P1_P2, base, P2);
//
//	// P1_nP2 = P1 * conj(P2)
//	Fp12::mul(P1_nP2, base, P2, true);
//
//#if defined(__CUDA_ARCH__)
//	if (threadIdx.x == 0 && blockIdx.x == 0) {
//		printf("GPU Trace step 2\n");
//	}
//#endif
//
//	// 3. Init with MSB
//	len--;
//
//
//
//	// Initialization
//	{
//		auto pair = jsf.get(len);
//		OpSelect sel = get_op(pair.d1, pair.d2, base, P2, P1_P2, P1_nP2);
//#if defined(__CUDA_ARCH__)
//		if (sel.ptr == nullptr) {
//			printf("GPU Trace null\n");
//		}
//#endif
//		vec_copy(&out, sel.ptr, sizeof(Fp12));
//		if (sel.inv) Fp12::conjugate(out);
//	}
//
//#if defined(__CUDA_ARCH__)
//	if (threadIdx.x == 0 && blockIdx.x == 0) {
//		printf("GPU Trace Init: Len=%d, Out[0]=%llx\n", len, (unsigned long long)out[0][0][0][0]);
//	}
//#endif
//	
//	    // 4. Main Loop
//	    for (len--; len >= 0; --len) {
//	        Fp12::cyclotomic_sqr(out, out);
//	
//	        auto pair = jsf.get(len);
//	        int u = pair.d1;
//	        int v = pair.d2;
//	
//	        if (u == 0 && v == 0)
//	        {
//	#if defined(__CUDA_ARCH__)
//	            if (threadIdx.x == 0 && blockIdx.x == 0) {
//	                printf("GPU Trace Step %d (SQR only): Out[0]=%llx\n", len, (unsigned long long)out[0][0][0][0]);
//	            }
//	#endif
//	            continue;
//	        }
//	        // Minimized divergence: 
//	        // No stack allocations, no "conjugate" function call.
//	        // Just a boolean flag passed to mul.
//	        OpSelect sel = get_op(u, v, base, P2, P1_P2, P1_nP2);
//	        Fp12::mul(out, out, *sel.ptr, sel.inv);
//	
//	#if defined(__CUDA_ARCH__)
//	        if (threadIdx.x == 0 && blockIdx.x == 0) {
//	            printf("GPU Trace Step %d (MUL u=%d v=%d): Out[0]=%llx\n", len, u, v, (unsigned long long)out[0][0][0][0]);
//	        }
//	#endif
//	    }
//}

// Wrapper for Optimized JSF Exp
__global__ void kernel_cyclotomic_exp_glv_jsf(Fp12* out, const Fp12* base, const GLV_JSF_Exp* scheds) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= BATCH_SIZE) return;

	// Use modulo to cycle through available schedules (e.g. 8) 
	// or just assume 1:1 mapping if bench sets up distinct schedules per thread.
	// For this bench, we map 1:1.
	Fp12::cyclotomic_exp_glv_jsf(out[tid], base[tid], scheds[tid]);
}


// --- Fixtures for Data Setup ---

class Fp12CudaFixture : public benchmark::Fixture {
protected:
	Fp12* d_bases = nullptr;
	Fp12* d_out = nullptr;
	byte* d_scalars_flat = nullptr;
	GLV_JSF_Exp* d_schedules = nullptr;

	void SetUp(const ::benchmark::State& state) override {
		// Guard against multiple SetUp calls
		if (d_bases != nullptr) return;

		// Safe stack size for FP12 operations
		//cudaDeviceSetLimit(cudaLimitStackSize, 16384);

		// Allocations
		CUDA_CHECK(cudaMalloc(&d_bases, BATCH_SIZE * sizeof(Fp12)));
		CUDA_CHECK(cudaMalloc(&d_out, BATCH_SIZE * sizeof(Fp12)));

		// Host Data gen. Perform a pairing to get a valid GT element
		std::vector<Fp12> h_bases(BATCH_SIZE);

		Fp12 seed;
		vec_copy(&seed, &Fp12::one(), sizeof(seed));
		// Make seed non-trivial
		seed[0][0][1][0] = 2;
		seed[0][1][0][0] = 3;
		seed[1][2][1][0] = 4;
		//memset(&seed, 0xcc, sizeof(seed));

		Fp12 base;
		Fp12::final_exp(base, seed);
		BLST_ASSERT(Fp12::is_cyclotomic(base));

		// Simple distinct fake data
		for (int i = 0; i < BATCH_SIZE; i++) {
			h_bases[i] = base;
			Fp12::cyclotomic_sqr(base, base);
		}

		CUDA_CHECK(cudaMemcpy(d_bases, h_bases.data(), BATCH_SIZE * sizeof(Fp12), cudaMemcpyHostToDevice));
	}

	void TearDown(const ::benchmark::State& state) override {
		if (d_bases) { cudaFree(d_bases); d_bases = nullptr; }
		if (d_out) { cudaFree(d_out); d_out = nullptr; }
		if (d_scalars_flat) { cudaFree(d_scalars_flat); d_scalars_flat = nullptr; }
		if (d_schedules) { cudaFree(d_schedules); d_schedules = nullptr; }
	}

	// Helper to setup scalars specific to the test
	void SetupScalars_Basic() {
		std::vector<byte> h_scalars(BATCH_SIZE * 32);
		for (int i = 0; i < BATCH_SIZE; i++) {
			// Fill with pattern
			memset(&h_scalars[i * 32], 0xAA + (i % 5), 32);
		}
		CUDA_CHECK(cudaMalloc(&d_scalars_flat, BATCH_SIZE * 32));
		CUDA_CHECK(cudaMemcpy(d_scalars_flat, h_scalars.data(), BATCH_SIZE * 32, cudaMemcpyHostToDevice));
	}

	void SetupScalars_JSF() {
		std::vector<GLV_JSF_Exp> h_scheds(BATCH_SIZE);

		// Random generator for FULL 256-bit comparable loads
		// (Split into two ~128-bit scalars)
		std::random_device rd;
		std::mt19937_64 gen(rd());
		std::uniform_int_distribution<uint64_t> dis;

		// Create 8 template schedules with HEAVY workloads
		GLV_JSF_Exp templates[8];

		for (int t = 0; t < 8; t++) {
			// Generate random ~128-bit values for k1, k2
			limb_t k1[4] = { dis(gen), dis(gen), 0, 0 };
			limb_t k2[4] = { dis(gen), dis(gen), 0, 0 };

			// Recode on Host
			templates[t] = Fp12::glv_recode_jsf(k1, k2);

			// Sanity check (ensure length is reasonable, expected ~128-130)
			if (templates[t].sched_len < 120) {
				printf("WARNING: Generated short schedule len: %d\n", templates[t].sched_len);
			}
		}

		// Distribute
		for (int i = 0; i < BATCH_SIZE; i++) {
			h_scheds[i] = templates[i % 8];
		}

		CUDA_CHECK(cudaMalloc(&d_schedules, BATCH_SIZE * sizeof(GLV_JSF_Exp)));
		CUDA_CHECK(cudaMemcpy(d_schedules, h_scheds.data(), BATCH_SIZE * sizeof(GLV_JSF_Exp), cudaMemcpyHostToDevice));
	}

public:
	Fp12CudaFixture() : d_scalars_flat(nullptr), d_schedules(nullptr) {}
};


// --- Benchmarks ---

BENCHMARK_DEFINE_F(Fp12CudaFixture, BM_Cuda_CyclotomicExp_Basic)(benchmark::State& state) {
	SetupScalars_Basic();

	for (auto _ : state) {
		kernel_cyclotomic_exp_basic << <BATCH_SIZE / 128, 128 >> > (d_out, d_bases, d_scalars_flat);
		// CRITICAL: Synchronize to measure GPU time, not launch time
		cudaDeviceSynchronize();
	}
	// Items = Total Ops = Batches * Iterations
	state.SetItemsProcessed(state.iterations() * BATCH_SIZE);
}
BENCHMARK_REGISTER_F(Fp12CudaFixture, BM_Cuda_CyclotomicExp_Basic)
->Unit(benchmark::kMillisecond)
->UseRealTime(); // Essential for GPU timing


BENCHMARK_DEFINE_F(Fp12CudaFixture, BM_Cuda_CyclotomicExp_GLV_JSF)(benchmark::State& state) {
	SetupScalars_JSF();

	uint64_t h_verification = 0;
	for (auto _ : state) {
		kernel_cyclotomic_exp_glv_jsf << <BATCH_SIZE / 128, 128 >> > (d_out, d_bases, d_schedules);
		// CRITICAL: Synchronize to measure GPU time, not launch time
		cudaDeviceSynchronize();

		//// Force dependency: Copy back 1st byte of 1st result
		//// This ensures the kernel MUST run and write to d_out
		//CUDA_CHECK(cudaMemcpy(&h_verification, d_out, sizeof(uint64_t), cudaMemcpyDeviceToHost));

		//// Prevent host compiler optimization
		//benchmark::DoNotOptimize(h_verification);
	}
	state.SetItemsProcessed(state.iterations() * BATCH_SIZE);
}
BENCHMARK_REGISTER_F(Fp12CudaFixture, BM_Cuda_CyclotomicExp_GLV_JSF)
->Unit(benchmark::kMillisecond)
->UseRealTime();







































// ============================================================================
// GLV-4 EXPERIMENTAL IMPLEMENTATIONS
// ============================================================================

// 1. STANDARD TABLE IMPLEMENTATION (High Stack Usage)
// Logic copied from Fp12.cpp but adapted for pre-split inputs to avoid internal deps
BLST_HD void cyclotomic_exp_glv4_table_impl(Fp12& out, const Fp12& base, const int64_t k[4])
{
	// 1. Prepare Basis Points {P0, P1, P2, P3}
	uint64_t abs_k[4];
	Fp12 basis[4];

	basis[0] = base; // P0
	Fp12::frobenius_map(basis[1], base, 1); // P1
	Fp12::frobenius_map(basis[2], base, 2); // P2
	Fp12::frobenius_map(basis[3], base, 3); // P3

	// Handle signs
	for (int i = 0; i < 4; ++i) {
		if (k[i] < 0) {
			abs_k[i] = (uint64_t)(-k[i]);
			Fp12::conjugate(basis[i]);
		}
		else {
			abs_k[i] = (uint64_t)k[i];
		}
	}

	// 2. Precompute Lookup Table (16 points) - THE MEMORY HOG
	Fp12 T[16];

	// T[0] = 1
	vec_copy(&T[0], &Fp12::one(), sizeof(Fp12));

	// Powers of 2 indices
	T[1] = basis[0];
	T[2] = basis[1];
	T[4] = basis[2];
	T[8] = basis[3];

	// Fill rest
	Fp12::mul(T[3], T[1], T[2]);  // 1|2
	Fp12::mul(T[5], T[1], T[4]);  // 1|4
	Fp12::mul(T[6], T[2], T[4]);  // 2|4
	Fp12::mul(T[7], T[3], T[4]);  // 1|2|4

	Fp12::mul(T[9], T[1], T[8]);  // 1|8
	Fp12::mul(T[10], T[2], T[8]); // 2|8
	Fp12::mul(T[11], T[3], T[8]); // 1|2|8
	Fp12::mul(T[12], T[4], T[8]); // 4|8
	Fp12::mul(T[13], T[5], T[8]); // 1|4|8
	Fp12::mul(T[14], T[6], T[8]); // 2|4|8
	Fp12::mul(T[15], T[7], T[8]); // 1|2|4|8

	// 3. Main Loop
	int msb = 63;
	uint64_t combined = abs_k[0] | abs_k[1] | abs_k[2] | abs_k[3];
	while (msb >= 0 && ((combined >> msb) & 1) == 0) msb--;

	if (msb < 0) {
		vec_copy(&out, &Fp12::one(), sizeof(Fp12));
		return;
	}

	// Init with first window
	unsigned idx = 0;
	for (int i = 0; i < 4; ++i) if ((abs_k[i] >> msb) & 1) idx |= (1 << i);
	out = T[idx];

	for (int j = msb - 1; j >= 0; --j) {
		Fp12::cyclotomic_sqr(out, out);

		idx = 0;
		for (int i = 0; i < 4; ++i) if ((abs_k[i] >> j) & 1) idx |= (1 << i);

		if (idx != 0) Fp12::mul(out, out, T[idx]);
	}
}


// 2. OPTIMIZED ON-FLY IMPLEMENTATION (Low Stack Usage)
BLST_HD void cyclotomic_exp_glv4_onfly_impl(Fp12& out, const Fp12& base, const int64_t k[4])
{
	// 1. Prepare Basis locally (Stack: 4 * 576B = 2.3KB vs 10KB)
	Fp12 basis[4];
	basis[0] = base;
	Fp12::frobenius_map(basis[1], base, 1);
	Fp12::frobenius_map(basis[2], base, 2);
	Fp12::frobenius_map(basis[3], base, 3);

	uint64_t abs_k[4];

	for (int i = 0; i < 4; ++i) {
		if (k[i] < 0) {
			abs_k[i] = (uint64_t)(-k[i]);
			Fp12::conjugate(basis[i]);
		}
		else {
			abs_k[i] = (uint64_t)k[i];
		}
	}

	// 2. Main Loop
	int msb = 63;
	uint64_t combined = abs_k[0] | abs_k[1] | abs_k[2] | abs_k[3];
	while (msb >= 0 && ((combined >> msb) & 1) == 0) msb--;

	if (msb < 0) {
		vec_copy(&out, &Fp12::one(), sizeof(Fp12));
		return;
	}

	// Init Accumulator (Manual Unroll for first step)
	{
		bool first = true;
		for (int i = 0; i < 4; ++i) {
			if ((abs_k[i] >> msb) & 1) {
				if (first) {
					out = basis[i];
					first = false;
				}
				else {
					Fp12::mul(out, out, basis[i]);
				}
			}
		}
		// If msb>=0, at least one bit is set, so out is initialized
	}

	// Temp for step accumulation
	Fp12 term;

	for (int j = msb - 1; j >= 0; --j) {
		Fp12::cyclotomic_sqr(out, out);

		// Replace Table lookup with sparse multiplication sequence
		bool term_set = false;

		for (int i = 0; i < 4; ++i) {
			if ((abs_k[i] >> j) & 1) {
				if (!term_set) {
					term = basis[i];
					term_set = true;
				}
				else {
					Fp12::mul(term, term, basis[i]);
				}
			}
		}

		if (term_set) {
			Fp12::mul(out, out, term);
		}
	}
}


// --- GLV4 Kernels ---

// Structure for passing split scalars conveniently
struct GLV4_Scalars { int64_t k[4]; };

__global__ void kernel_cyclotomic_exp_glv4_table(Fp12* out, const Fp12* base, const GLV4_Scalars* scalars) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= BATCH_SIZE) return;

	// Use the scalars passed from host
	cyclotomic_exp_glv4_table_impl(out[tid], base[tid], scalars[tid].k);
}

__global__ void kernel_cyclotomic_exp_glv4_onfly(Fp12* out, const Fp12* base, const GLV4_Scalars* scalars) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= BATCH_SIZE) return;

	cyclotomic_exp_glv4_onfly_impl(out[tid], base[tid], scalars[tid].k);

	// Diagnostics for correctness (Output 0 check)
	if (tid == 0) {
		out[tid][0][0][0][0] += 0; // Ensure write
	}
}


// --- Benchmark Class Extensions ---
// Extending existing fixture with GLV4 setup

class Fp12CudaGLV4ExpFixture : public Fp12CudaFixture {
protected:
	GLV4_Scalars* d_glv4_scalars;

	void SetUp(const ::benchmark::State& state) override {
		Fp12CudaFixture::SetUp(state); // Call parent setup for bases

		// setup random scalars
		std::vector<GLV4_Scalars> h_scalars(BATCH_SIZE);
		std::random_device rd;
		std::mt19937_64 gen(rd());
		std::uniform_int_distribution<int64_t> dis; // Signed for GLV

		for (int i = 0; i < BATCH_SIZE; i++) {
			// Generate random 64-bit splits 
			// (Simulating GLV decomp output, approx 64-bits each)
			h_scalars[i].k[0] = dis(gen);
			h_scalars[i].k[1] = dis(gen);
			h_scalars[i].k[2] = dis(gen);
			h_scalars[i].k[3] = dis(gen);
		}

		CUDA_CHECK(cudaMalloc(&d_glv4_scalars, BATCH_SIZE * sizeof(GLV4_Scalars)));
		CUDA_CHECK(cudaMemcpy(d_glv4_scalars, h_scalars.data(), BATCH_SIZE * sizeof(GLV4_Scalars), cudaMemcpyHostToDevice));
	}

	void TearDown(const ::benchmark::State& state) override {
		cudaFree(d_glv4_scalars);
		Fp12CudaFixture::TearDown(state);
	}
};


// --- GLV4 Benchmarks ---

BENCHMARK_DEFINE_F(Fp12CudaGLV4ExpFixture, BM_Cuda_CyclotomicExp_GLV4_Table)(benchmark::State& state) {
	// 16KB Stack needed for array of 16 Fp12s + recursion overhead
	//cudaError_t err = cudaDeviceSetLimit(cudaLimitStackSize, 16384);
	//if (err != cudaSuccess) {
	//	state.SkipWithError("Could not set 16KB stack for Table Bench");
	//	return;
	//}

	uint64_t verify_sink;
	auto fsiz = sizeof(verify_sink);
	for (auto _ : state) {
		kernel_cyclotomic_exp_glv4_table << <BATCH_SIZE / 128, 128 >> > (d_out, d_bases, d_glv4_scalars);
		CUDA_CHECK(cudaMemcpy(&verify_sink, d_out, fsiz, cudaMemcpyDeviceToHost));
		benchmark::DoNotOptimize(verify_sink);
	}
	state.SetItemsProcessed(state.iterations() * BATCH_SIZE);
}
BENCHMARK_REGISTER_F(Fp12CudaGLV4ExpFixture, BM_Cuda_CyclotomicExp_GLV4_Table)
->Unit(benchmark::kMillisecond)->UseRealTime();


































// ============================================================================
// ========================= GLV-4 CUDA BENCHMARKS ============================
// ============================================================================

// 1. Kernel that performs runtime splitting (scalar -> GLV4 split -> Exp)
// Accepts a raw 256-bit scalar per thread (4 limbs).
__global__ void kernel_cyclotomic_exp_glv4_scalar(Fp12* out, const Fp12* base, const limb_t* scalars) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= BATCH_SIZE) return;

	// Each thread gets a pointer to its own 4-limb scalar
	// Note: We assume pointers in 'scalars' are properly aligned for limb_t access
	const limb_t* k = scalars + (tid * 4);

	// Calls the overload: cyclotomic_exp_glv4(Fp12& out, const Fp12& base, const limb_t* scalar)
	// This performs the GLV decomposition on the GPU.
	Fp12::cyclotomic_exp_glv4(out[tid], base[tid], k);
}

// 2. Kernel that uses offline pre-splitting (GLV4_Exp -> Exp)
// Accepts a precompiled GLV4_Exp struct per thread.
__global__ void kernel_cyclotomic_exp_glv4_props(Fp12* out, const Fp12* base, const GLV4_Exp* props) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= BATCH_SIZE) return;

	// Calls the overload: cyclotomic_exp_glv4(Fp12& out, const Fp12& base, const GLV4_Exp& props)
	// This skips decomposition and runs the schedule directly.
	Fp12::cyclotomic_exp_glv4(out[tid], base[tid], props[tid]);
}

// Specialized Fixture for GLV4 data setup

// Specialized Fixture for GLV4 data setup
class Fp12CudaGLV4Fixture : public Fp12CudaFixture {
protected:
	limb_t* d_scalars_limbs = nullptr;
	GLV4_Exp* d_glv4_props = nullptr;

public:

	void SetUp(const ::benchmark::State& state) override {
		if (d_scalars_limbs != nullptr || d_glv4_props != nullptr) return;

		// 1. MUST call parent SetUp first to allocate d_bases/d_out
		Fp12CudaFixture::SetUp(state);

		// 2. Set Stack Size
		// GLV4 uses 16 * sizeof(Fp12) ~ 9KB + basis ~ 2KB + overhead.
		// 16KB is sufficient.
			CUDA_CHECK(cudaDeviceSetLimit(cudaLimitStackSize, 16384));
	}

	void TearDown(const ::benchmark::State& state) override {
		// Free derived class resources FIRST
		if (d_scalars_limbs) { cudaFree(d_scalars_limbs); d_scalars_limbs = nullptr; }
		if (d_glv4_props) { cudaFree(d_glv4_props); d_glv4_props = nullptr; }

		// Then call parent TearDown
		Fp12CudaFixture::TearDown(state);
	}

	// Helper to generate a valid Fr scalar (reduced mod r)
	void generate_fr_scalar(limb_t* k) {
		static std::random_device rd;
		static std::mt19937_64 gen(rd());
		static std::uniform_int_distribution<uint64_t> dis;

		// Generate random bytes
		byte raw[32];
		uint64_t* p = (uint64_t*)raw;
		for (int i = 0; i < 4; ++i) p[i] = dis(gen);

		// Use Fr::from_le_bytes to ensure reduction mod r
		BLST_ERROR err;
		Fr scalar;
		// Create a "point" representation to ensure it's a valid field element
		// Fr::from_scalar performs the reduction if needed or checks bounds
		// Actually, blst_fr_from_scalar just copies, we want random valid Fr.
		// Let's use the raw bytes -> Fr conversion which usually handles reduction or masking.
		// If Fr isn't available, we assume the user has blst::Fr.

		// Simplest stable way without knowing internal Fr layout:
		// Mask the top bit to ensure it is likely smaller than r (r is ~255 bits for BLS12-381)
		// valid r: 0x73eda753... (255 bits). 
		// Masking top byte to 0x70 ensures < r for benchmark validity without full reduction logic.
		raw[31] &= 0x70;

		// Copy to limb_t array
		std::memcpy(k, raw, 32);
	}

	// Setup raw scalars for runtime split benchmark
	void SetupScalars_GLV4_Scalar() {
		if (d_scalars_limbs) return;

		std::vector<limb_t> h_scalars(BATCH_SIZE * 4);

		for (int i = 0; i < BATCH_SIZE; ++i) {
			generate_fr_scalar(&h_scalars[i * 4]);
		}

		CUDA_CHECK(cudaMalloc(&d_scalars_limbs, BATCH_SIZE * 4 * sizeof(limb_t)));
		CUDA_CHECK(cudaMemcpy(d_scalars_limbs, h_scalars.data(), h_scalars.size() * sizeof(limb_t), cudaMemcpyHostToDevice));
	}

	// Setup precompiled properties for offline split benchmark
	void SetupScalars_GLV4_Props() {
		if (d_glv4_props) return;

		std::vector<GLV4_Exp> h_props(BATCH_SIZE);

		// Compile random scalars on CPU to simulate offline processing
		for (int i = 0; i < BATCH_SIZE; ++i) {
			limb_t k[4];
			generate_fr_scalar(k);
			h_props[i] = Fp12::glv4_compile(k);
		}

		CUDA_CHECK(cudaMalloc(&d_glv4_props, BATCH_SIZE * sizeof(GLV4_Exp)));
		CUDA_CHECK(cudaMemcpy(d_glv4_props, h_props.data(), h_props.size() * sizeof(GLV4_Exp), cudaMemcpyHostToDevice));
	}
};

// Benchmark: Runtime Scalar Split + Exp
BENCHMARK_DEFINE_F(Fp12CudaGLV4Fixture, BM_Cuda_CyclotomicExp_GLV4_Scalar)(benchmark::State& state) {
	SetupScalars_GLV4_Scalar();

	for (auto _ : state) {
		kernel_cyclotomic_exp_glv4_scalar << <BATCH_SIZE / 128, 128 >> > (d_out, d_bases, d_scalars_limbs);
		cudaDeviceSynchronize();
	}
	state.SetItemsProcessed(state.iterations() * BATCH_SIZE);
}
BENCHMARK_REGISTER_F(Fp12CudaGLV4Fixture, BM_Cuda_CyclotomicExp_GLV4_Scalar)
->Unit(benchmark::kMillisecond)
->UseRealTime();

// Benchmark: Precompiled Split + Exp
BENCHMARK_DEFINE_F(Fp12CudaGLV4Fixture, BM_Cuda_CyclotomicExp_GLV4_Props)(benchmark::State& state) {
	SetupScalars_GLV4_Props();

	for (auto _ : state) {
		kernel_cyclotomic_exp_glv4_props << <BATCH_SIZE / 128, 128 >> > (d_out, d_bases, d_glv4_props);
		cudaDeviceSynchronize();
	}
	state.SetItemsProcessed(state.iterations() * BATCH_SIZE);
}
BENCHMARK_REGISTER_F(Fp12CudaGLV4Fixture, BM_Cuda_CyclotomicExp_GLV4_Props)
->Unit(benchmark::kMillisecond)
->UseRealTime();





























// ============================================================================
// GPU-OPTIMIZED GLV-4 IMPLEMENTATION (Uniform Control Flow)
// ============================================================================

/// GPU-optimized GLV4 schedule with fixed iteration count.
/// All threads execute exactly 64 iterations with identical control flow.
struct GLV4_Exp_GPU {
	uint8_t sign_mask;        ///< Bit i set = basis[i] was negated
	uint8_t pad[3];           ///< Align to 4 bytes
	uint8_t table_idx[64];    ///< Fixed 64 entries: 0=identity, 1-15=T[idx]
};

/// Compile a 256-bit scalar into a GPU-friendly fixed-length schedule.
/// @param k_in  Input scalar as 4 limbs (256 bits)
/// @return GPU-optimized schedule with exactly 64 steps
BLST_HD GLV4_Exp_GPU glv4_compile_gpu(const limb_t* k_in) {
	GLV4_Exp_GPU props;
	props.sign_mask = 0;
	memset(props.table_idx, 0, sizeof(props.table_idx));

	// 1. Decompose scalar using GLV-4 splitting
	int64_t k[4];
	Fp12::gt_split_scalar_4d(k, k_in);

	// 2. Extract signs and absolute values
	uint64_t abs_k[4];
	for (int i = 0; i < 4; i++) {
		if (k[i] < 0) {
			props.sign_mask |= (1 << i);
			abs_k[i] = (uint64_t)(-k[i]);
		}
		else {
			abs_k[i] = (uint64_t)k[i];
		}
	}

	// 3. Build fixed 64-entry schedule (bit 63 down to bit 0)
	// Each entry is a 4-bit index into T[0..15]
	// T[0] = identity (multiply is no-op when using T[0])
	for (int bit = 63; bit >= 0; --bit) {
		unsigned idx = 0;
		if ((abs_k[0] >> bit) & 1) idx |= 1;
		if ((abs_k[1] >> bit) & 1) idx |= 2;
		if ((abs_k[2] >> bit) & 1) idx |= 4;
		if ((abs_k[3] >> bit) & 1) idx |= 8;

		// Store in reverse order so we process MSB first
		props.table_idx[63 - bit] = (uint8_t)idx;
	}

	return props;
}

/// GPU-optimized exponentiation with uniform control flow.
/// Every thread executes exactly the same number of operations.
/// @param out   Output Fp12 element
/// @param base  Input base (must be in cyclotomic subgroup)
/// @param props Precompiled GPU schedule
BLST_HD void cyclotomic_exp_glv4_gpu_impl(Fp12& out, const Fp12& base, const GLV4_Exp_GPU& props) {
	// 1. Build basis {phi^0(base), phi^1(base), phi^2(base), phi^3(base)}
	//Fp12 basis[4];
	Fp12 T[16];
	T[1] = base;
	Fp12::frobenius_map(T[2], base, 1);
	Fp12::frobenius_map(T[4], base, 2);
	Fp12::frobenius_map(T[8], base, 3);

	// Apply signs from decomposition
	if (props.sign_mask & 1) Fp12::conjugate(T[1]);
	if (props.sign_mask & 2) Fp12::conjugate(T[2]);
	if (props.sign_mask & 4) Fp12::conjugate(T[4]);
	if (props.sign_mask & 8) Fp12::conjugate(T[8]);

	// 2. Build lookup table T[0..15]
	// CRITICAL: T[0] = identity, so mul(x, T[0]) is effectively x * 1 = x
	vec_copy(&T[0], &Fp12::one(), sizeof(Fp12));  // T[0] = 1 (identity)

	//T[1] = basis[0];
	//T[2] = basis[1];
	//T[4] = basis[2];
	//T[8] = basis[3];

	Fp12::mul(T[3], T[1], T[2]);   // basis[0] * basis[1]
	Fp12::mul(T[5], T[1], T[4]);   // basis[0] * basis[2]
	Fp12::mul(T[6], T[2], T[4]);   // basis[1] * basis[2]
	Fp12::mul(T[7], T[3], T[4]);   // basis[0] * basis[1] * basis[2]
	Fp12::mul(T[9], T[1], T[8]);   // basis[0] * basis[3]
	Fp12::mul(T[10], T[2], T[8]);   // basis[1] * basis[3]
	Fp12::mul(T[11], T[3], T[8]);   // basis[0] * basis[1] * basis[3]
	Fp12::mul(T[12], T[4], T[8]);   // basis[2] * basis[3]
	Fp12::mul(T[13], T[5], T[8]);   // basis[0] * basis[2] * basis[3]
	Fp12::mul(T[14], T[6], T[8]);   // basis[1] * basis[2] * basis[3]
	Fp12::mul(T[15], T[7], T[8]);   // all four

	// 3. Initialize accumulator with first table entry
	// Skip leading zeros by finding first non-zero index
	int start = 0;
	while (start < 64 && props.table_idx[start] == 0) start++;

	if (start >= 64) {
		// All zeros = scalar was 0
		vec_copy(&out, &Fp12::one(), sizeof(Fp12));
		return;
	}

	out = T[props.table_idx[start]];

	// 4. UNIFORM MAIN LOOP: Fixed iteration count from start+1 to 63
	// Every thread does the same number of squarings and multiplies
	// (though multiply by T[0]=identity is effectively a no-op mathematically)
	for (int i = start + 1; i < 64; ++i) {
		Fp12::cyclotomic_sqr(out, out);

		// ALWAYS multiply - T[0] is identity so this is safe
		// This eliminates the branch that causes warp divergence
		Fp12::mul(out, out, T[props.table_idx[i]]);
	}
}

/// Alternative: Fully uniform version (no early-start optimization)
/// Guaranteed identical control flow for all threads regardless of scalar value.
BLST_HD void cyclotomic_exp_glv4_gpu_uniform_impl(Fp12& out, const Fp12& base, const GLV4_Exp_GPU& props) {
	// 1. Build basis
	Fp12 basis[4];
	basis[0] = base;
	Fp12::frobenius_map(basis[1], base, 1);
	Fp12::frobenius_map(basis[2], base, 2);
	Fp12::frobenius_map(basis[3], base, 3);

	if (props.sign_mask & 1) Fp12::conjugate(basis[0]);
	if (props.sign_mask & 2) Fp12::conjugate(basis[1]);
	if (props.sign_mask & 4) Fp12::conjugate(basis[2]);
	if (props.sign_mask & 8) Fp12::conjugate(basis[3]);

	// 2. Build lookup table
	Fp12 T[16];
	vec_copy(&T[0], &Fp12::one(), sizeof(Fp12));

	T[1] = basis[0];
	T[2] = basis[1];
	T[4] = basis[2];
	T[8] = basis[3];

	Fp12::mul(T[3], T[1], T[2]);
	Fp12::mul(T[5], T[1], T[4]);
	Fp12::mul(T[6], T[2], T[4]);
	Fp12::mul(T[7], T[3], T[4]);
	Fp12::mul(T[9], T[1], T[8]);
	Fp12::mul(T[10], T[2], T[8]);
	Fp12::mul(T[11], T[3], T[8]);
	Fp12::mul(T[12], T[4], T[8]);
	Fp12::mul(T[13], T[5], T[8]);
	Fp12::mul(T[14], T[6], T[8]);
	Fp12::mul(T[15], T[7], T[8]);

	// 3. FULLY UNIFORM: Start with identity, always do 64 iterations
	vec_copy(&out, &Fp12::one(), sizeof(Fp12));

	// Unroll hint for better instruction-level parallelism
#pragma unroll 4
	for (int i = 0; i < 64; ++i) {
		Fp12::cyclotomic_sqr(out, out);
		Fp12::mul(out, out, T[props.table_idx[i]]);
	}
}

// --- CUDA Kernels ---

__global__ void kernel_cyclotomic_exp_glv4_gpu(Fp12* out, const Fp12* base, const GLV4_Exp_GPU* props) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= BATCH_SIZE) return;

	cyclotomic_exp_glv4_gpu_impl(out[tid], base[tid], props[tid]);
}

__global__ void kernel_cyclotomic_exp_glv4_gpu_uniform(Fp12* out, const Fp12* base, const GLV4_Exp_GPU* props) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= BATCH_SIZE) return;

	cyclotomic_exp_glv4_gpu_uniform_impl(out[tid], base[tid], props[tid]);
}

// Kernel with runtime compilation (for comparison)
__global__ void kernel_cyclotomic_exp_glv4_gpu_runtime(Fp12* out, const Fp12* base, const limb_t* scalars) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= BATCH_SIZE) return;

	const limb_t* k = scalars + (tid * 4);
	GLV4_Exp_GPU props = glv4_compile_gpu(k);
	cyclotomic_exp_glv4_gpu_impl(out[tid], base[tid], props);
}

// --- Benchmark Fixture ---

class Fp12CudaGLV4GPUFixture : public Fp12CudaFixture {
protected:
	limb_t* d_scalars_limbs = nullptr;
	GLV4_Exp_GPU* d_glv4_props_gpu = nullptr;

public:
	void SetUp(const ::benchmark::State& state) override {
		if (d_scalars_limbs != nullptr || d_glv4_props_gpu != nullptr) return;

		Fp12CudaFixture::SetUp(state);

		// GLV4 with 16-element table needs ~10KB stack
		CUDA_CHECK(cudaDeviceSetLimit(cudaLimitStackSize, 16384));
	}

	void TearDown(const ::benchmark::State& state) override {
		if (d_scalars_limbs) { cudaFree(d_scalars_limbs); d_scalars_limbs = nullptr; }
		if (d_glv4_props_gpu) { cudaFree(d_glv4_props_gpu); d_glv4_props_gpu = nullptr; }
		Fp12CudaFixture::TearDown(state);
	}

	void generate_fr_scalar_local(limb_t* k) {
		static std::random_device rd;
		static std::mt19937_64 gen(rd());
		static std::uniform_int_distribution<uint64_t> dis;

		byte raw[32];
		uint64_t* p = (uint64_t*)raw;
		for (int i = 0; i < 4; ++i) p[i] = dis(gen);
		raw[31] &= 0x70;  // Ensure < r
		std::memcpy(k, raw, 32);
	}

	void SetupScalars_GPU() {
		if (d_scalars_limbs) return;

		std::vector<limb_t> h_scalars(BATCH_SIZE * 4);
		for (int i = 0; i < BATCH_SIZE; ++i) {
			generate_fr_scalar_local(&h_scalars[i * 4]);
		}

		CUDA_CHECK(cudaMalloc(&d_scalars_limbs, BATCH_SIZE * 4 * sizeof(limb_t)));
		CUDA_CHECK(cudaMemcpy(d_scalars_limbs, h_scalars.data(),
			h_scalars.size() * sizeof(limb_t), cudaMemcpyHostToDevice));
	}

	void SetupProps_GPU() {
		if (d_glv4_props_gpu) return;

		std::vector<GLV4_Exp_GPU> h_props(BATCH_SIZE);

		for (int i = 0; i < BATCH_SIZE; ++i) {
			limb_t k[4];
			generate_fr_scalar_local(k);
			h_props[i] = glv4_compile_gpu(k);
		}

		CUDA_CHECK(cudaMalloc(&d_glv4_props_gpu, BATCH_SIZE * sizeof(GLV4_Exp_GPU)));
		CUDA_CHECK(cudaMemcpy(d_glv4_props_gpu, h_props.data(),
			h_props.size() * sizeof(GLV4_Exp_GPU), cudaMemcpyHostToDevice));
	}
};

// --- Benchmarks ---

// Precompiled GPU schedule (recommended for fixed scalars)
BENCHMARK_DEFINE_F(Fp12CudaGLV4GPUFixture, BM_Cuda_CyclotomicExp_GLV4_GPU_Precompiled)(benchmark::State& state) {
	SetupProps_GPU();

	for (auto _ : state) {
		kernel_cyclotomic_exp_glv4_gpu << <BATCH_SIZE / 128, 128 >> > (d_out, d_bases, d_glv4_props_gpu);
		cudaDeviceSynchronize();
	}
	state.SetItemsProcessed(state.iterations() * BATCH_SIZE);
}
BENCHMARK_REGISTER_F(Fp12CudaGLV4GPUFixture, BM_Cuda_CyclotomicExp_GLV4_GPU_Precompiled)
->Unit(benchmark::kMillisecond)
->UseRealTime();

// Fully uniform version (maximum warp coherence)
BENCHMARK_DEFINE_F(Fp12CudaGLV4GPUFixture, BM_Cuda_CyclotomicExp_GLV4_GPU_Uniform)(benchmark::State& state) {
	SetupProps_GPU();

	for (auto _ : state) {
		kernel_cyclotomic_exp_glv4_gpu_uniform << <BATCH_SIZE / 128, 128 >> > (d_out, d_bases, d_glv4_props_gpu);
		cudaDeviceSynchronize();
	}
	state.SetItemsProcessed(state.iterations() * BATCH_SIZE);
}
BENCHMARK_REGISTER_F(Fp12CudaGLV4GPUFixture, BM_Cuda_CyclotomicExp_GLV4_GPU_Uniform)
->Unit(benchmark::kMillisecond)
->UseRealTime();

// Runtime compilation on GPU (for variable scalars)
BENCHMARK_DEFINE_F(Fp12CudaGLV4GPUFixture, BM_Cuda_CyclotomicExp_GLV4_GPU_Runtime)(benchmark::State& state) {
	SetupScalars_GPU();

	for (auto _ : state) {
		kernel_cyclotomic_exp_glv4_gpu_runtime << <BATCH_SIZE / 128, 128 >> > (d_out, d_bases, d_scalars_limbs);
		cudaDeviceSynchronize();
	}
	state.SetItemsProcessed(state.iterations() * BATCH_SIZE);
}
BENCHMARK_REGISTER_F(Fp12CudaGLV4GPUFixture, BM_Cuda_CyclotomicExp_GLV4_GPU_Runtime)
->Unit(benchmark::kMillisecond)
->UseRealTime();

// ============================================================================
// REVISED GPU OPTIMIZATION: Remove restrictive launch bounds
// ============================================================================

// Remove __launch_bounds__ or use very conservative values
// The Fp12 operations are register-heavy by nature

/// Low-register version: compute table entries on-the-fly
/// Uses 4 basis elements instead of 16 table entries
BLST_HD void cyclotomic_exp_glv4_lowreg_impl(Fp12& out, const Fp12& base, const limb_t* scalar) {
	// 1. Decompose
	int64_t k[4];
	Fp12::gt_split_scalar_4d(k, scalar);

	// 2. Build basis only (4 elements vs 16)
	Fp12 basis[4];
	uint64_t abs_k[4];

	basis[0] = base;
	Fp12::frobenius_map(basis[1], base, 1);
	Fp12::frobenius_map(basis[2], base, 2);
	Fp12::frobenius_map(basis[3], base, 3);

	for (int i = 0; i < 4; ++i) {
		if (k[i] < 0) {
			abs_k[i] = (uint64_t)(-k[i]);
			Fp12::conjugate(basis[i]);
		}
		else {
			abs_k[i] = (uint64_t)k[i];
		}
	}

	// 3. Find MSB
	uint64_t combined = abs_k[0] | abs_k[1] | abs_k[2] | abs_k[3];
	int msb = 63;
	while (msb >= 0 && ((combined >> msb) & 1) == 0) msb--;

	if (msb < 0) {
		vec_copy(&out, &Fp12::one(), sizeof(Fp12));
		return;
	}

	// 4. Initialize with first non-zero combination
	{
		unsigned idx = 0;
		for (int i = 0; i < 4; ++i)
			if ((abs_k[i] >> msb) & 1) idx |= (1 << i);

		// Build T[idx] on the fly
		bool first = true;
		for (int i = 0; i < 4; ++i) {
			if (idx & (1 << i)) {
				if (first) { out = basis[i]; first = false; }
				else { Fp12::mul(out, out, basis[i]); }
			}
		}
	}

	// 5. Main loop - compute table entries on-demand
	for (int j = msb - 1; j >= 0; --j) {
		Fp12::cyclotomic_sqr(out, out);

		unsigned idx = 0;
		for (int i = 0; i < 4; ++i)
			if ((abs_k[i] >> j) & 1) idx |= (1 << i);

		if (idx != 0) {
			// Build T[idx] on the fly (max 3 extra muls, avg ~1.5)
			Fp12 term;
			bool first = true;
			for (int i = 0; i < 4; ++i) {
				if (idx & (1 << i)) {
					if (first) { term = basis[i]; first = false; }
					else { Fp12::mul(term, term, basis[i]); }
				}
			}
			Fp12::mul(out, out, term);
		}
	}
}

// NO launch_bounds - let compiler decide based on register usage
__global__ void kernel_cyclotomic_exp_glv4_lowreg(Fp12* out, const Fp12* base, const limb_t* scalars) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= BATCH_SIZE) return;

	const limb_t* k = scalars + (tid * 4);
	cyclotomic_exp_glv4_lowreg_impl(out[tid], base[tid], k);
}

// Simple wrapper - no launch_bounds
__global__ void kernel_cyclotomic_exp_glv4_smem(Fp12* out, const Fp12* base, const limb_t* scalars) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= BATCH_SIZE) return;

	const limb_t* k = scalars + (tid * 4);
	Fp12::cyclotomic_exp_glv4(out[tid], base[tid], k);
}

// --- Benchmark ---
BENCHMARK_DEFINE_F(Fp12CudaGLV4GPUFixture, BM_Cuda_CyclotomicExp_GLV4_LowReg)(benchmark::State& state) {
	SetupScalars_GPU();

	for (auto _ : state) {
		kernel_cyclotomic_exp_glv4_lowreg << <BATCH_SIZE / 128, 128 >> > (d_out, d_bases, d_scalars_limbs);
		cudaDeviceSynchronize();
	}
	state.SetItemsProcessed(state.iterations() * BATCH_SIZE);
}
BENCHMARK_REGISTER_F(Fp12CudaGLV4GPUFixture, BM_Cuda_CyclotomicExp_GLV4_LowReg)
->Unit(benchmark::kMillisecond)
->UseRealTime();