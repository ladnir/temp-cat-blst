#include <stdint.h>

#include "blst/common/defines.h"


template<typename Field>
BLST_DEVICE BLST_CPP_INLINE
void ptx_add_mod_u64(
	uint64_t* BLST_restrict r,
	const uint64_t* BLST_restrict a,
	const uint64_t* BLST_restrict b)
{

	constexpr auto n = sizeof(Field) / sizeof(uint64_t);
	constexpr auto BLST_restrict p = Field::modulus().data;

	if constexpr (n == 6)
	{

		uint64_t x0;
		uint64_t x1;
		uint64_t x2;
		uint64_t x3;
		uint64_t x4;
		uint64_t x5;

		uint64_t s0;
		uint64_t s1;
		uint64_t s2;
		uint64_t s3;
		uint64_t s4;
		uint64_t s5;

		uint64_t t0;
		uint64_t t1;
		uint64_t t2;
		uint64_t t3;
		uint64_t t4;
		uint64_t t5;


		// ===== a + b with carry -> x, carry1 in a GPR ===================================
		uint32_t carry1 = 0;
		asm volatile(
			"add.cc.u64  %0,  %7,  %13;\n\t"
			"addc.cc.u64 %1,  %8,  %14;\n\t"
			"addc.cc.u64 %2,  %9,  %15;\n\t"
			"addc.cc.u64 %3,  %10, %16;\n\t"
			"addc.cc.u64 %4,  %11, %17;\n\t"
			"addc.u64    %5,  %12, %18;\n\t"
			"mov.u32     %6,  0;\n\t"
			"addc.u32    %6,  %6,  %6;\n\t"   // carry1 = CC
			: "=&l"(x0),
			"=&l"(x1),
			"=&l"(x2),
			"=&l"(x3),
			"=&l"(x4),
			"=&l"(x5),
			"+r"(carry1)
			: "l"(a[0]), "l"(a[1]), "l"(a[2]), "l"(a[3]), "l"(a[4]), "l"(a[5]),
			"l"(b[0]), "l"(b[1]), "l"(b[2]), "l"(b[3]), "l"(b[4]), "l"(b[5])
			);

		t0 = x0;
		t1 = x1;
		t2 = x2;
		t3 = x3;
		t4 = x4;
		t5 = x5;

		// ===== x - p -> s (we’ll compute ge(x,p) separately) ============================
		asm volatile(
			"sub.cc.u64  %0,  %6,  %12;\n\t"
			"subc.cc.u64 %1,  %7,  %13;\n\t"
			"subc.cc.u64 %2,  %8,  %14;\n\t"
			"subc.cc.u64 %3,  %9,  %15;\n\t"
			"subc.cc.u64 %4,  %10, %16;\n\t"
			"subc.u64    %5,  %11, %17;\n\t"
			: "=&l"(s0),
			"=&l"(s1),
			"=&l"(s2),
			"=&l"(s3),
			"=&l"(s4),
			"=&l"(s5)
			: "l"(x0), "l"(x1), "l"(x2), "l"(x3), "l"(x4), "l"(x5),
			"l"(p[0]), "l"(p[1]), "l"(p[2]), "l"(p[3]), "l"(p[4]), "l"(p[5])
			);

		// ===== need_reduce = (carry1 != 0) || (x >= p) ==================================
		uint32_t need_reduce;
		asm volatile(
			// Inputs here: %1=carry1, %2..%7 = t0..t5, %8..%13 = p0..p5
			".reg .pred peq;\n\t"
			".reg .pred peq_next;\n\t"
			".reg .pred pgt;\n\t"
			".reg .pred ptmp;\n\t"
			".reg .pred pge;\n\t"

			"setp.gt.u64 pgt, %7,  %13;\n\t"   // compare high limb first: t5 > p5 ?
			"setp.eq.u64 peq, %7,  %13;\n\t"

			"setp.gt.u64 ptmp, %6, %12;\n\t"
			"and.pred    ptmp, ptmp, peq;\n\t"
			"or.pred     pgt,  pgt,  ptmp;\n\t"
			"setp.eq.u64 peq_next, %6, %12;\n\t"
			"and.pred    peq, peq, peq_next;\n\t"

			"setp.gt.u64 ptmp, %5, %11;\n\t"
			"and.pred    ptmp, ptmp, peq;\n\t"
			"or.pred     pgt,  pgt,  ptmp;\n\t"
			"setp.eq.u64 peq_next, %5, %11;\n\t"
			"and.pred    peq, peq, peq_next;\n\t"

			"setp.gt.u64 ptmp, %4, %10;\n\t"
			"and.pred    ptmp, ptmp, peq;\n\t"
			"or.pred     pgt,  pgt,  ptmp;\n\t"
			"setp.eq.u64 peq_next, %4, %10;\n\t"
			"and.pred    peq, peq, peq_next;\n\t"

			"setp.gt.u64 ptmp, %3,  %9;\n\t"
			"and.pred    ptmp, ptmp, peq;\n\t"
			"or.pred     pgt,  pgt,  ptmp;\n\t"
			"setp.eq.u64 peq_next, %3,  %9;\n\t"
			"and.pred    peq, peq, peq_next;\n\t"

			"setp.gt.u64 ptmp, %2,  %8;\n\t"
			"and.pred    ptmp, ptmp, peq;\n\t"
			"or.pred     pgt,  pgt,  ptmp;\n\t"
			"setp.eq.u64 peq_next, %2,  %8;\n\t"
			"and.pred    peq, peq, peq_next;\n\t"

			"or.pred  pge, pgt, peq;\n\t"
			"selp.u32 %0, 1, 0, pge;\n\t"      // %0 = ge ? 1 : 0
			"or.b32   %0, %0, %1;\n\t"         // need_reduce |= carry1
			: "=&r"(need_reduce)
			: "r"(carry1),
			"l"(t0), "l"(t1), "l"(t2), "l"(t3), "l"(t4), "l"(t5),
			"l"(p[0]), "l"(p[1]), "l"(p[2]), "l"(p[3]), "l"(p[4]), "l"(p[5])
			);

		// ===== blend: r = need_reduce ? s : t ==========================================
		uint64_t o0;
		uint64_t o1;
		uint64_t o2;
		uint64_t o3;
		uint64_t o4;
		uint64_t o5;

		asm volatile(
			".reg .pred psel;\n\t"
			"setp.ne.u32 psel, %12, 0;\n\t"
			"selp.u64 %0, %6,  %13, psel;\n\t"
			"selp.u64 %1, %7,  %14, psel;\n\t"
			"selp.u64 %2, %8,  %15, psel;\n\t"
			"selp.u64 %3, %9,  %16, psel;\n\t"
			"selp.u64 %4, %10, %17, psel;\n\t"
			"selp.u64 %5, %11, %18, psel;\n\t"
			: "=&l"(o0),
			"=&l"(o1),
			"=&l"(o2),
			"=&l"(o3),
			"=&l"(o4),
			"=&l"(o5)
			: "l"(s0), "l"(s1), "l"(s2), "l"(s3), "l"(s4), "l"(s5),
			"r"(need_reduce),
			"l"(t0), "l"(t1), "l"(t2), "l"(t3), "l"(t4), "l"(t5)
			);

		r[0] = o0;
		r[1] = o1;
		r[2] = o2;
		r[3] = o3;
		r[4] = o4;
		r[5] = o5;
	}
	else if constexpr (n == 4)
	{
		uint64_t x0;
		uint64_t x1;
		uint64_t x2;
		uint64_t x3;

		uint64_t s0;
		uint64_t s1;
		uint64_t s2;
		uint64_t s3;

		uint64_t t0;
		uint64_t t1;
		uint64_t t2;
		uint64_t t3;

		uint32_t carry1;
		carry1 = 0;

		// x = a + b ; carry1 = final carry bit
		asm volatile(
				"add.cc.u64  %0,  %5,  %9;\n\t"
				"addc.cc.u64 %1,  %6,  %10;\n\t"
				"addc.cc.u64 %2,  %7,  %11;\n\t"
				"addc.u64    %3,  %8,  %12;\n\t"
				"mov.u32     %4,  0;\n\t"
				"addc.u32    %4,  %4,  %4;\n\t"
				: "=&l"(x0),
				"=&l"(x1),
				"=&l"(x2),
				"=&l"(x3),
				"+r"(carry1)
				: "l"(a[0]), "l"(a[1]), "l"(a[2]), "l"(a[3]),
				"l"(b[0]), "l"(b[1]), "l"(b[2]), "l"(b[3])
				);

		t0 = x0;
		t1 = x1;
		t2 = x2;
		t3 = x3;

		// s = x - p
		asm volatile(
				"sub.cc.u64  %0,  %4,  %8;\n\t"
				"subc.cc.u64 %1,  %5,  %9;\n\t"
				"subc.cc.u64 %2,  %6,  %10;\n\t"
				"subc.u64    %3,  %7,  %11;\n\t"
				: "=&l"(s0),
				"=&l"(s1),
				"=&l"(s2),
				"=&l"(s3)
				: "l"(x0), "l"(x1), "l"(x2), "l"(x3),
				"l"(p[0]), "l"(p[1]), "l"(p[2]), "l"(p[3])
				);

		// need_reduce = carry1 || (x >= p)
		uint32_t need_reduce;

		asm volatile(
				".reg .pred peq;\n\t"
				".reg .pred peq2;\n\t"
				".reg .pred pgt;\n\t"
				".reg .pred pt;\n\t"
				".reg .pred pge;\n\t"

				// compare from most significant limb: t3..t0 vs p3..p0
				"setp.gt.u64 pgt, %5,  %9;\n\t"
				"setp.eq.u64 peq, %5,  %9;\n\t"

				"setp.gt.u64 pt,  %4,  %8;\n\t"
				"and.pred    pt,  pt,  peq;\n\t"
				"or.pred     pgt, pgt, pt;\n\t"
				"setp.eq.u64 peq2,%4,  %8;\n\t"
				"and.pred    peq, peq, peq2;\n\t"

				"setp.gt.u64 pt,  %3,  %7;\n\t"
				"and.pred    pt,  pt,  peq;\n\t"
				"or.pred     pgt, pgt, pt;\n\t"
				"setp.eq.u64 peq2,%3,  %7;\n\t"
				"and.pred    peq, peq, peq2;\n\t"

				"setp.gt.u64 pt,  %2,  %6;\n\t"
				"and.pred    pt,  pt,  peq;\n\t"
				"or.pred     pgt, pgt, pt;\n\t"
				"setp.eq.u64 peq2,%2,  %6;\n\t"
				"and.pred    peq, peq, peq2;\n\t"

				"or.pred     pge, pgt, peq;\n\t"
				"selp.u32    %0, 1, 0, pge;\n\t"
				"or.b32      %0, %0, %1;\n\t"
				: "=&r"(need_reduce)
				: "r"(carry1),
				"l"(t0), "l"(t1), "l"(t2), "l"(t3),
				"l"(p[0]), "l"(p[1]), "l"(p[2]), "l"(p[3])
				);

		// r = need_reduce ? s : t
		uint64_t o0;
		uint64_t o1;
		uint64_t o2;
		uint64_t o3;

		asm volatile(
				".reg .pred psel;\n\t"
				"setp.ne.u32 psel, %8, 0;\n\t"
				"selp.u64 %0, %4, %9,  psel;\n\t"
				"selp.u64 %1, %5, %10, psel;\n\t"
				"selp.u64 %2, %6, %11, psel;\n\t"
				"selp.u64 %3, %7, %12, psel;\n\t"
				: "=&l"(o0),
				"=&l"(o1),
				"=&l"(o2),
				"=&l"(o3)
				: "l"(s0), "l"(s1), "l"(s2), "l"(s3),
				"r"(need_reduce),
				"l"(t0), "l"(t1), "l"(t2), "l"(t3)
				);

		r[0] = o0;
		r[1] = o1;
		r[2] = o2;
		r[3] = o3;
	}
	else
	{
		// static assert to trigger a compile-time error for unsupported Field sizes
		static_assert(n == 6 || n == 4, "ptx_add_mod_u64 only supports 256 or 384-bit fields");
	}
}


#ifdef BLST_ADD_MOD_384_U32_PTX

__device__ __forceinline__
void add_mod_384_u32_ptx(
	uint32_t* BLST_restrict r,   // 12 limbs, little-endian
	const uint32_t* BLST_restrict a,
	const uint32_t* BLST_restrict b,
	const uint32_t* BLST_restrict p)
{
	uint32_t x0;
	uint32_t x1;
	uint32_t x2;
	uint32_t x3;
	uint32_t x4;
	uint32_t x5;
	uint32_t x6;
	uint32_t x7;
	uint32_t x8;
	uint32_t x9;
	uint32_t x10;
	uint32_t x11;

	uint32_t s0;
	uint32_t s1;
	uint32_t s2;
	uint32_t s3;
	uint32_t s4;
	uint32_t s5;
	uint32_t s6;
	uint32_t s7;
	uint32_t s8;
	uint32_t s9;
	uint32_t s10;
	uint32_t s11;

	uint32_t t0;
	uint32_t t1;
	uint32_t t2;
	uint32_t t3;
	uint32_t t4;
	uint32_t t5;
	uint32_t t6;
	uint32_t t7;
	uint32_t t8;
	uint32_t t9;
	uint32_t t10;
	uint32_t t11;

	uint32_t carry1;
	carry1 = 0;

	// x = a + b, capture final carry into carry1
	asm volatile
		(
			"add.cc.u32  %0,  %12, %24;\n\t"
			"addc.cc.u32 %1,  %13, %25;\n\t"
			"addc.cc.u32 %2,  %14, %26;\n\t"
			"addc.cc.u32 %3,  %15, %27;\n\t"
			"addc.cc.u32 %4,  %16, %28;\n\t"
			"addc.cc.u32 %5,  %17, %29;\n\t"
			"addc.cc.u32 %6,  %18, %30;\n\t"
			"addc.cc.u32 %7,  %19, %31;\n\t"
			"addc.cc.u32 %8,  %20, %32;\n\t"
			"addc.cc.u32 %9,  %21, %33;\n\t"
			"addc.cc.u32 %10, %22, %34;\n\t"
			"addc.u32    %11, %23, %35;\n\t"
			"mov.u32     %36, 0;\n\t"
			"addc.u32    %36, %36, %36;\n\t"  // carry1 = CC
			: "=&r"(x0),
			"=&r"(x1),
			"=&r"(x2),
			"=&r"(x3),
			"=&r"(x4),
			"=&r"(x5),
			"=&r"(x6),
			"=&r"(x7),
			"=&r"(x8),
			"=&r"(x9),
			"=&r"(x10),
			"=&r"(x11),
			"+r"(carry1)
			: "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(a[3]),
			"r"(a[4]), "r"(a[5]), "r"(a[6]), "r"(a[7]),
			"r"(a[8]), "r"(a[9]), "r"(a[10]), "r"(a[11]),
			"r"(b[0]), "r"(b[1]), "r"(b[2]), "r"(b[3]),
			"r"(b[4]), "r"(b[5]), "r"(b[6]), "r"(b[7]),
			"r"(b[8]), "r"(b[9]), "r"(b[10]), "r"(b[11])
			);

	t0 = x0;
	t1 = x1;
	t2 = x2;
	t3 = x3;
	t4 = x4;
	t5 = x5;
	t6 = x6;
	t7 = x7;
	t8 = x8;
	t9 = x9;
	t10 = x10;
	t11 = x11;

	// s = x - p
	asm volatile
		(
			"sub.cc.u32  %0,  %12, %24;\n\t"
			"subc.cc.u32 %1,  %13, %25;\n\t"
			"subc.cc.u32 %2,  %14, %26;\n\t"
			"subc.cc.u32 %3,  %15, %27;\n\t"
			"subc.cc.u32 %4,  %16, %28;\n\t"
			"subc.cc.u32 %5,  %17, %29;\n\t"
			"subc.cc.u32 %6,  %18, %30;\n\t"
			"subc.cc.u32 %7,  %19, %31;\n\t"
			"subc.cc.u32 %8,  %20, %32;\n\t"
			"subc.cc.u32 %9,  %21, %33;\n\t"
			"subc.cc.u32 %10, %22, %34;\n\t"
			"subc.u32    %11, %23, %35;\n\t"
			: "=&r"(s0),
			"=&r"(s1),
			"=&r"(s2),
			"=&r"(s3),
			"=&r"(s4),
			"=&r"(s5),
			"=&r"(s6),
			"=&r"(s7),
			"=&r"(s8),
			"=&r"(s9),
			"=&r"(s10),
			"=&r"(s11)
			: "r"(x0), "r"(x1), "r"(x2), "r"(x3),
			"r"(x4), "r"(x5), "r"(x6), "r"(x7),
			"r"(x8), "r"(x9), "r"(x10), "r"(x11),
			"r"(p[0]), "r"(p[1]), "r"(p[2]), "r"(p[3]),
			"r"(p[4]), "r"(p[5]), "r"(p[6]), "r"(p[7]),
			"r"(p[8]), "r"(p[9]), "r"(p[10]), "r"(p[11])
			);

	// need_reduce = carry1 || (x >= p)
	// Lex compare high → low: (x11..x0) vs (p11..p0), all predicated and branchless.
	uint32_t need_reduce;
	asm volatile
		(
			".reg .pred peq;\n\t"
			".reg .pred peq2;\n\t"
			".reg .pred pgt;\n\t"
			".reg .pred pt;\n\t"
			".reg .pred pge;\n\t"

			"setp.gt.u32 pgt, %13, %25;\n\t"
			"setp.eq.u32 peq, %13, %25;\n\t"

			"setp.gt.u32 pt, %12, %24;\n\t"
			"and.pred    pt, pt, peq;\n\t"
			"or.pred     pgt, pgt, pt;\n\t"
			"setp.eq.u32 peq2, %12, %24;\n\t"
			"and.pred    peq, peq, peq2;\n\t"

			"setp.gt.u32 pt, %11, %23;\n\t"
			"and.pred    pt, pt, peq;\n\t"
			"or.pred     pgt, pgt, pt;\n\t"
			"setp.eq.u32 peq2, %11, %23;\n\t"
			"and.pred    peq, peq, peq2;\n\t"

			"setp.gt.u32 pt, %10, %22;\n\t"
			"and.pred    pt, pt, peq;\n\t"
			"or.pred     pgt, pgt, pt;\n\t"
			"setp.eq.u32 peq2, %10, %22;\n\t"
			"and.pred    peq, peq, peq2;\n\t"

			"setp.gt.u32 pt, %9, %21;\n\t"
			"and.pred    pt, pt, peq;\n\t"
			"or.pred     pgt, pgt, pt;\n\t"
			"setp.eq.u32 peq2, %9, %21;\n\t"
			"and.pred    peq, peq, peq2;\n\t"

			"setp.gt.u32 pt, %8, %20;\n\t"
			"and.pred    pt, pt, peq;\n\t"
			"or.pred     pgt, pgt, pt;\n\t"
			"setp.eq.u32 peq2, %8, %20;\n\t"
			"and.pred    peq, peq, peq2;\n\t"

			"setp.gt.u32 pt, %7, %19;\n\t"
			"and.pred    pt, pt, peq;\n\t"
			"or.pred     pgt, pgt, pt;\n\t"
			"setp.eq.u32 peq2, %7, %19;\n\t"
			"and.pred    peq, peq, peq2;\n\t"

			"setp.gt.u32 pt, %6, %18;\n\t"
			"and.pred    pt, pt, peq;\n\t"
			"or.pred     pgt, pgt, pt;\n\t"
			"setp.eq.u32 peq2, %6, %18;\n\t"
			"and.pred    peq, peq, peq2;\n\t"

			"setp.gt.u32 pt, %5, %17;\n\t"
			"and.pred    pt, pt, peq;\n\t"
			"or.pred     pgt, pgt, pt;\n\t"
			"setp.eq.u32 peq2, %5, %17;\n\t"
			"and.pred    peq, peq, peq2;\n\t"

			"setp.gt.u32 pt, %4, %16;\n\t"
			"and.pred    pt, pt, peq;\n\t"
			"or.pred     pgt, pgt, pt;\n\t"
			"setp.eq.u32 peq2, %4, %16;\n\t"
			"and.pred    peq, peq, peq2;\n\t"

			"setp.gt.u32 pt, %3, %15;\n\t"
			"and.pred    pt, pt, peq;\n\t"
			"or.pred     pgt, pgt, pt;\n\t"
			"setp.eq.u32 peq2, %3, %15;\n\t"
			"and.pred    peq, peq, peq2;\n\t"

			"setp.gt.u32 pt, %2, %14;\n\t"
			"and.pred    pt, pt, peq;\n\t"
			"or.pred     pgt, pgt, pt;\n\t"
			"setp.eq.u32 peq2, %2, %14;\n\t"
			"and.pred    peq, peq, peq2;\n\t"

			"or.pred     pge, pgt, peq;\n\t"
			"selp.u32    %0, 1, 0, pge;\n\t"
			"or.b32      %0, %0, %1;\n\t"
			: "=&r"(need_reduce)
			: "r"(carry1),
			"r"(t0), "r"(t1), "r"(t2), "r"(t3),
			"r"(t4), "r"(t5), "r"(t6), "r"(t7),
			"r"(t8), "r"(t9), "r"(t10), "r"(t11),
			"r"(p[0]), "r"(p[1]), "r"(p[2]), "r"(p[3]),
			"r"(p[4]), "r"(p[5]), "r"(p[6]), "r"(p[7]),
			"r"(p[8]), "r"(p[9]), "r"(p[10]), "r"(p[11])
			);

	uint32_t o0;
	uint32_t o1;
	uint32_t o2;
	uint32_t o3;
	uint32_t o4;
	uint32_t o5;
	uint32_t o6;
	uint32_t o7;
	uint32_t o8;
	uint32_t o9;
	uint32_t o10;
	uint32_t o11;

	asm volatile
		(
			".reg .pred psel;\n\t"
			"setp.ne.u32 psel, %25, 0;\n\t"
			"selp.u32 %0,  %13, %2,  psel;\n\t"
			"selp.u32 %1,  %14, %3,  psel;\n\t"
			"selp.u32 %2,  %15, %4,  psel;\n\t"
			"selp.u32 %3,  %16, %5,  psel;\n\t"
			"selp.u32 %4,  %17, %6,  psel;\n\t"
			"selp.u32 %5,  %18, %7,  psel;\n\t"
			"selp.u32 %6,  %19, %8,  psel;\n\t"
			"selp.u32 %7,  %20, %9,  psel;\n\t"
			"selp.u32 %8,  %21, %10, psel;\n\t"
			"selp.u32 %9,  %22, %11, psel;\n\t"
			"selp.u32 %10, %23, %12, psel;\n\t"
			"selp.u32 %11, %24, %1,  psel;\n\t"
			: "=&r"(o0),
			"=&r"(o1),
			"=&r"(o2),
			"=&r"(o3),
			"=&r"(o4),
			"=&r"(o5),
			"=&r"(o6),
			"=&r"(o7),
			"=&r"(o8),
			"=&r"(o9),
			"=&r"(o10),
			"=&r"(o11)
			: "r"(t0), "r"(s0),
			"r"(t1), "r"(s1),
			"r"(t2), "r"(s2),
			"r"(t3), "r"(s3),
			"r"(t4), "r"(s4),
			"r"(t5), "r"(s5),
			"r"(t6), "r"(s6),
			"r"(t7), "r"(s7),
			"r"(t8), "r"(s8),
			"r"(t9), "r"(s9),
			"r"(t10), "r"(s10),
			"r"(t11), "r"(s11),
			"r"(need_reduce)
			);

	r[0] = o0;
	r[1] = o1;
	r[2] = o2;
	r[3] = o3;
	r[4] = o4;
	r[5] = o5;
	r[6] = o6;
	r[7] = o7;
	r[8] = o8;
	r[9] = o9;
	r[10] = o10;
	r[11] = o11;
}

#endif // if BLST_ADD_MOD_384_U32_PTX