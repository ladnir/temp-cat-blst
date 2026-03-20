#include <stdint.h>

#include "blst/common/defines.h"

#include "blst/backend/backend_portable.h"

namespace blst
{

	// write limb array to buffer as hex string for debugging
	BLST_HD BLST_INLINE
		void print_limb(const char* buffer, const uint64_t* v, size_t n) {
		auto bytes = n * sizeof(uint64_t);
		char* buf_ptr = (char*)buffer;
		char* src_ptr = (char*)v + bytes - 1;
		for (size_t i = 0; i < bytes; i++)
		{
			unsigned char b = *src_ptr--;
			auto lo = b & 15;
			auto hi = b >> 4;
			*buf_ptr++ = hi > 9 ? 'A' + (hi - 10) : '0' + hi;
			*buf_ptr++ = lo > 9 ? 'A' + (lo - 10) : '0' + lo;
			if ((i % 32) == 31)
				*buf_ptr++ = ' ';
		}
		*buf_ptr = 0;
	}

	template<typename Field>
	BLST_DEVICE BLST_CPP_INLINE
		void ptx_add_mod_u64(
			uint64_t* r,
			const uint64_t* a,
			const uint64_t* b)
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
				"add.cc.u64   %0,  %7,  %13;\n\t"
				"addc.cc.u64  %1,  %8,  %14;\n\t"
				"addc.cc.u64  %2,  %9,  %15;\n\t"
				"addc.cc.u64  %3,  %10, %16;\n\t"
				"addc.cc.u64  %4,  %11, %17;\n\t"
				"addc.cc.u64  %5,  %12, %18;\n\t"   // final must be .cc
				"mov.u32      %6,  0;\n\t"
				"addc.u32     %6,  %6,  %6;\n\t"    // carry1 = CC
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
			uint32_t noborrow;
			asm volatile(
				// %0..%5  -> s0..s5 (64-bit outputs)
				// %6      -> noborrow (32-bit output)  **NOT** used by the sub lines
				// %7..%12 -> x0..x5 (64-bit inputs)
				// %13..%18-> p0..p5 (64-bit inputs)

				"sub.cc.u64   %0,  %7,  %13;\n\t"
				"subc.cc.u64  %1,  %8,  %14;\n\t"
				"subc.cc.u64  %2,  %9,  %15;\n\t"
				"subc.cc.u64  %3,  %10, %16;\n\t"
				"subc.cc.u64  %4,  %11, %17;\n\t"
				"subc.cc.u64  %5,  %12, %18;\n\t"   // CC now = !borrow
				"mov.u32      %6,  0;\n\t"
				"addc.u32     %6,  %6,  %6;\n\t"    // noborrow = CC (1=no borrow)
				: "=&l"(s0),
				"=&l"(s1),
				"=&l"(s2),
				"=&l"(s3),
				"=&l"(s4),
				"=&l"(s5),
				"=&r"(noborrow)                   // %6 — not referenced above
				: "l"(x0), "l"(x1), "l"(x2), "l"(x3), "l"(x4), "l"(x5),
				"l"(p[0]), "l"(p[1]), "l"(p[2]), "l"(p[3]), "l"(p[4]), "l"(p[5])
				);

			// ===== need_reduce = (carry1 != 0) || (x >= p) ==================================
			uint32_t need_reduce;
			asm volatile("or.b32 %0, %1, %2;" : "=&r"(need_reduce) : "r"(carry1), "r"(noborrow));


			//// ===== blend: r = need_reduce ? s : t ==========================================
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


			// ===== a + b with carry -> x, carry1 in a GPR ===================================
			uint32_t carry1 = 0;
			asm volatile(
				"add.cc.u64   %0,  %5,  %9;\n\t"
				"addc.cc.u64  %1,  %6,  %10;\n\t"
				"addc.cc.u64  %2,  %7,  %11;\n\t"
				"addc.cc.u64  %3,  %8,  %12;\n\t"  // last limb must be .cc to update CC
				"mov.u32      %4,  0;\n\t"
				"addc.u32     %4,  %4,  %4;\n\t"   // carry1 = CC
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

			// ===== x - p -> s (we’ll compute ge(x,p) separately) ============================
			uint32_t noborrow;
			asm volatile(
				// %0..%3  -> s0..s3 (64-bit outputs)
				// %4      -> noborrow (32-bit output)  **not used** in sub lines
				// %5..%8  -> x0..x3 (64-bit inputs)
				// %9..%12 -> p0..p3 (64-bit inputs)

				"sub.cc.u64   %0, %5,  %9;\n\t"
				"subc.cc.u64  %1, %6,  %10;\n\t"
				"subc.cc.u64  %2, %7,  %11;\n\t"
				"subc.cc.u64  %3, %8,  %12;\n\t"   // final must be .cc
				"mov.u32      %4, 0;\n\t"
				"addc.u32     %4, %4, %4;\n\t"     // noborrow = CC (1 = no borrow)
				: "=&l"(s0),
				"=&l"(s1),
				"=&l"(s2),
				"=&l"(s3),
				"=&r"(noborrow)
				: "l"(x0), "l"(x1), "l"(x2), "l"(x3),
				"l"(p[0]), "l"(p[1]), "l"(p[2]), "l"(p[3])
				);


			// ===== need_reduce = (carry1 != 0) || (x >= p) ==================================
			uint32_t need_reduce;
			asm volatile("or.b32 %0, %1, %2;" : "=&r"(need_reduce) : "r"(carry1), "r"(noborrow));


			// ===== blend: r = need_reduce ? s : t ==========================================
			uint64_t o0;
			uint64_t o1;
			uint64_t o2;
			uint64_t o3;

			asm volatile(
				".reg .pred psel;\n\t"
				"setp.ne.u32 psel, %8, 0;\n\t"
				"selp.u64 %0, %4,  %9,  psel;\n\t"
				"selp.u64 %1, %5,  %10, psel;\n\t"
				"selp.u64 %2, %6,  %11, psel;\n\t"
				"selp.u64 %3, %7,  %12, psel;\n\t"
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
//
//#ifndef NDEBUG
//		static_assert(std::is_same_v<blst::limb_t, uint64_t>, "limb_t must be uint64_t for ptx_add_mod_u64");
//		uint64_t tt[n];
//		blst::cc_add_mod_n<n>(tt, a, b, (uint64_t*)Field::modulus().data);
//		for (size_t i = 0; i < n; i++) {
//			if (r[i] != tt[i])
//			{
//				printf("ptx_add_mod_u64<%d> error: \nact = ", (int)n);
//				for (size_t j = 0; j < n; j++)
//				{
//					// decimal, fixed width with leading zeros
//					printf("%016llx ", (unsigned long long)r[j]);
//				}
//				printf("\nexp = ");
//				for (size_t j = 0; j < n; j++)
//				{
//					// decimal
//					printf("%016llx ", (unsigned long long)tt[j]);
//				}
//				printf("\n   a = ");
//				for (size_t j = 0; j < n; j++)
//				{
//					// decimal
//					printf("%016llx ", (unsigned long long)a[j]);
//				}
//				printf("\n   b = ");
//				for (size_t j = 0; j < n; j++)
//				{
//					// decimal
//					printf("%016llx ", (unsigned long long)b[j]);
//				}
//				printf("\n ");
//			}
//			BLST_ASSERT(r[i] == tt[i]);
//		}
//#endif

	}
	template<typename Field>
	BLST_DEVICE BLST_CPP_INLINE
		void ptx_sub_mod_u64(
			uint64_t* r,
			const uint64_t* a,
			const uint64_t* b)
	{
		constexpr int n = sizeof(Field) / sizeof(uint64_t);
		constexpr auto BLST_restrict p = Field::modulus().data;

		if constexpr (n == 6)
		{
			uint64_t r0;
			uint64_t r1;
			uint64_t r2;
			uint64_t r3;
			uint64_t r4;
			uint64_t r5;

			// ---- Pass 1: r = a - b ; capture final borrow (via CC -> NOT-borrow) ----
			// After the last subc, CC = !borrow. We extract CC with addc.u32 tmp,0,0 and then XOR 1.
			uint32_t borrow;

			// outputs: r0=%0 r1=%1 r2=%2 r3=%3 r4=%4 r5=%5 borrow=%6
			// inputs:  a0=%7 a1=%8 a2=%9 a3=%10 a4=%11 a5=%12  b0=%13 b1=%14 b2=%15 b3=%16 b4=%17 b5=%18
			asm volatile(
				"sub.cc.u64      %0,  %7,  %13;\n\t"
				"subc.cc.u64     %1,  %8,  %14;\n\t"
				"subc.cc.u64     %2,  %9,  %15;\n\t"
				"subc.cc.u64     %3,  %10, %16;\n\t"
				"subc.cc.u64     %4,  %11, %17;\n\t"
				"subc.cc.u64     %5,  %12, %18;\n\t"
				"addc.u32        %6,  0,   0;\n\t"   // %6 = !borrow
				"xor.b32         %6,  %6,  1;\n\t"   // %6 = borrow
				: "=&l"(r0),
				"=&l"(r1),
				"=&l"(r2),
				"=&l"(r3),
				"=&l"(r4),
				"=&l"(r5),
				"=&r"(borrow)
				: "l"(a[0]),
				"l"(a[1]),
				"l"(a[2]),
				"l"(a[3]),
				"l"(a[4]),
				"l"(a[5]),
				"l"(b[0]),
				"l"(b[1]),
				"l"(b[2]),
				"l"(b[3]),
				"l"(b[4]),
				"l"(b[5])
				);

			limb_t mask_, carry_, borrow_;
			limb_t ret[n];
			{
				size_t i;
				for (borrow_ = 0, i = 0; i < n; i++) {
					wide_limb_t wide;
					cc_sub_borrow_wide_limb_limb(&wide, a[i], b[i], borrow_);
					ret[i] = wide.lo;
					borrow_ = wide.hi & 1;
				}
				BLST_ASSERT(borrow_ == borrow);
				BLST_ASSERT(ret[0] == r0);
				BLST_ASSERT(ret[1] == r1);
				BLST_ASSERT(ret[2] == r2);
				BLST_ASSERT(ret[3] == r3);
				BLST_ASSERT(ret[4] == r4);
				BLST_ASSERT(ret[5] == r5);
			}

			// ---- Pass 2: if (borrow) r += p ; predicated, carry-chained; final add need not set CC ----
			asm volatile(
				".reg .pred pb;\n\t"
				"setp.ne.u32     pb, %12, 0;\n\t"

				"@pb add.cc.u64  %0,  %0,  %6;\n\t"
				"@pb addc.cc.u64 %1,  %1,  %7;\n\t"
				"@pb addc.cc.u64 %2,  %2,  %8;\n\t"
				"@pb addc.cc.u64 %3,  %3,  %9;\n\t"
				"@pb addc.cc.u64 %4,  %4,  %10;\n\t"
				"@pb addc.u64     %5,  %5,  %11;\n\t"
				: "+l"(r0),
				"+l"(r1),
				"+l"(r2),
				"+l"(r3),
				"+l"(r4),
				"+l"(r5)
				: "l"(p[0]),
				"l"(p[1]),
				"l"(p[2]),
				"l"(p[3]),
				"l"(p[4]),
				"l"(p[5]),
				"r"(borrow)
				);

			{
				mask_ = 0 - borrow_;
				BLST_launder(mask_);
				limb_t i;
				for (carry_ = 0, i = 0; i < n; i++) {
					wide_limb_t wide;
					cc_add_carry_wide_limb_limb(&wide, ret[i], p[i] & mask_, carry_);
					ret[i] = wide.lo;
					carry_ = wide.hi;
				}

				BLST_ASSERT(borrow_ == borrow);
				BLST_ASSERT(ret[0] == r0);
				BLST_ASSERT(ret[1] == r1);
				BLST_ASSERT(ret[2] == r2);
				BLST_ASSERT(ret[3] == r3);
				BLST_ASSERT(ret[4] == r4);
				BLST_ASSERT(ret[5] == r5);
			}

			r[0] = r0;
			r[1] = r1;
			r[2] = r2;
			r[3] = r3;
			r[4] = r4;
			r[5] = r5;
		}
		else if constexpr (n == 4)
		{
			uint64_t r0;
			uint64_t r1;
			uint64_t r2;
			uint64_t r3;

			uint32_t borrow;

			// ---- Pass 1: r = a - b ; get final borrow from CC ----
			asm volatile(
				"sub.cc.u64      %0,  %5,  %9;\n\t"
				"subc.cc.u64     %1,  %6,  %10;\n\t"
				"subc.cc.u64     %2,  %7,  %11;\n\t"
				"subc.cc.u64     %3,  %8,  %12;\n\t"
				"addc.u32        %4,  0,   0;\n\t"   // %4 = !borrow
				"xor.b32         %4,  %4,  1;\n\t"   // %4 = borrow
				: "=&l"(r0),
				"=&l"(r1),
				"=&l"(r2),
				"=&l"(r3),
				"=&r"(borrow)
				: "l"(a[0]),
				"l"(a[1]),
				"l"(a[2]),
				"l"(a[3]),
				"l"(b[0]),
				"l"(b[1]),
				"l"(b[2]),
				"l"(b[3])
				);

			// ---- Pass 2: if (borrow) r += p ; predicated, carry-chained ----
			asm volatile(
				".reg .pred pb;\n\t"
				"setp.ne.u32     pb, %8, 0;\n\t"

				"@pb add.cc.u64  %0,  %0,  %4;\n\t"
				"@pb addc.cc.u64 %1,  %1,  %5;\n\t"
				"@pb addc.cc.u64 %2,  %2,  %6;\n\t"
				"@pb addc.u64     %3,  %3,  %7;\n\t"
				: "+l"(r0),
				"+l"(r1),
				"+l"(r2),
				"+l"(r3)
				: "l"(p[0]),
				"l"(p[1]),
				"l"(p[2]),
				"l"(p[3]),
				"r"(borrow)
				);

			r[0] = r0;
			r[1] = r1;
			r[2] = r2;
			r[3] = r3;
		}
		else
		{
			// static assert to trigger a compile-time error for unsupported Field sizes
			static_assert(n == 6 || n == 4, "ptx_sub_mod_u64 only supports 256 or 384-bit fields");
		}


		//#ifndef NDEBUG
		//	uint64_t tt[n];
		//	blst::cc_sub_mod_n<n>(tt, a, b, Field::modulus().data);
		//	for (size_t i = 0; i < n; i++) {
		//		BLST_ASSERT(r[i] == tt[i]);
		//	}
		//#endif
	}
	//
	//
	//template<typename Field>
	//BLST_DEVICE BLST_CPP_INLINE
	//void ptx_sub_mod_u64(
	//	uint64_t* BLST_restrict r,
	//	const uint64_t* BLST_restrict a,
	//	const uint64_t* BLST_restrict b)
	//{
	//	constexpr auto n = sizeof(Field) / sizeof(uint64_t);
	//	constexpr auto BLST_restrict p = Field::modulus().data;
	//
	//	if constexpr (n == 6)
	//	{
	//		// x = a - b
	//		uint64_t t0;
	//		uint64_t t1;
	//		uint64_t t2;
	//		uint64_t t3;
	//		uint64_t t4;
	//		uint64_t t5;
	//		asm volatile(
	//			"sub.cc.u64  %0,  %6,  %12;\n\t"
	//			"subc.cc.u64 %1,  %7,  %13;\n\t"
	//			"subc.cc.u64 %2,  %8,  %14;\n\t"
	//			"subc.cc.u64 %3,  %9,  %15;\n\t"
	//			"subc.cc.u64 %4,  %10, %16;\n\t"
	//			"subc.u64    %5,  %11, %17;\n\t"
	//			: "=&l"(t0),
	//			"=&l"(t1),
	//			"=&l"(t2),
	//			"=&l"(t3),
	//			"=&l"(t4),
	//			"=&l"(t5)
	//			: "l"(a[0]), "l"(a[1]), "l"(a[2]), "l"(a[3]), "l"(a[4]), "l"(a[5]),
	//			"l"(b[0]), "l"(b[1]), "l"(b[2]), "l"(b[3]), "l"(b[4]), "l"(b[5])
	//			);
	//
	//		// s = x + p
	//		uint64_t s0;
	//		uint64_t s1;
	//		uint64_t s2;
	//		uint64_t s3;
	//		uint64_t s4;
	//		uint64_t s5;
	//		asm volatile(
	//			"add.cc.u64  %0,  %6,  %12;\n\t"
	//			"addc.cc.u64 %1,  %7,  %13;\n\t"
	//			"addc.cc.u64 %2,  %8,  %14;\n\t"
	//			"addc.cc.u64 %3,  %9,  %15;\n\t"
	//			"addc.cc.u64 %4,  %10, %16;\n\t"
	//			"addc.u64    %5,  %11, %17;\n\t"
	//			: "=&l"(s0),
	//			"=&l"(s1),
	//			"=&l"(s2),
	//			"=&l"(s3),
	//			"=&l"(s4),
	//			"=&l"(s5)
	//			: "l"(t0), "l"(t1), "l"(t2), "l"(t3), "l"(t4), "l"(t5),
	//			"l"(p[0]), "l"(p[1]), "l"(p[2]), "l"(p[3]), "l"(p[4]), "l"(p[5])
	//			);
	//
	//		// need_add = (a < b)
	//		uint32_t need_add;
	//		asm volatile(
	//			".reg .pred plt;\n\t"
	//			".reg .pred peq;\n\t"
	//			".reg .pred peq_next;\n\t"
	//			".reg .pred pt;\n\t"
	//
	//			// Compare from most significant limb: a5..a0 vs b5..b0
	//			"setp.lt.u64 plt, %6,  %12;\n\t"
	//			"setp.eq.u64 peq, %6,  %12;\n\t"
	//
	//			"setp.lt.u64 pt,  %5,  %11;\n\t"
	//			"and.pred    pt,  pt,  peq;\n\t"
	//			"or.pred     plt, plt, pt;\n\t"
	//			"setp.eq.u64 peq_next, %5, %11;\n\t"
	//			"and.pred    peq, peq, peq_next;\n\t"
	//
	//			"setp.lt.u64 pt,  %4,  %10;\n\t"
	//			"and.pred    pt,  pt,  peq;\n\t"
	//			"or.pred     plt, plt, pt;\n\t"
	//			"setp.eq.u64 peq_next, %4, %10;\n\t"
	//			"and.pred    peq, peq, peq_next;\n\t"
	//
	//			"setp.lt.u64 pt,  %3,  %9;\n\t"
	//			"and.pred    pt,  pt,  peq;\n\t"
	//			"or.pred     plt, plt, pt;\n\t"
	//			"setp.eq.u64 peq_next, %3, %9;\n\t"
	//			"and.pred    peq, peq, peq_next;\n\t"
	//
	//			"setp.lt.u64 pt,  %2,  %8;\n\t"
	//			"and.pred    pt,  pt,  peq;\n\t"
	//			"or.pred     plt, plt, pt;\n\t"
	//			"setp.eq.u64 peq_next, %2, %8;\n\t"
	//			"and.pred    peq, peq, peq_next;\n\t"
	//
	//			"setp.lt.u64 pt,  %1,  %7;\n\t"
	//			"and.pred    pt,  pt,  peq;\n\t"
	//			"or.pred     plt, plt, pt;\n\t"
	//
	//			"selp.u32    %0, 1, 0, plt;\n\t"
	//			: "=&r"(need_add)
	//			: "l"(a[0]), "l"(a[1]), "l"(a[2]), "l"(a[3]), "l"(a[4]), "l"(a[5]),
	//			"l"(b[0]), "l"(b[1]), "l"(b[2]), "l"(b[3]), "l"(b[4]), "l"(b[5])
	//			);
	//
	//		// r = need_add ? s : t
	//		uint64_t o0;
	//		uint64_t o1;
	//		uint64_t o2;
	//		uint64_t o3;
	//		uint64_t o4;
	//		uint64_t o5;
	//
	//		asm volatile(
	//			".reg .pred psel;\n\t"
	//			"setp.ne.u32 psel, %12, 0;\n\t"
	//			"selp.u64 %0, %6,  %13, psel;\n\t"
	//			"selp.u64 %1, %7,  %14, psel;\n\t"
	//			"selp.u64 %2, %8,  %15, psel;\n\t"
	//			"selp.u64 %3, %9,  %16, psel;\n\t"
	//			"selp.u64 %4, %10, %17, psel;\n\t"
	//			"selp.u64 %5, %11, %18, psel;\n\t"
	//			: "=&l"(o0),
	//			"=&l"(o1),
	//			"=&l"(o2),
	//			"=&l"(o3),
	//			"=&l"(o4),
	//			"=&l"(o5)
	//			: "l"(s0), "l"(s1), "l"(s2), "l"(s3), "l"(s4), "l"(s5),
	//			"r"(need_add),
	//			"l"(t0), "l"(t1), "l"(t2), "l"(t3), "l"(t4), "l"(t5)
	//			);
	//
	//		r[0] = o0;
	//		r[1] = o1;
	//		r[2] = o2;
	//		r[3] = o3;
	//		r[4] = o4;
	//		r[5] = o5;
	//	}
	//	else if constexpr (n == 4)
	//	{
	//
	//		// x = a - b
	//		uint64_t t0;
	//		uint64_t t1;
	//		uint64_t t2;
	//		uint64_t t3;
	//
	//		asm volatile(
	//			"sub.cc.u64  %0,  %4,  %8;\n\t"
	//			"subc.cc.u64 %1,  %5,  %9;\n\t"
	//			"subc.cc.u64 %2,  %6,  %10;\n\t"
	//			"subc.u64    %3,  %7,  %11;\n\t"
	//			: "=&l"(t0),
	//			"=&l"(t1),
	//			"=&l"(t2),
	//			"=&l"(t3)
	//			: "l"(a[0]), "l"(a[1]), "l"(a[2]), "l"(a[3]),
	//			"l"(b[0]), "l"(b[1]), "l"(b[2]), "l"(b[3])
	//			);
	//
	//		// s = x + p
	//		uint64_t s0;
	//		uint64_t s1;
	//		uint64_t s2;
	//		uint64_t s3;
	//
	//		asm volatile(
	//			"add.cc.u64  %0,  %4,  %8;\n\t"
	//			"addc.cc.u64 %1,  %5,  %9;\n\t"
	//			"addc.cc.u64 %2,  %6,  %10;\n\t"
	//			"addc.u64    %3,  %7,  %11;\n\t"
	//			: "=&l"(s0),
	//			"=&l"(s1),
	//			"=&l"(s2),
	//			"=&l"(s3)
	//			: "l"(t0), "l"(t1), "l"(t2), "l"(t3),
	//			"l"(p[0]), "l"(p[1]), "l"(p[2]), "l"(p[3])
	//			);
	//
	//		uint32_t need_add;
	//
	//		// need_add = (a < b)  (lexicographic compare, MS limb first)
	//		asm volatile(
	//			".reg .pred plt;\n\t"
	//			".reg .pred peq;\n\t"
	//			".reg .pred peq_next;\n\t"
	//			".reg .pred ptmp;\n\t"
	//
	//			// a3 ? b3
	//			"setp.lt.u64 plt, %4, %8;\n\t"
	//			"setp.eq.u64 peq, %4, %8;\n\t"
	//
	//			// a2 ? b2 (only matters if a3==b3)
	//			"setp.lt.u64 ptmp, %3, %7;\n\t"
	//			"and.pred    ptmp, ptmp, peq;\n\t"
	//			"or.pred     plt,  plt,  ptmp;\n\t"
	//			"setp.eq.u64 peq_next, %3, %7;\n\t"
	//			"and.pred    peq, peq, peq_next;\n\t"
	//
	//			// a1 ? b1
	//			"setp.lt.u64 ptmp, %2, %6;\n\t"
	//			"and.pred    ptmp, ptmp, peq;\n\t"
	//			"or.pred     plt,  plt,  ptmp;\n\t"
	//			"setp.eq.u64 peq_next, %2, %6;\n\t"
	//			"and.pred    peq, peq, peq_next;\n\t"
	//
	//			// a0 ? b0
	//			"setp.lt.u64 ptmp, %1, %5;\n\t"
	//			"and.pred    ptmp, ptmp, peq;\n\t"
	//			"or.pred     plt,  plt,  ptmp;\n\t"
	//
	//			// need_add = plt ? 1 : 0
	//			"selp.u32    %0, 1, 0, plt;\n\t"
	//			: "=&r"(need_add)
	//			: "l"(a[0]), "l"(a[1]), "l"(a[2]), "l"(a[3]),
	//			"l"(b[0]), "l"(b[1]), "l"(b[2]), "l"(b[3])
	//			);
	//
	//		// r = need_add ? s : t
	//		uint64_t o0;
	//		uint64_t o1;
	//		uint64_t o2;
	//		uint64_t o3;
	//
	//		asm volatile(
	//			".reg .pred psel;\n\t"
	//			"setp.ne.u32 psel, %8, 0;\n\t"
	//			"selp.u64 %0, %4, %9,  psel;\n\t"
	//			"selp.u64 %1, %5, %10, psel;\n\t"
	//			"selp.u64 %2, %6, %11, psel;\n\t"
	//			"selp.u64 %3, %7, %12, psel;\n\t"
	//			: "=&l"(o0),
	//			"=&l"(o1),
	//			"=&l"(o2),
	//			"=&l"(o3)
	//			: "l"(s0), "l"(s1), "l"(s2), "l"(s3),
	//			"r"(need_add),
	//			"l"(t0), "l"(t1), "l"(t2), "l"(t3)
	//			);
	//
	//		r[0] = o0;
	//		r[1] = o1;
	//		r[2] = o2;
	//		r[3] = o3;
	//	}
	//	else
	//	{
	//		static_assert(n == 6 || n == 4, "ptx_sub_mod_u64 only supports 256 or 384-bit fields");
	//	}
	//
	//#ifndef NDEBUG
	//	uint64_t tt[n];
	//	blst::cc_sub_mod_n<n>(tt, a, b, Field::modulus().data);
	//	for (size_t i = 0; i < n; i++) {
	//		BLST_ASSERT(r[i] == tt[i]);
	//	}
	//#endif
	//}
}