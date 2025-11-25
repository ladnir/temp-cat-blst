#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>


namespace blst {


    //void shuffle_down(
    //    uint32_t out[],
    //    const uint32_t in[],
    //    int n,
    //    int shift);

    //void shuffle_up(
    //    uint32_t out[],
    //    const uint32_t in[],
    //    int n,
    //    int shift);

  //  template<int m, typename Limb = uint64_t>
  //  void compress5to2(
  //      const Limb* x0,
  //      const Limb* x1,
  //      const Limb* x2,
  //      const Limb* x3,
  //      const Limb* x4,
  //      Limb* s,
  //      Limb* c
  //  ) {
  //      static_assert(std::is_unsigned<Limb>::value, "Limb must be unsigned");
  //      constexpr int W = int(sizeof(Limb) * 8);

  //      Limb spill1 = 0;  // for (n1<<2): two-bit cross-limb spill (values 0..3 in the low bits)
  //      Limb spill0 = 0;  // for (n0<<1): one-bit cross-limb spill   (values 0/1  in the low bit)
		//Limb spillc = 0;  // for (c<<1): one-bit cross-limb spill   (values 0/1  in the low bit)

  //      for (std::size_t i = 0; i < m; i++) {
  //          Limb a = x0[i];
  //          Limb b = x1[i];
  //          Limb c0 = x2[i];
  //          Limb d = x3[i];
  //          Limb e = x4[i];

  //          // ---- Round 1: 5 -> 3 (X0..X4 = S0 + 2*n0 + 4*n1), all per-bit within the limb ----

  //          // (a,b,c0) -> sA + 2*cA
  //          Limb sA = a ^ b ^ c0;
  //          Limb cA = (a & b) | (b & c0) | (a & c0);

  //          // (d,e,0)  -> sB + 2*cB
  //          Limb sB = d ^ e;
  //          Limb cB = d & e;

  //          // sA + sB -> S0 + 2*cAB
  //          Limb S0 = sA ^ sB;
  //          Limb cAB = sA & sB;

  //          // N = cA + cB + cAB in {0,1,2,3}; encode as n0 (bit0), n1 (bit1)
  //          Limb n0 = (cA ^ cB) ^ cAB;                           // parity of N
  //          Limb n1 = (cA & cB) ^ (cA & cAB) ^ (cB & cAB);       // bit 1 of N

  //          // ---- Round 2: 3 -> 2 CSA on S0, (n0<<1), (n1<<2) with cross-limb stitching ----

  //          // Prepare aligned inputs with prior spills OR'ed into the LSBs.
  //          Limb A = S0;
  //          Limb B = (n0 << 1) | spill0;
  //          Limb C = (n1 << 2) | spill1;

  //          // Full-adder (bitwise): S = A^B^C,  K = (A&B) | (A&C) | (B&C)
  //          Limb S = (A ^ B) ^ C;
  //          Limb K = (A & B) | (A & C) | (B & C);

  //          s[i] = S;
  //          c[i] = (K << 1) | spillc;
  //          spillc = K >> (W - 1);

  //          // Compute next spills from the left shifts we just did
  //          spill0 = (n0 >> (W - 1));
  //          spill1 = (n1 >> (W - 2));  // two-bit value (0..3) to drop into next limb's low bits
  //      }

  //      // Note: If you need the exact total (not mod 2^n), you can recover the two
  //      // overflow bits as:
  //      //   top_bit_n       = (c[m-1] >> (W-1)) ^ (spill0 & 1) ^ (spill1 & 1);
  //      //   top_bit_n_plus1 = (spill1 >> 1) ^ (((c[m-1] >> (W-1)) & (spill0 & 1)) | ((c[m-1] >> (W-1)) & (spill1 & 1)) | ((spill0 & 1) & (spill1 & 1)));
  //      // or expose them as optional out-params. They’re omitted here by request.
  //  }


}
