#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>


namespace blst {


    void shuffle_down(
        uint32_t out[],
        const uint32_t in[],
        int n,
        int shift);

    void shuffle_up(
        uint32_t out[],
        const uint32_t in[],
        int n,
        int shift);

    // 5->2 compressor over m limbs of type Limb (default 64-bit).
    // Inputs : cin  (a one limb carry in. can have value [0,1,2,3])
    // Inputs : x0, x1, x2, x3, x4  (each m-limb n-bit vector)
	// Outputs: s, c  
    // 
    // s + c = x0 + x1 + x2 + x3 + x4 + cin
    //
    // Internals:
    //   Round 1 (5->3, column-local):  X0+...+X4 = S0 + 2*n0 + 4*n1
    //   Round 2 (3->2 CSA):            (S0), (n0<<1), (n1<<2)  ->  s + c
    //
    // No '+' is used inside a limb; everything is XOR/AND/OR/<< with cross-limb “spill” stitching.
    // The two possible overflow bits (at 2^n and 2^(n+1)) are dropped by design in this signature.
    template<int m, typename Limb = uint64_t>
    void compress5to2(
        Limb cin,
        const Limb* x0,
        const Limb* x1,
        const Limb* x2,
        const Limb* x3,
        const Limb* x4,
        Limb* s,
        Limb* c
    ) {
        static_assert(std::is_unsigned<Limb>::value, "Limb must be unsigned");
        constexpr int W = int(sizeof(Limb) * 8);

		BLST_ASSERT(cin < 4);  // cin must be in {0,1,2,3}
        // -------- Round 1: 5->3 local; compute S0, n0, n1 and pack outgoing spills --------
        Limb S0[m];
        Limb n0[m];
        Limb n1[m];

        // pack3[i] holds 3 bits: [bit2=b2 | bit1=b1 | bit0=b0]
        // b0 = (n0 >> (W-1))       (spill for <<1)
        // b1 = (n1 >> (W-2)) & 1   (LSB spill for <<2)
        // b2 = (n1 >> (W-1)) & 1   (MSB spill for <<2)
        uint32_t pack3[m];

        for (int i = 0; i < m; i++) {
            Limb a = x0[i];
            Limb b = x1[i];
            Limb c0 = x2[i];
            Limb d = x3[i];
            Limb e = x4[i];

            // (a,b,c0) -> sA + 2*cA
            Limb sA = a ^ b ^ c0;
            Limb cA = (a & b) | (b & c0) | (a & c0);

            // (d,e,0) -> sB + 2*cB
            Limb sB = d ^ e;
            Limb cB = d & e;

            // sA + sB -> S0 + 2*cAB
            Limb s0 = sA ^ sB;
            Limb cAB = sA & sB;

            // N = cA + cB + cAB in {0,1,2,3}
            Limb n0_i = (cA ^ cB) ^ cAB;                       // bit0(N)
            Limb n1_i = (cA & cB) ^ (cA & cAB) ^ (cB & cAB);   // bit1(N)

            S0[i] = s0;
            n0[i] = n0_i;
            n1[i] = n1_i;

            uint32_t b0 = uint32_t(n0_i >> (W - 1));
            uint32_t b1 = uint32_t(n1_i >> (W - 2));
            pack3[i] = (b1 << 1) | b0;
        }

        // Communicate the 3 spill bits up by one limb (shift=1)
        uint32_t in3[m];
        shuffle_up(in3, pack3, m, /*shift=*/1);

        // -------- Form CSA inputs with received spills and do 3:2 CSA --------
        // A = S0
        // B = (n0<<1) | in_b0
        // C = (n1<<2) | (in_b1) | (in_b2<<1)
        // Then S = A^B^C; K = (A&B)|(A&C)|(B&C)
        Limb K_row[m];
        for (int i = 0; i < m; i++) {
            uint32_t in = i ? in3[i] : (cin << 1);
            uint32_t b0 = in & 1;
            uint32_t b1 = in >> 1;
            std::cout << "i " << i << " cin " << b0 << b1 << std::endl;

            Limb A = S0[i];
            Limb B = (n0[i] << 1) | Limb(b0);
            Limb C = (n1[i] << 2) | Limb(b1);

            Limb S = (A ^ B) ^ C;
            Limb K = (A & B) | (A & C) | (B & C);

            s[i] = S;      // final sum row
            K_row[i] = K;      // unshifted carry row (needs <<1 across limbs to become 'c')
        
            pack3[i] = K_row[i] >> (W - 1);  // single spill bit
        }
        shuffle_up(in3, pack3, m, /*shift=*/1);

        for (int i = 0; i < m; i++) {
			auto in = i ? in3[i] : 0;
            c[i] = (K_row[i] << 1) | in;
        }

        {
            Limb exp[m+1];
            exp[0] = cin;
            Limb cc = 0;
            for (u64 i = 0; i < m; ++i)
            {
                wide_limb_t w{ 0,0 };
                cc_add_wide_wide_limb(&w, w, exp[i]);
                cc_add_wide_wide_limb(&w, w, x0[i]);
                cc_add_wide_wide_limb(&w, w, x1[i]);
                cc_add_wide_wide_limb(&w, w, x2[i]);
                cc_add_wide_wide_limb(&w, w, x3[i]);
                cc_add_wide_wide_limb(&w, w, x4[i]);

                exp[i] = w.lo;
                exp[i+1] = w.hi;

                
                //auto act = s[i] + c[i];
                cc_add_wide_limb_limb(&w, s[i], c[i]);
                cc_add_wide_wide_limb(&w, w, cc);
                cc = w.hi;
                if (exp[i] != w.lo)
                {
                    std::cout << "bad at " << i << std::endl;
                }
            }
        }
    }

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
