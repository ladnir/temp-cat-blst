
#pragma once

#include "blst/scalar/curve/point.h"
#include "blst/scalar/fields/hash_to_field.h"


namespace blst::scalar
{





	template<typename Field>
	BLST_HD BLST_INLINE void map_field_times_Zz(Field map[], const Field isogeny_map[],
		const Field Zz_powers[], size_t n)
	{
		while (n--)
			Field::mul(map[n], isogeny_map[n], Zz_powers[n]);
	}


	template<typename Field>
	BLST_HD BLST_INLINE void map_field(Field& acc, const Field& x, const Field map[], size_t n)
	{
		while (n--) {
			Field::mul(acc, acc, x);
			Field::add(acc, acc, map[n]);
		}
	}


	template<typename Group>
	BLST_HD BLST_INLINE void isogeny_map_to(Group& out, const Group& p);


	template<>
	BLST_HD BLST_INLINE void isogeny_map_to<G2Point>(G2Point& out, const G2Point& p)
	{
		/*
		 * x = x_num / x_den, where
		 * x_num = k_(1,3) * x'^3 + k_(1,2) * x'^2 + k_(1,1) * x' + k_(1,0)
		 * ...
		 */
		const Fp2 isogeny_map_x_num[] = {    /* (k_(1,*)<<384) % P   */
		 {Fp{ TO_LIMB_T(0x47f671c71ce05e62), TO_LIMB_T(0x06dd57071206393e),
			TO_LIMB_T(0x7c80cd2af3fd71a2), TO_LIMB_T(0x048103ea9e6cd062),
			TO_LIMB_T(0xc54516acc8d037f6), TO_LIMB_T(0x13808f550920ea41) },
		  Fp{ TO_LIMB_T(0x47f671c71ce05e62), TO_LIMB_T(0x06dd57071206393e),
			TO_LIMB_T(0x7c80cd2af3fd71a2), TO_LIMB_T(0x048103ea9e6cd062),
			TO_LIMB_T(0xc54516acc8d037f6), TO_LIMB_T(0x13808f550920ea41) }},
		 {Fp{ 0 },
		  Fp{ TO_LIMB_T(0x5fe55555554c71d0), TO_LIMB_T(0x873fffdd236aaaa3),
			TO_LIMB_T(0x6a6b4619b26ef918), TO_LIMB_T(0x21c2888408874945),
			TO_LIMB_T(0x2836cda7028cabc5), TO_LIMB_T(0x0ac73310a7fd5abd) }},
		 {Fp{ TO_LIMB_T(0x0a0c5555555971c3), TO_LIMB_T(0xdb0c00101f9eaaae),
			TO_LIMB_T(0xb1fb2f941d797997), TO_LIMB_T(0xd3960742ef416e1c),
			TO_LIMB_T(0xb70040e2c20556f4), TO_LIMB_T(0x149d7861e581393b) },
		  Fp{ TO_LIMB_T(0xaff2aaaaaaa638e8), TO_LIMB_T(0x439fffee91b55551),
			TO_LIMB_T(0xb535a30cd9377c8c), TO_LIMB_T(0x90e144420443a4a2),
			TO_LIMB_T(0x941b66d3814655e2), TO_LIMB_T(0x0563998853fead5e) }},
		 {Fp{ TO_LIMB_T(0x40aac71c71c725ed), TO_LIMB_T(0x190955557a84e38e),
			TO_LIMB_T(0xd817050a8f41abc3), TO_LIMB_T(0xd86485d4c87f6fb1),
			TO_LIMB_T(0x696eb479f885d059), TO_LIMB_T(0x198e1a74328002d2) },
		  Fp{ 0 }}
		};
		/* ...
		 * x_den = x'^2 + k_(2,1) * x' + k_(2,0)
		 */
		const Fp2 isogeny_map_x_den[] = {    /* (k_(2,*)<<384) % P   */
		 {Fp{ 0 },
		  Fp{ TO_LIMB_T(0x1f3affffff13ab97), TO_LIMB_T(0xf25bfc611da3ff3e),
			TO_LIMB_T(0xca3757cb3819b208), TO_LIMB_T(0x3e6427366f8cec18),
			TO_LIMB_T(0x03977bc86095b089), TO_LIMB_T(0x04f69db13f39a952) }},
		 {Fp{ TO_LIMB_T(0x447600000027552e), TO_LIMB_T(0xdcb8009a43480020),
			TO_LIMB_T(0x6f7ee9ce4a6e8b59), TO_LIMB_T(0xb10330b7c0a95bc6),
			TO_LIMB_T(0x6140b1fcfb1e54b7), TO_LIMB_T(0x0381be097f0bb4e1) },
		  Fp{ TO_LIMB_T(0x7588ffffffd8557d), TO_LIMB_T(0x41f3ff646e0bffdf),
			TO_LIMB_T(0xf7b1e8d2ac426aca), TO_LIMB_T(0xb3741acd32dbb6f8),
			TO_LIMB_T(0xe9daf5b9482d581f), TO_LIMB_T(0x167f53e0ba7431b8) }}
		};
		/*
		 * y = y' * y_num / y_den, where
		 * y_num = k_(3,3) * x'^3 + k_(3,2) * x'^2 + k_(3,1) * x' + k_(3,0)
		 * ...
		 */
		const Fp2 isogeny_map_y_num[] = {    /* (k_(3,*)<<384) % P   */
		 {Fp{ TO_LIMB_T(0x96d8f684bdfc77be), TO_LIMB_T(0xb530e4f43b66d0e2),
			TO_LIMB_T(0x184a88ff379652fd), TO_LIMB_T(0x57cb23ecfae804e1),
			TO_LIMB_T(0x0fd2e39eada3eba9), TO_LIMB_T(0x08c8055e31c5d5c3) },
		  Fp{ TO_LIMB_T(0x96d8f684bdfc77be), TO_LIMB_T(0xb530e4f43b66d0e2),
			TO_LIMB_T(0x184a88ff379652fd), TO_LIMB_T(0x57cb23ecfae804e1),
			TO_LIMB_T(0x0fd2e39eada3eba9), TO_LIMB_T(0x08c8055e31c5d5c3) }},
		 {Fp{ 0 },
		  Fp{ TO_LIMB_T(0xbf0a71c71c91b406), TO_LIMB_T(0x4d6d55d28b7638fd),
			TO_LIMB_T(0x9d82f98e5f205aee), TO_LIMB_T(0xa27aa27b1d1a18d5),
			TO_LIMB_T(0x02c3b2b2d2938e86), TO_LIMB_T(0x0c7d13420b09807f) }},
		 {Fp{ TO_LIMB_T(0xd7f9555555531c74), TO_LIMB_T(0x21cffff748daaaa8),
			TO_LIMB_T(0x5a9ad1866c9bbe46), TO_LIMB_T(0x4870a2210221d251),
			TO_LIMB_T(0x4a0db369c0a32af1), TO_LIMB_T(0x02b1ccc429ff56af) },
		  Fp{ TO_LIMB_T(0xe205aaaaaaac8e37), TO_LIMB_T(0xfcdc000768795556),
			TO_LIMB_T(0x0c96011a8a1537dd), TO_LIMB_T(0x1c06a963f163406e),
			TO_LIMB_T(0x010df44c82a881e6), TO_LIMB_T(0x174f45260f808feb) }},
		 {Fp{ TO_LIMB_T(0xa470bda12f67f35c), TO_LIMB_T(0xc0fe38e23327b425),
			TO_LIMB_T(0xc9d3d0f2c6f0678d), TO_LIMB_T(0x1c55c9935b5a982e),
			TO_LIMB_T(0x27f6c0e2f0746764), TO_LIMB_T(0x117c5e6e28aa9054) },
		  Fp{ 0 }}
		};
		/* ...
		 * y_den = x'^3 + k_(4,2) * x'^2 + k_(4,1) * x' + k_(4,0)
		 */
		const Fp2 isogeny_map_y_den[] = {    /* (k_(4,*)<<384) % P   */
		 {Fp{ TO_LIMB_T(0x0162fffffa765adf), TO_LIMB_T(0x8f7bea480083fb75),
			TO_LIMB_T(0x561b3c2259e93611), TO_LIMB_T(0x11e19fc1a9c875d5),
			TO_LIMB_T(0xca713efc00367660), TO_LIMB_T(0x03c6a03d41da1151) },
		  Fp{ TO_LIMB_T(0x0162fffffa765adf), TO_LIMB_T(0x8f7bea480083fb75),
			TO_LIMB_T(0x561b3c2259e93611), TO_LIMB_T(0x11e19fc1a9c875d5),
			TO_LIMB_T(0xca713efc00367660), TO_LIMB_T(0x03c6a03d41da1151) }},
		 {Fp{ 0 },
		  Fp{ TO_LIMB_T(0x5db0fffffd3b02c5), TO_LIMB_T(0xd713f52358ebfdba),
			TO_LIMB_T(0x5ea60761a84d161a), TO_LIMB_T(0xbb2c75a34ea6c44a),
			TO_LIMB_T(0x0ac6735921c1119b), TO_LIMB_T(0x0ee3d913bdacfbf6) }},
		 {Fp{ TO_LIMB_T(0x66b10000003affc5), TO_LIMB_T(0xcb1400e764ec0030),
			TO_LIMB_T(0xa73e5eb56fa5d106), TO_LIMB_T(0x8984c913a0fe09a9),
			TO_LIMB_T(0x11e10afb78ad7f13), TO_LIMB_T(0x05429d0e3e918f52) },
		  Fp{ TO_LIMB_T(0x534dffffffc4aae6), TO_LIMB_T(0x5397ff174c67ffcf),
			TO_LIMB_T(0xbff273eb870b251d), TO_LIMB_T(0xdaf2827152870915),
			TO_LIMB_T(0x393a9cbaca9e2dc3), TO_LIMB_T(0x14be74dbfaee5748) }}
		};
		Fp2 Zz_powers[3], map[3], xn, xd, yn, yd;

		/* lay down Z^2 powers in descending order                          */
		Fp2::sqr(Zz_powers[2], p.z);                       /* ZZ^1          */
		Fp2::sqr(Zz_powers[1], Zz_powers[2]);               /* ZZ^2  1+1     */
		Fp2::mul(Zz_powers[0], Zz_powers[2], Zz_powers[1]); /* ZZ^3  2+1     */

		map_field_times_Zz<Fp2>(map, isogeny_map_x_num, Zz_powers, 3);
		Fp2::mul(xn, p.x, isogeny_map_x_num[3]);
		Fp2::add(xn, xn, map[2]);
		map_field<Fp2>(xn, p.x, map, 2);

		map_field_times_Zz<Fp2>(map, isogeny_map_x_den, Zz_powers + 1, 2);
		Fp2::add(xd, p.x, map[1]);
		map_field<Fp2>(xd, p.x, map, 1);
		Fp2::mul(xd, xd, Zz_powers[2]);      /* xd *= Z^2                    */

		map_field_times_Zz<Fp2>(map, isogeny_map_y_num, Zz_powers, 3);
		Fp2::mul(yn, p.x, isogeny_map_y_num[3]);
		Fp2::add(yn, yn, map[2]);
		map_field<Fp2>(yn, p.x, map, 2);
		Fp2::mul(yn, yn, p.y);              /* yn *= Y                      */

		map_field_times_Zz<Fp2>(map, isogeny_map_y_den, Zz_powers, 3);
		Fp2::add(yd, p.x, map[2]);
		map_field<Fp2>(yd, p.x, map, 2);
		Fp2::mul(Zz_powers[2], Zz_powers[2], p.z);
		Fp2::mul(yd, yd, Zz_powers[2]);      /* yd *= Z^3                    */

		/* convert (xn, xd, yn, yd) to Jacobian coordinates                 */
		Fp2::mul(out.z, xd, yd);            /* Z = xd * yd                  */
		Fp2::mul(out.x, xn, yd);
		Fp2::mul(out.x, out.x, out.z);    /* X = xn * xd * yd^2           */
		Fp2::sqr(out.y, out.z);
		Fp2::mul(out.y, out.y, xd);
		Fp2::mul(out.y, out.y, yn);        /* Y = yn * xd^3 * yd^2         */
	}




	template<>
	BLST_HD BLST_INLINE void isogeny_map_to<G1Point>(G1Point& out, const G1Point& p)
	{
		/*
		 * x = x_num / x_den, where
		 * x_num = k_(1,11) * x'^11 + k_(1,10) * x'^10 + k_(1,9) * x'^9 +
		 *         ... + k_(1,0)
		 * ...
		 */
		const Fp isogeny_map_x_num[] = { /*  (k_(1,*)<<384) % P  */
		  { TO_LIMB_T(0x4d18b6f3af00131c), TO_LIMB_T(0x19fa219793fee28c),
			TO_LIMB_T(0x3f2885f1467f19ae), TO_LIMB_T(0x23dcea34f2ffb304),
			TO_LIMB_T(0xd15b58d2ffc00054), TO_LIMB_T(0x0913be200a20bef4)  },
		  { TO_LIMB_T(0x898985385cdbbd8b), TO_LIMB_T(0x3c79e43cc7d966aa),
			TO_LIMB_T(0x1597e193f4cd233a), TO_LIMB_T(0x8637ef1e4d6623ad),
			TO_LIMB_T(0x11b22deed20d827b), TO_LIMB_T(0x07097bc5998784ad)  },
		  { TO_LIMB_T(0xa542583a480b664b), TO_LIMB_T(0xfc7169c026e568c6),
			TO_LIMB_T(0x5ba2ef314ed8b5a6), TO_LIMB_T(0x5b5491c05102f0e7),
			TO_LIMB_T(0xdf6e99707d2a0079), TO_LIMB_T(0x0784151ed7605524)  },
		  { TO_LIMB_T(0x494e212870f72741), TO_LIMB_T(0xab9be52fbda43021),
			TO_LIMB_T(0x26f5577994e34c3d), TO_LIMB_T(0x049dfee82aefbd60),
			TO_LIMB_T(0x65dadd7828505289), TO_LIMB_T(0x0e93d431ea011aeb)  },
		  { TO_LIMB_T(0x90ee774bd6a74d45), TO_LIMB_T(0x7ada1c8a41bfb185),
			TO_LIMB_T(0x0f1a8953b325f464), TO_LIMB_T(0x104c24211be4805c),
			TO_LIMB_T(0x169139d319ea7a8f), TO_LIMB_T(0x09f20ead8e532bf6)  },
		  { TO_LIMB_T(0x6ddd93e2f43626b7), TO_LIMB_T(0xa5482c9aa1ccd7bd),
			TO_LIMB_T(0x143245631883f4bd), TO_LIMB_T(0x2e0a94ccf77ec0db),
			TO_LIMB_T(0xb0282d480e56489f), TO_LIMB_T(0x18f4bfcbb4368929)  },
		  { TO_LIMB_T(0x23c5f0c953402dfd), TO_LIMB_T(0x7a43ff6958ce4fe9),
			TO_LIMB_T(0x2c390d3d2da5df63), TO_LIMB_T(0xd0df5c98e1f9d70f),
			TO_LIMB_T(0xffd89869a572b297), TO_LIMB_T(0x1277ffc72f25e8fe)  },
		  { TO_LIMB_T(0x79f4f0490f06a8a6), TO_LIMB_T(0x85f894a88030fd81),
			TO_LIMB_T(0x12da3054b18b6410), TO_LIMB_T(0xe2a57f6505880d65),
			TO_LIMB_T(0xbba074f260e400f1), TO_LIMB_T(0x08b76279f621d028)  },
		  { TO_LIMB_T(0xe67245ba78d5b00b), TO_LIMB_T(0x8456ba9a1f186475),
			TO_LIMB_T(0x7888bff6e6b33bb4), TO_LIMB_T(0xe21585b9a30f86cb),
			TO_LIMB_T(0x05a69cdcef55feee), TO_LIMB_T(0x09e699dd9adfa5ac)  },
		  { TO_LIMB_T(0x0de5c357bff57107), TO_LIMB_T(0x0a0db4ae6b1a10b2),
			TO_LIMB_T(0xe256bb67b3b3cd8d), TO_LIMB_T(0x8ad456574e9db24f),
			TO_LIMB_T(0x0443915f50fd4179), TO_LIMB_T(0x098c4bf7de8b6375)  },
		  { TO_LIMB_T(0xe6b0617e7dd929c7), TO_LIMB_T(0xfe6e37d442537375),
			TO_LIMB_T(0x1dafdeda137a489e), TO_LIMB_T(0xe4efd1ad3f767ceb),
			TO_LIMB_T(0x4a51d8667f0fe1cf), TO_LIMB_T(0x054fdf4bbf1d821c)  },
		  { TO_LIMB_T(0x72db2a50658d767b), TO_LIMB_T(0x8abf91faa257b3d5),
			TO_LIMB_T(0xe969d6833764ab47), TO_LIMB_T(0x464170142a1009eb),
			TO_LIMB_T(0xb14f01aadb30be2f), TO_LIMB_T(0x18ae6a856f40715d)  }
		};
		/* ...
		 * x_den = x'^10 + k_(2,9) * x'^9 + k_(2,8) * x'^8 + ... + k_(2,0)
		 */
		const Fp isogeny_map_x_den[] = { /*  (k_(2,*)<<384) % P  */
		  { TO_LIMB_T(0xb962a077fdb0f945), TO_LIMB_T(0xa6a9740fefda13a0),
			TO_LIMB_T(0xc14d568c3ed6c544), TO_LIMB_T(0xb43fc37b908b133e),
			TO_LIMB_T(0x9c0b3ac929599016), TO_LIMB_T(0x0165aa6c93ad115f)  },
		  { TO_LIMB_T(0x23279a3ba506c1d9), TO_LIMB_T(0x92cfca0a9465176a),
			TO_LIMB_T(0x3b294ab13755f0ff), TO_LIMB_T(0x116dda1c5070ae93),
			TO_LIMB_T(0xed4530924cec2045), TO_LIMB_T(0x083383d6ed81f1ce)  },
		  { TO_LIMB_T(0x9885c2a6449fecfc), TO_LIMB_T(0x4a2b54ccd37733f0),
			TO_LIMB_T(0x17da9ffd8738c142), TO_LIMB_T(0xa0fba72732b3fafd),
			TO_LIMB_T(0xff364f36e54b6812), TO_LIMB_T(0x0f29c13c660523e2)  },
		  { TO_LIMB_T(0xe349cc118278f041), TO_LIMB_T(0xd487228f2f3204fb),
			TO_LIMB_T(0xc9d325849ade5150), TO_LIMB_T(0x43a92bd69c15c2df),
			TO_LIMB_T(0x1c2c7844bc417be4), TO_LIMB_T(0x12025184f407440c)  },
		  { TO_LIMB_T(0x587f65ae6acb057b), TO_LIMB_T(0x1444ef325140201f),
			TO_LIMB_T(0xfbf995e71270da49), TO_LIMB_T(0xccda066072436a42),
			TO_LIMB_T(0x7408904f0f186bb2), TO_LIMB_T(0x13b93c63edf6c015)  },
		  { TO_LIMB_T(0xfb918622cd141920), TO_LIMB_T(0x4a4c64423ecaddb4),
			TO_LIMB_T(0x0beb232927f7fb26), TO_LIMB_T(0x30f94df6f83a3dc2),
			TO_LIMB_T(0xaeedd424d780f388), TO_LIMB_T(0x06cc402dd594bbeb)  },
		  { TO_LIMB_T(0xd41f761151b23f8f), TO_LIMB_T(0x32a92465435719b3),
			TO_LIMB_T(0x64f436e888c62cb9), TO_LIMB_T(0xdf70a9a1f757c6e4),
			TO_LIMB_T(0x6933a38d5b594c81), TO_LIMB_T(0x0c6f7f7237b46606)  },
		  { TO_LIMB_T(0x693c08747876c8f7), TO_LIMB_T(0x22c9850bf9cf80f0),
			TO_LIMB_T(0x8e9071dab950c124), TO_LIMB_T(0x89bc62d61c7baf23),
			TO_LIMB_T(0xbc6be2d8dad57c23), TO_LIMB_T(0x17916987aa14a122)  },
		  { TO_LIMB_T(0x1be3ff439c1316fd), TO_LIMB_T(0x9965243a7571dfa7),
			TO_LIMB_T(0xc7f7f62962f5cd81), TO_LIMB_T(0x32c6aa9af394361c),
			TO_LIMB_T(0xbbc2ee18e1c227f4), TO_LIMB_T(0x0c102cbac531bb34)  },
		  { TO_LIMB_T(0x997614c97bacbf07), TO_LIMB_T(0x61f86372b99192c0),
			TO_LIMB_T(0x5b8c95fc14353fc3), TO_LIMB_T(0xca2b066c2a87492f),
			TO_LIMB_T(0x16178f5bbf698711), TO_LIMB_T(0x12a6dcd7f0f4e0e8)  }
		};
		/*
		 * y = y' * y_num / y_den, where
		 * y_num = k_(3,15) * x'^15 + k_(3,14) * x'^14 + k_(3,13) * x'^13 +
		 *         ... + k_(3,0)
		 * ...
		 */
		const Fp isogeny_map_y_num[] = { /*  (k_(3,*)<<384) % P  */
		  { TO_LIMB_T(0x2b567ff3e2837267), TO_LIMB_T(0x1d4d9e57b958a767),
			TO_LIMB_T(0xce028fea04bd7373), TO_LIMB_T(0xcc31a30a0b6cd3df),
			TO_LIMB_T(0x7d7b18a682692693), TO_LIMB_T(0x0d300744d42a0310)  },
		  { TO_LIMB_T(0x99c2555fa542493f), TO_LIMB_T(0xfe7f53cc4874f878),
			TO_LIMB_T(0x5df0608b8f97608a), TO_LIMB_T(0x14e03832052b49c8),
			TO_LIMB_T(0x706326a6957dd5a4), TO_LIMB_T(0x0a8dadd9c2414555)  },
		  { TO_LIMB_T(0x13d942922a5cf63a), TO_LIMB_T(0x357e33e36e261e7d),
			TO_LIMB_T(0xcf05a27c8456088d), TO_LIMB_T(0x0000bd1de7ba50f0),
			TO_LIMB_T(0x83d0c7532f8c1fde), TO_LIMB_T(0x13f70bf38bbf2905)  },
		  { TO_LIMB_T(0x5c57fd95bfafbdbb), TO_LIMB_T(0x28a359a65e541707),
			TO_LIMB_T(0x3983ceb4f6360b6d), TO_LIMB_T(0xafe19ff6f97e6d53),
			TO_LIMB_T(0xb3468f4550192bf7), TO_LIMB_T(0x0bb6cde49d8ba257)  },
		  { TO_LIMB_T(0x590b62c7ff8a513f), TO_LIMB_T(0x314b4ce372cacefd),
			TO_LIMB_T(0x6bef32ce94b8a800), TO_LIMB_T(0x6ddf84a095713d5f),
			TO_LIMB_T(0x64eace4cb0982191), TO_LIMB_T(0x0386213c651b888d)  },
		  { TO_LIMB_T(0xa5310a31111bbcdd), TO_LIMB_T(0xa14ac0f5da148982),
			TO_LIMB_T(0xf9ad9cc95423d2e9), TO_LIMB_T(0xaa6ec095283ee4a7),
			TO_LIMB_T(0xcf5b1f022e1c9107), TO_LIMB_T(0x01fddf5aed881793)  },
		  { TO_LIMB_T(0x65a572b0d7a7d950), TO_LIMB_T(0xe25c2d8183473a19),
			TO_LIMB_T(0xc2fcebe7cb877dbd), TO_LIMB_T(0x05b2d36c769a89b0),
			TO_LIMB_T(0xba12961be86e9efb), TO_LIMB_T(0x07eb1b29c1dfde1f)  },
		  { TO_LIMB_T(0x93e09572f7c4cd24), TO_LIMB_T(0x364e929076795091),
			TO_LIMB_T(0x8569467e68af51b5), TO_LIMB_T(0xa47da89439f5340f),
			TO_LIMB_T(0xf4fa918082e44d64), TO_LIMB_T(0x0ad52ba3e6695a79)  },
		  { TO_LIMB_T(0x911429844e0d5f54), TO_LIMB_T(0xd03f51a3516bb233),
			TO_LIMB_T(0x3d587e5640536e66), TO_LIMB_T(0xfa86d2a3a9a73482),
			TO_LIMB_T(0xa90ed5adf1ed5537), TO_LIMB_T(0x149c9c326a5e7393)  },
		  { TO_LIMB_T(0x462bbeb03c12921a), TO_LIMB_T(0xdc9af5fa0a274a17),
			TO_LIMB_T(0x9a558ebde836ebed), TO_LIMB_T(0x649ef8f11a4fae46),
			TO_LIMB_T(0x8100e1652b3cdc62), TO_LIMB_T(0x1862bd62c291dacb)  },
		  { TO_LIMB_T(0x05c9b8ca89f12c26), TO_LIMB_T(0x0194160fa9b9ac4f),
			TO_LIMB_T(0x6a643d5a6879fa2c), TO_LIMB_T(0x14665bdd8846e19d),
			TO_LIMB_T(0xbb1d0d53af3ff6bf), TO_LIMB_T(0x12c7e1c3b28962e5)  },
		  { TO_LIMB_T(0xb55ebf900b8a3e17), TO_LIMB_T(0xfedc77ec1a9201c4),
			TO_LIMB_T(0x1f07db10ea1a4df4), TO_LIMB_T(0x0dfbd15dc41a594d),
			TO_LIMB_T(0x389547f2334a5391), TO_LIMB_T(0x02419f98165871a4)  },
		  { TO_LIMB_T(0xb416af000745fc20), TO_LIMB_T(0x8e563e9d1ea6d0f5),
			TO_LIMB_T(0x7c763e17763a0652), TO_LIMB_T(0x01458ef0159ebbef),
			TO_LIMB_T(0x8346fe421f96bb13), TO_LIMB_T(0x0d2d7b829ce324d2)  },
		  { TO_LIMB_T(0x93096bb538d64615), TO_LIMB_T(0x6f2a2619951d823a),
			TO_LIMB_T(0x8f66b3ea59514fa4), TO_LIMB_T(0xf563e63704f7092f),
			TO_LIMB_T(0x724b136c4cf2d9fa), TO_LIMB_T(0x046959cfcfd0bf49)  },
		  { TO_LIMB_T(0xea748d4b6e405346), TO_LIMB_T(0x91e9079c2c02d58f),
			TO_LIMB_T(0x41064965946d9b59), TO_LIMB_T(0xa06731f1d2bbe1ee),
			TO_LIMB_T(0x07f897e267a33f1b), TO_LIMB_T(0x1017290919210e5f)  },
		  { TO_LIMB_T(0x872aa6c17d985097), TO_LIMB_T(0xeecc53161264562a),
			TO_LIMB_T(0x07afe37afff55002), TO_LIMB_T(0x54759078e5be6838),
			TO_LIMB_T(0xc4b92d15db8acca8), TO_LIMB_T(0x106d87d1b51d13b9)  }
		};
		/* ...
		 * y_den = x'^15 + k_(4,14) * x'^14 + k_(4,13) * x'^13 + ... + k_(4,0)
		 */
		const Fp isogeny_map_y_den[] = { /*  (k_(4,*)<<384) % P  */
		  { TO_LIMB_T(0xeb6c359d47e52b1c), TO_LIMB_T(0x18ef5f8a10634d60),
			TO_LIMB_T(0xddfa71a0889d5b7e), TO_LIMB_T(0x723e71dcc5fc1323),
			TO_LIMB_T(0x52f45700b70d5c69), TO_LIMB_T(0x0a8b981ee47691f1)  },
		  { TO_LIMB_T(0x616a3c4f5535b9fb), TO_LIMB_T(0x6f5f037395dbd911),
			TO_LIMB_T(0xf25f4cc5e35c65da), TO_LIMB_T(0x3e50dffea3c62658),
			TO_LIMB_T(0x6a33dca523560776), TO_LIMB_T(0x0fadeff77b6bfe3e)  },
		  { TO_LIMB_T(0x2be9b66df470059c), TO_LIMB_T(0x24a2c159a3d36742),
			TO_LIMB_T(0x115dbe7ad10c2a37), TO_LIMB_T(0xb6634a652ee5884d),
			TO_LIMB_T(0x04fe8bb2b8d81af4), TO_LIMB_T(0x01c2a7a256fe9c41)  },
		  { TO_LIMB_T(0xf27bf8ef3b75a386), TO_LIMB_T(0x898b367476c9073f),
			TO_LIMB_T(0x24482e6b8c2f4e5f), TO_LIMB_T(0xc8e0bbd6fe110806),
			TO_LIMB_T(0x59b0c17f7631448a), TO_LIMB_T(0x11037cd58b3dbfbd)  },
		  { TO_LIMB_T(0x31c7912ea267eec6), TO_LIMB_T(0x1dbf6f1c5fcdb700),
			TO_LIMB_T(0xd30d4fe3ba86fdb1), TO_LIMB_T(0x3cae528fbee9a2a4),
			TO_LIMB_T(0xb1cce69b6aa9ad9a), TO_LIMB_T(0x044393bb632d94fb)  },
		  { TO_LIMB_T(0xc66ef6efeeb5c7e8), TO_LIMB_T(0x9824c289dd72bb55),
			TO_LIMB_T(0x71b1a4d2f119981d), TO_LIMB_T(0x104fc1aafb0919cc),
			TO_LIMB_T(0x0e49df01d942a628), TO_LIMB_T(0x096c3a09773272d4)  },
		  { TO_LIMB_T(0x9abc11eb5fadeff4), TO_LIMB_T(0x32dca50a885728f0),
			TO_LIMB_T(0xfb1fa3721569734c), TO_LIMB_T(0xc4b76271ea6506b3),
			TO_LIMB_T(0xd466a75599ce728e), TO_LIMB_T(0x0c81d4645f4cb6ed)  },
		  { TO_LIMB_T(0x4199f10e5b8be45b), TO_LIMB_T(0xda64e495b1e87930),
			TO_LIMB_T(0xcb353efe9b33e4ff), TO_LIMB_T(0x9e9efb24aa6424c6),
			TO_LIMB_T(0xf08d33680a237465), TO_LIMB_T(0x0d3378023e4c7406)  },
		  { TO_LIMB_T(0x7eb4ae92ec74d3a5), TO_LIMB_T(0xc341b4aa9fac3497),
			TO_LIMB_T(0x5be603899e907687), TO_LIMB_T(0x03bfd9cca75cbdeb),
			TO_LIMB_T(0x564c2935a96bfa93), TO_LIMB_T(0x0ef3c33371e2fdb5)  },
		  { TO_LIMB_T(0x7ee91fd449f6ac2e), TO_LIMB_T(0xe5d5bd5cb9357a30),
			TO_LIMB_T(0x773a8ca5196b1380), TO_LIMB_T(0xd0fda172174ed023),
			TO_LIMB_T(0x6cb95e0fa776aead), TO_LIMB_T(0x0d22d5a40cec7cff)  },
		  { TO_LIMB_T(0xf727e09285fd8519), TO_LIMB_T(0xdc9d55a83017897b),
			TO_LIMB_T(0x7549d8bd057894ae), TO_LIMB_T(0x178419613d90d8f8),
			TO_LIMB_T(0xfce95ebdeb5b490a), TO_LIMB_T(0x0467ffaef23fc49e)  },
		  { TO_LIMB_T(0xc1769e6a7c385f1b), TO_LIMB_T(0x79bc930deac01c03),
			TO_LIMB_T(0x5461c75a23ede3b5), TO_LIMB_T(0x6e20829e5c230c45),
			TO_LIMB_T(0x828e0f1e772a53cd), TO_LIMB_T(0x116aefa749127bff)  },
		  { TO_LIMB_T(0x101c10bf2744c10a), TO_LIMB_T(0xbbf18d053a6a3154),
			TO_LIMB_T(0xa0ecf39ef026f602), TO_LIMB_T(0xfc009d4996dc5153),
			TO_LIMB_T(0xb9000209d5bd08d3), TO_LIMB_T(0x189e5fe4470cd73c)  },
		  { TO_LIMB_T(0x7ebd546ca1575ed2), TO_LIMB_T(0xe47d5a981d081b55),
			TO_LIMB_T(0x57b2b625b6d4ca21), TO_LIMB_T(0xb0a1ba04228520cc),
			TO_LIMB_T(0x98738983c2107ff3), TO_LIMB_T(0x13dddbc4799d81d6)  },
		  { TO_LIMB_T(0x09319f2e39834935), TO_LIMB_T(0x039e952cbdb05c21),
			TO_LIMB_T(0x55ba77a9a2f76493), TO_LIMB_T(0xfd04e3dfc6086467),
			TO_LIMB_T(0xfb95832e7d78742e), TO_LIMB_T(0x0ef9c24eccaf5e0e)  }
		};
		Fp Zz_powers[15], map[15], xn, xd, yn, yd;

		/* lay down Z^2 powers in descending order                          */
		Fp::sqr(Zz_powers[14], p.z);                        /* ZZ^1         */
#ifdef __OPTIMIZE_SIZE__
		for (size_t i = 14; i > 0; i--)
			Fp::mul(Zz_powers[i - 1], Zz_powers[i], Zz_powers[14]);
#else
		Fp::sqr(Zz_powers[13], Zz_powers[14]);               /* ZZ^2  1+1    */
		Fp::mul(Zz_powers[12], Zz_powers[14], Zz_powers[13]);/* ZZ^3  2+1    */
		Fp::sqr(Zz_powers[11], Zz_powers[13]);               /* ZZ^4  2+2    */
		Fp::mul(Zz_powers[10], Zz_powers[13], Zz_powers[12]);/* ZZ^5  2+3    */
		Fp::sqr(Zz_powers[9], Zz_powers[12]);               /* ZZ^6  3+3    */
		Fp::mul(Zz_powers[8], Zz_powers[12], Zz_powers[11]);/* ZZ^7  3+4    */
		Fp::sqr(Zz_powers[7], Zz_powers[11]);               /* ZZ^8  4+4    */
		Fp::mul(Zz_powers[6], Zz_powers[11], Zz_powers[10]);/* ZZ^9  4+5    */
		Fp::sqr(Zz_powers[5], Zz_powers[10]);               /* ZZ^10 5+5    */
		Fp::mul(Zz_powers[4], Zz_powers[10], Zz_powers[9]); /* ZZ^11 5+6    */
		Fp::sqr(Zz_powers[3], Zz_powers[9]);                /* ZZ^12 6+6    */
		Fp::mul(Zz_powers[2], Zz_powers[9], Zz_powers[8]); /* ZZ^13 6+7    */
		Fp::sqr(Zz_powers[1], Zz_powers[8]);                /* ZZ^14 7+7    */
		Fp::mul(Zz_powers[0], Zz_powers[8], Zz_powers[7]); /* ZZ^15 7+8    */
#endif

		map_field_times_Zz<Fp>(map, isogeny_map_x_num, Zz_powers + 4, 11);
		Fp::mul(xn, p.x, isogeny_map_x_num[11]);
		Fp::add(xn, xn, map[10]);
		map_field<Fp>(xn, p.x, map, 10);

		map_field_times_Zz<Fp>(map, isogeny_map_x_den, Zz_powers + 5, 10);
		Fp::add(xd, p.x, map[9]);
		map_field<Fp>(xd, p.x, map, 9);
		Fp::mul(xd, xd, Zz_powers[14]);      /* xd *= Z^2                    */

		map_field_times_Zz<Fp>(map, isogeny_map_y_num, Zz_powers, 15);
		Fp::mul(yn, p.x, isogeny_map_y_num[15]);
		Fp::add(yn, yn, map[14]);
		map_field<Fp>(yn, p.x, map, 14);
		Fp::mul(yn, yn, p.y);               /* yn *= Y                      */

		map_field_times_Zz<Fp>(map, isogeny_map_y_den, Zz_powers, 15);
		Fp::add(yd, p.x, map[14]);
		map_field<Fp>(yd, p.x, map, 14);
		Fp::mul(Zz_powers[14], Zz_powers[14], p.z);
		Fp::mul(yd, yd, Zz_powers[14]);      /* yd *= Z^3                    */

		/* convert (xn, xd, yn, yd) to Jacobian coordinates                 */
		Fp::mul(out.z, xd, yd);             /* Z = xd * yd                  */
		Fp::mul(out.x, xn, yd);
		Fp::mul(out.x, out.x, out.z);     /* X = xn * xd * yd^2           */
		Fp::sqr(out.y, out.z);
		Fp::mul(out.y, out.y, xd);
		Fp::mul(out.y, out.y, yn);         /* Y = yn * xd^3 * yd^2         */
	}



	template<typename Group>
	BLST_HD BLST_INLINE void map_to_isogenous(Group& p, const typename Group::Field_type& u);

	template<>
	BLST_HD BLST_INLINE void map_to_isogenous<G1Point>(G1Point& p, const Fp& u)
	{
		const Fp minus_A = { /* P - A */
			TO_LIMB_T(0x8a9955f1650a005a), TO_LIMB_T(0x9865b3d192cfe93c),
			TO_LIMB_T(0xaed3ed0f3ef3c441), TO_LIMB_T(0x3c962ef33d92c442),
			TO_LIMB_T(0x22e438dbd74f94a2), TO_LIMB_T(0x04acbc265478c915)
		};
		const Fp Z = {       /* (11<<384) % P */
			TO_LIMB_T(0x886c00000023ffdc), TO_LIMB_T(0x0f70008d3090001d),
			TO_LIMB_T(0x77672417ed5828c3), TO_LIMB_T(0x9dac23e943dc1740),
			TO_LIMB_T(0x50553f1b9c131521), TO_LIMB_T(0x078c712fbe0ab6e8)
		};
		const Fp sqrt_minus_ZZZ = {
			TO_LIMB_T(0x43b571cad3215f1f), TO_LIMB_T(0xccb460ef1c702dc2),
			TO_LIMB_T(0x742d884f4f97100b), TO_LIMB_T(0xdb2c3e3238a3382b),
			TO_LIMB_T(0xe40f3fa13fce8f88), TO_LIMB_T(0x0073a2af9892a2ff)
		};
		const Fp ZxA = {
			TO_LIMB_T(0x7f674ea0a8915178), TO_LIMB_T(0xb0f945fc13b8fa65),
			TO_LIMB_T(0x4b46759a38e87d76), TO_LIMB_T(0x2e7a929641bbb6a1),
			TO_LIMB_T(0x1668ddfa462bf6b6), TO_LIMB_T(0x00960e2ed1cf294c)
		};
		Fp uu, tv2, x2n, gx1, gxd, y2;
# define xn     p.x
		auto& y = p.y;
# define xd     p.z
# define x1n    xn
# define y1     y
# define Zuu    x2n
# define tv4    y1
#define sgn0_fp(a) (Fp::sgn0_pty_mont(a) & 1)
		bool_t e1, e2;

		/*
		 * as per map_to_curve() from poc/sswu_opt.sage at
		 * https://github.com/cfrg/draft-irtf-cfrg-hash-to-curve
		 */
		 /* x numerator variants                                             */
		Fp::sqr(uu, u);                      /* uu = u^2                     */
		Fp::mul(Zuu, Z, uu);                 /* Zuu = Z * uu                 */
		Fp::sqr(tv2, Zuu);                   /* tv2 = Zuu^2                  */
		Fp::add(tv2, tv2, Zuu);              /* tv2 = tv2 + Zuu              */
		Fp::add(x1n, tv2, get_BLST_Rx().p);   /* x1n = tv2 + 1                */
		Fp::mul(x1n, x1n, get_BLST_Bprime_E1());        /* x1n = x1n * B                */
		Fp::mul(x2n, Zuu, x1n);              /* x2n = Zuu * x1n              */

		/* x denumenator                                                    */
		Fp::mul(xd, minus_A, tv2);           /* xd = -A * tv2                */
		e1 = vec_is_zero(xd.l, sizeof(xd));   /* e1 = xd == 0                 */
		vec_select(xd.l, ZxA.l, xd.l, sizeof(xd), e1);    /*              # If xd == 0, set xd = Z*A */

		/* y numerators variants                                            */
		Fp::sqr(tv2, xd);                    /* tv2 = xd^2                   */
		Fp::mul(gxd, xd, tv2);               /* gxd = xd^3                   */
		Fp::mul(tv2, get_BLST_Aprime_E1(), tv2);        /* tv2 = A * tv2                */
		Fp::sqr(gx1, x1n);                   /* gx1 = x1n^2                  */
		Fp::add(gx1, gx1, tv2);              /* gx1 = gx1 + tv2      # x1n^2 + A*xd^2 */
		Fp::mul(gx1, gx1, x1n);              /* gx1 = gx1 * x1n      # x1n^3 + A*x1n*xd^2 */
		Fp::mul(tv2, get_BLST_Bprime_E1(), gxd);        /* tv2 = B * gxd                */
		Fp::add(gx1, gx1, tv2);              /* gx1 = gx1 + tv2      # x1^3 + A*x1*xd^2 + B*xd^3 */
		Fp::sqr(tv4, gxd);                   /* tv4 = gxd^2                  */
		Fp::mul(tv2, gx1, gxd);              /* tv2 = gx1 * gxd              */
		Fp::mul(tv4, tv4, tv2);              /* tv4 = tv4 * tv2      # gx1*gxd^3 */
		e2 = Fp::recip_sqrt(y1, tv4);        /* y1 = tv4^c1          # (gx1*gxd^3)^((p-3)/4) */
		Fp::mul(y1, y1, tv2);                /* y1 = y1 * tv2        # gx1*gxd*y1 */
		Fp::mul(y2, y1, sqrt_minus_ZZZ);     /* y2 = y1 * c2         # y2 = y1*sqrt(-Z^3) */
		Fp::mul(y2, y2, uu);                 /* y2 = y2 * uu                 */
		Fp::mul(y2, y2, u);                  /* y2 = y2 * u                  */

		/* choose numerators                                                */
		vec_select(xn.l, x1n.l, x2n.l, sizeof(xn), e2);   /* xn = e2 ? x1n : x2n  */
		vec_select(y.l, y1.l, y2.l, sizeof(y), e2);       /* y  = e2 ? y1 : y2    */

		e1 = sgn0_fp(u);
		e2 = sgn0_fp(y);
		Fp::cneg(y, y, e1 ^ e2);               /* fix sign of y                */
		/* return (xn, xd, y, 1)        */

	/* convert (xn, xd, y, 1) to Jacobian projective coordinates        */
		Fp::mul(p.x, xn, xd);               /* X = xn * xd                  */
		Fp::mul(p.y, y, gxd);               /* Y = y * xd^3                 */
#ifndef xd
		vec_copy(p.z, xd, sizeof(xd));     /* Z = xd                       */
#else
# undef xn
# undef xd
# undef x1n
# undef y1
# undef Zuu
# undef tv4
#endif
#undef sgn0_fp
	}




	template<>
	BLST_HD BLST_INLINE void map_to_isogenous<G2Point>(G2Point& p, const Fp2& u)
	{
		const Fp2 minus_A = {
		  Fp{ 0 },
		  Fp{ TO_LIMB_T(0xd4c4fffffcec5869), TO_LIMB_T(0x1da3f3eed25bfd79),
			TO_LIMB_T(0x7fa833c5136fff67), TO_LIMB_T(0x59261433cd540cbd),
			TO_LIMB_T(0x48450f5f2b84682c), TO_LIMB_T(0x07e05d00bf959233) }
		};
		const Fp2 Z = {              /* -2 - i */
		  Fp{ TO_LIMB_T(0x87ebfffffff9555c), TO_LIMB_T(0x656fffe5da8ffffa),
			TO_LIMB_T(0x0fd0749345d33ad2), TO_LIMB_T(0xd951e663066576f4),
			TO_LIMB_T(0xde291a3d41e980d3), TO_LIMB_T(0x0815664c7dfe040d) },
		  Fp{ TO_LIMB_T(0x43f5fffffffcaaae), TO_LIMB_T(0x32b7fff2ed47fffd),
			TO_LIMB_T(0x07e83a49a2e99d69), TO_LIMB_T(0xeca8f3318332bb7a),
			TO_LIMB_T(0xef148d1ea0f4c069), TO_LIMB_T(0x040ab3263eff0206) }
		};
		const Fp2 ZxA = {            /* 240 - 480*i */
		  Fp{ TO_LIMB_T(0xe53a000003135242), TO_LIMB_T(0x01080c0fdef80285),
			TO_LIMB_T(0xe7889edbe340f6bd), TO_LIMB_T(0x0b51375126310601),
			TO_LIMB_T(0x02d6985717c744ab), TO_LIMB_T(0x1220b4e979ea5467) },
		  Fp{ TO_LIMB_T(0xa989fffff9d8b0d2), TO_LIMB_T(0x3b47e7dda4b7faf3),
			TO_LIMB_T(0xff50678a26dffece), TO_LIMB_T(0xb24c28679aa8197a),
			TO_LIMB_T(0x908a1ebe5708d058), TO_LIMB_T(0x0fc0ba017f2b2466) }
		};
		Fp2 uu, tv2, tv4, x2n, gx1, gxd, y2;
# define xn     p.x
		auto& y = p.y;
# define xd     p.z
# define x1n    xn
# define y1     y
# define Zuu    x2n
#define sgn0_fp2(a) (Fp2::sgn0_pty_mont(a) & 1)
		bool_t e1, e2;

		/*
		 * as per map_to_curve() from poc/sswu_opt.sage at
		 * https://github.com/cfrg/draft-irtf-cfrg-hash-to-curve
		 * with 9mod16 twists...
		 */
		 /* x numerator variants                                             */
		Fp2::sqr(uu, u);                           /* uu = u^2                     */
		Fp2::mul(Zuu, Z, uu);                      /* Zuu = Z * uu                 */
		Fp2::sqr(tv2, Zuu);                        /* tv2 = Zuu^2                  */
		Fp2::add(tv2, tv2, Zuu);                   /* tv2 = tv2 + Zuu              */
		Fp2::add(x1n, tv2, *(Fp2*)get_BLST_Rx().p2);     /* x1n = tv2 + 1                */
		Fp2::mul(x1n, x1n, get_BLST_Bprime_E2()); /* x1n = x1n * B                */
		Fp2::mul(x2n, Zuu, x1n);                   /* x2n = Zuu * x1n              */

		/* x denumenator                                                    */
		Fp2::mul(xd, minus_A, tv2);          /* xd = -A * tv2                */
		e1 = vec_is_zero(&(xd), sizeof(xd));   /* e1 = xd == 0                 */
		vec_select(&(xd), &ZxA, &(xd), sizeof(xd), e1);    /*              # If xd == 0, set xd = Z*A */

		/* y numerators variants                                            */
		Fp2::sqr(tv2, xd);                   /* tv2 = xd^2                   */
		Fp2::mul(gxd, xd, tv2);              /* gxd = xd^3                   */
		Fp2::mul(tv2, get_BLST_Aprime_E2(), tv2);       /* tv2 = A * tv2                */
		Fp2::sqr(gx1, x1n);                  /* gx1 = x1n^2                  */
		Fp2::add(gx1, gx1, tv2);             /* gx1 = gx1 + tv2      # x1n^2 + A*xd^2 */
		Fp2::mul(gx1, gx1, x1n);             /* gx1 = gx1 * x1n      # x1n^3 + A*x1n*xd^2 */
		Fp2::mul(tv2, get_BLST_Bprime_E2(), gxd);       /* tv2 = B * gxd                */
		Fp2::add(gx1, gx1, tv2);             /* gx1 = gx1 + tv2      # x1^3 + A*x1*xd^2 + B*xd^3 */
		Fp2::sqr(tv4, gxd);                  /* tv4 = gxd^2                  */
		Fp2::mul(tv2, gx1, gxd);             /* tv2 = gx1 * gxd              */
		Fp2::mul(tv4, tv4, tv2);             /* tv4 = tv4 * tv2      # gx1*gxd^3 */
		e2 = Fp2::recip_sqrt(y1, tv4);        /* y1 = tv4^c1          # (gx1*gxd^3)^((p^2-9)/16) */
		Fp2::mul(y1, y1, tv2);               /* y1 = y1 * tv2        # gx1*gxd*y1 */
		Fp2::mul(y2, y1, uu);                /* y2 = y1 * uu                 */
		Fp2::mul(y2, y2, u);                 /* y2 = y2 * u                  */

		/* choose numerators                                                */
		vec_select(&(xn), &(x1n), &x2n, sizeof(xn), e2);   /* xn = e2 ? x1n : x2n  */
		vec_select(&y, &y1, &y2, sizeof(y), e2);       /* y  = e2 ? y1 : y2    */

		e1 = sgn0_fp2(u);
		e2 = sgn0_fp2(y);
		Fp2::cneg(y, y, e1 ^ e2);              /* fix sign of y                */
		/* return (xn, xd, y, 1)        */

	/* convert (xn, xd, y, 1) to Jacobian projective coordinates        */
		Fp2::mul(p.x, xn, xd);              /* X = xn * xd                  */
		Fp2::mul(p.y, y, gxd);              /* Y = y * xd^3                 */
#ifndef xd
		vec_copy(p.z, xd, sizeof(xd));     /* Z = xd                       */
#else
# undef xn
# undef xd
# undef x1n
# undef y1
# undef Zuu
# undef tv4
#endif
#undef sgn0_fp2
	}


	template<typename Group>
	BLST_HD BLST_INLINE void clear_cofactor(Group& out, const Group& p)
	{
		if constexpr (std::is_same_v<Group, G1Point>) {
			/* clear the cofactor by multiplying |p| by 1-z, 0xd201000000010001 */
			times_minus_z<Group>(out, p);
			Group::dadd(out, out, p, NULL);
		}
		else
		{

			G2Point t0, t1;

			/* A.Budroni, F.Pintore, "Efficient hash maps to G2 on BLS curves"  */
			G2Point::doubling(out, p);           /* out = 2P                     */
			psi(out, out);                      /* out = Ψ(2P)                  */
			psi(out, out);                      /* out = Ψ²(2P)                 */

			vec_copy(&t0, &p, sizeof(t0));
			G2Point::cneg(t0, 1);             /* t0 = -P                      */
			psi(t1, t0);                      /* t1 = -Ψ(P)                   */
			G2Point::dadd(out, out, t0, NULL);/* out = Ψ²(2P) - P             */
			G2Point::dadd(out, out, t1, NULL);/* out = Ψ²(2P) - P - Ψ(P)      */

			times_minus_z<G2Point>(t0, p);    /* t0 = [-z]P                   */
			G2Point::dadd(t0, t0, p, NULL);  /* t0 = [-z + 1]P               */
			G2Point::dadd(t0, t0, t1, NULL);/* t0 = [-z + 1]P - Ψ(P)        */
			times_minus_z<G2Point>(t1, t0);  /* t1 = [z² - z]P + [z]Ψ(P)     */
			G2Point::dadd(out, out, t1, NULL);/* out = [z² - z - 1]P          */
			/*     + [z - 1]Ψ(P)            */
			/*     + Ψ²(2P)                 */

		}
	}

	template<typename Group>
	BLST_HD BLST_INLINE auto& get_Aprime()
	{
		if constexpr (std::is_same_v<Group, G1Point>) {
			return get_BLST_Aprime_E1();
		}
		else {
			return get_BLST_Aprime_E2();
		}
	}


	/*
	*  map_to_curve maps one or two field elements to a point on the
	*  curve in the specified group.  If two elements are provided,
	*  and these are the output of a hash function, then the
	*  group element is expected to be indistinguishable from a
	*  random element of the group.
	*
	 * |u|, |v| are expected to be in Montgomery representation
	 */
	template<typename Group>
	BLST_HD BLST_INLINE void map_to_curve(
		Group& out,
		const typename Group::Field_type& u,
		const typename Group::Field_type* v)
	{
		Group p;

		map_to_isogenous<Group>(p, u);

		if (v != NULL) {
			map_to_isogenous<Group>(out, *v);    /* borrow |out|                 */
			Group::dadd(p, p, out, &get_Aprime<Group>());
		}

		isogeny_map_to<Group>(p, p);          /* sprinkle isogenous powder    */

		clear_cofactor<Group>(out, p);
	}


	template<typename Group>
	BLST_HD BLST_INLINE void Encode_to_curve(Group& p, const unsigned char* msg, size_t msg_len,
		const unsigned char* DST, size_t DST_len,
		const unsigned char* aug, size_t aug_len)
	{
		using G = std::remove_cvref_t<Group>;
		typename G::Field_type u[1];
		constexpr int nelem = 1 + 1 * std::is_same_v<G, G2Point>;

		hash_to_field<nelem>((Fp*)u, aug, aug_len, msg, msg_len, DST, DST_len);
		map_to_curve<G>(p, u[0], NULL);
	}

	template<typename Group>
	BLST_HD BLST_INLINE void Hash_to_curve(Group& p, const unsigned char* msg, size_t msg_len,
		const unsigned char* DST, size_t DST_len,
		const unsigned char* aug, size_t aug_len)
	{
		typename Group::Field_type u[2];

		constexpr int nelem = 2 + 2 * std::is_same_v<Group, G2Point>;
		hash_to_field<nelem>((Fp*)u, aug, aug_len, msg, msg_len, DST, DST_len);
		map_to_curve<Group>(p, u[0], &u[1]);
	}

}
