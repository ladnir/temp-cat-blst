// Copyright Supranational LLC
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

macro_rules! pippenger_test_mod {
    (
        $test_mod:ident,
        $points:ident,
        $point:ty,
        $add_or_double:ident,
        $generator:ident,
        $mult:ident,
    ) => {
        mod $test_mod {
            use super::*;
            use rand_chacha::ChaCha20Rng;
            use rand_core::{RngCore, SeedableRng};

            #[test]
            fn test_mult() {
                const npoints: usize = 2000;
                const nbits: usize = 160;
                const nbytes: usize = (nbits + 7) / 8;

                let mut scalars = vec![0u8; nbytes * npoints];
                ChaCha20Rng::from_seed([0u8; 32]).fill_bytes(&mut scalars);

                let mut points = vec![<$point>::default(); npoints];

                let mut naive = <$point>::default();
                for i in 0..npoints {
                    unsafe {
                        let mut t = <$point>::default();
                        $mult(
                            &mut points[i],
                            $generator(),
                            scalars.as_ptr().add(i * nbytes),
                            core::cmp::min(32, nbits),
                        );
                        $mult(
                            &mut t,
                            &points[i],
                            scalars.as_ptr().add(i * nbytes),
                            nbits,
                        );
                        // The C boundary explicitly supports exact output/input
                        // aliasing. Raw pointers avoid manufacturing overlapping
                        // Rust references while exercising that contract.
                        $add_or_double(
                            core::ptr::addr_of_mut!(naive),
                            core::ptr::addr_of!(naive),
                            core::ptr::addr_of!(t),
                        );
                    }
                    if i < 27 {
                        let points = $points::from(&points[0..i + 1]);
                        assert_eq!(naive, points.mult(&scalars, nbits));
                    }
                }

                let points = $points::from(&points);

                assert_eq!(naive, points.mult(&scalars, nbits));
            }

            #[test]
            fn test_empty_operations() {
                let points = $points::from(&[]);

                assert_eq!(<$point>::default(), points.mult(&[], 255));
                assert_eq!(<$point>::default(), points.mult(&[], 0));
                assert_eq!(<$point>::default(), points.add());
            }

            #[test]
            fn test_zero_bit_mult() {
                let points = $points::from(&[unsafe { *$generator() }]);
                assert_eq!(<$point>::default(), points.mult(&[], 0));
            }

            #[test]
            fn test_scalar_width_boundaries() {
                let generator = unsafe { *$generator() };
                let points = $points::from(&[generator]);

                let low_nine_bits = [0x5au8, 0x01];
                let high_bits_set = [0x5au8, 0xff];
                assert_eq!(
                    points.mult(&low_nine_bits, 9),
                    points.mult(&high_bits_set, 9)
                );

                let mut scalar = [0u8; 32];
                ChaCha20Rng::from_seed([1u8; 32]).fill_bytes(&mut scalar);
                let expected = unsafe {
                    let mut out = <$point>::default();
                    $mult(
                        core::ptr::addr_of_mut!(out),
                        core::ptr::addr_of!(generator),
                        scalar.as_ptr(),
                        256,
                    );
                    out
                };
                assert_eq!(expected, points.mult(&scalar, 256));
            }

            #[test]
            #[should_panic(expected = "scalars length mismatch")]
            fn test_short_scalar_buffer_panics() {
                let generator = unsafe { *$generator() };
                let points = $points::from(&[generator, generator]);
                let _ = points.mult(&[0u8; 3], 9);
            }

            #[test]
            #[should_panic(expected = "scalar bit length exceeds 256 bits")]
            fn test_oversized_scalar_width_panics() {
                let points = $points::from(&[unsafe { *$generator() }]);
                let _ = points.mult(&[0u8; 33], 257);
            }

            #[test]
            fn test_add() {
                const npoints: usize = 2000;
                const nbits: usize = 32;
                const nbytes: usize = (nbits + 7) / 8;

                let mut scalars = vec![0u8; nbytes * npoints];
                ChaCha20Rng::from_seed([0u8; 32]).fill_bytes(&mut scalars);

                let mut points = vec![<$point>::default(); npoints];

                let mut naive = <$point>::default();
                for i in 0..npoints {
                    unsafe {
                        $mult(
                            &mut points[i],
                            $generator(),
                            scalars.as_ptr().add(i * nbytes),
                            32,
                        );
                        $add_or_double(
                            core::ptr::addr_of_mut!(naive),
                            core::ptr::addr_of!(naive),
                            core::ptr::addr_of!(points[i]),
                        );
                    }
                }

                let points = $points::from(&points);
                assert_eq!(naive, points.add());
            }

            #[test]
            fn test_validate_reports_infinity_deterministically() {
                let valid = $points::from(&[unsafe { *$generator() }]);
                assert_eq!(valid.as_slice().validate(), Ok(()));

                let infinity = $points::from(&[<$point>::default()]);
                assert_eq!(
                    infinity.as_slice().validate(),
                    Err(BLST_ERROR::BLST_PK_IS_INFINITY)
                );
            }
        }
    };
}
