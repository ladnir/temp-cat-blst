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
            use rand_core::{RngCore, SeedableRng};
            use rand_chacha::ChaCha20Rng;

            #[test]
            fn test_mult() {
                const npoints: usize = 2000;
                const nbits: usize = 160;
                const nbytes: usize = (nbits + 7) / 8;

                let mut scalars = Box::new([0u8; nbytes * npoints]);
                ChaCha20Rng::from_seed([0u8; 32]).fill_bytes(scalars.as_mut());

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
                        $add_or_double(&mut naive, &naive, &t);
                    }
                    if i < 27 {
                        let points = $points::from(&points[0..i + 1]);
                        assert_eq!(naive, points.mult(scalars.as_ref(), nbits));
                    }
                }

                let points = $points::from(&points);

                assert_eq!(naive, points.mult(scalars.as_ref(), nbits));
            }

            #[cfg(all(feature = "std", not(feature = "no-threads")))]
            #[test]
            fn test_oversubscribed_mult() {
                // Before 6271918, this combination entered breakdown() with a
                // three-bit window and divided by a zero-width subdivision.
                const npoints: usize = 32;
                const nbits: usize = 64;
                const nbytes: usize = (nbits + 7) / 8;
                const reported_ncpus: usize = 1024;

                let mut scalars = [0u8; nbytes * npoints];
                ChaCha20Rng::from_seed([1u8; 32]).fill_bytes(&mut scalars);

                let mut points = vec![<$point>::default(); npoints];
                let mut naive = <$point>::default();
                for i in 0..npoints {
                    unsafe {
                        let mut term = <$point>::default();
                        $mult(
                            &mut points[i],
                            $generator(),
                            scalars.as_ptr().add(i * nbytes),
                            32,
                        );
                        $mult(
                            &mut term,
                            &points[i],
                            scalars.as_ptr().add(i * nbytes),
                            nbits,
                        );
                        $add_or_double(&mut naive, &naive, &term);
                    }
                }

                let points = $points::from(&points);
                let actual = with_pippenger_ncpus(reported_ncpus, || {
                    points.mult(&scalars, nbits)
                });
                assert_eq!(naive, actual);
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
            fn test_add() {
                const npoints: usize = 2000;
                const nbits: usize = 32;
                const nbytes: usize = (nbits + 7) / 8;

                let mut scalars = Box::new([0u8; nbytes * npoints]);
                ChaCha20Rng::from_seed([0u8; 32]).fill_bytes(scalars.as_mut());

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
                        $add_or_double(&mut naive, &naive, &points[i]);
                    }
                }

                let points = $points::from(&points);
                assert_eq!(naive, points.add());
            }
        }
    };
}
