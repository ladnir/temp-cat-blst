// Copyright Supranational LLC
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

use core::ops::{Index, IndexMut};
use core::slice::SliceIndex;

#[inline]
fn scalar_layout(npoints: usize, nbits: usize) -> Option<(usize, usize)> {
    let nbytes = nbits / 8 + usize::from(nbits % 8 != 0);
    nbytes
        .checked_mul(npoints)
        .map(|scalar_bytes| (nbytes, scalar_bytes))
}

/// # Safety
///
/// Every element through `values.len()` must have been initialized as `T`.
#[inline]
unsafe fn assume_init_vec<T>(mut values: Vec<MaybeUninit<T>>) -> Vec<T> {
    let ptr = values.as_mut_ptr().cast::<T>();
    let len = values.len();
    let capacity = values.capacity();
    core::mem::forget(values);
    Vec::from_raw_parts(ptr, len, capacity)
}

macro_rules! pippenger_mult_impl {
    (
        $points:ident,
        $point:ty,
        $point_affine:ty,
        $to_affines:ident,
        $scratch_sizeof:ident,
        $multi_scalar_mult:ident,
        $add_or_double:ident,
        $double:ident,
        $test_mod:ident,
        $generator:ident,
        $mult:ident,
        $add:ident,
        $is_inf:ident,
        $in_group:ident,
    ) => {
        pub struct $points {
            points: Vec<$point_affine>,
        }

        impl<I: SliceIndex<[$point_affine]>> Index<I> for $points {
            type Output = I::Output;

            #[inline]
            fn index(&self, i: I) -> &Self::Output {
                &self.points[i]
            }
        }
        impl<I: SliceIndex<[$point_affine]>> IndexMut<I> for $points {
            #[inline]
            fn index_mut(&mut self, i: I) -> &mut Self::Output {
                &mut self.points[i]
            }
        }

        impl $points {
            #[inline]
            pub fn as_slice(&self) -> &[$point_affine] {
                self.points.as_slice()
            }

            pub fn from(points: &[$point]) -> Self {
                let npoints = points.len();
                if npoints == 0 {
                    return Self { points: Vec::new() };
                }

                let mut affine = Vec::with_capacity(npoints);
                affine.resize_with(npoints, MaybeUninit::uninit);
                let p: [*const $point; 2] = [points.as_ptr(), ptr::null()];
                unsafe {
                    $to_affines(
                        affine.as_mut_ptr().cast::<$point_affine>(),
                        &p[0],
                        npoints,
                    );
                    Self {
                        points: assume_init_vec(affine),
                    }
                }
            }

            /// Computes a variable-time MSM of these prevalidated subgroup
            /// points and public little-endian scalars. Each scalar occupies
            /// `ceil(nbits / 8)` consecutive bytes; `nbits` must not exceed
            /// 256. Empty points or zero bits return the identity. Panics when
            /// the scalar buffer is too short or its required size overflows.
            #[inline]
            pub fn mult(&self, scalars: &[u8], nbits: usize) -> $point {
                self.as_slice().mult(scalars, nbits)
            }

            #[inline]
            /// Adds these public subgroup points with a variable-time batch
            /// algorithm. Allocation failure inside the C boundary falls back
            /// to complete serial addition and does not change the result.
            pub fn add(&self) -> $point {
                self.as_slice().add()
            }
        }

        impl MultiPoint for [$point_affine] {
            type Output = $point;

            /// Computes a variable-time MSM over public scalar data.
            ///
            /// Every point must already be on-curve and in the prime-order
            /// subgroup. Scalars are packed contiguously in little-endian
            /// order, with `ceil(nbits / 8)` bytes per point. Bits above
            /// `nbits` are ignored and `nbits` must not exceed 256. Empty
            /// input or zero bits returns the identity. This method panics if
            /// the scalar buffer is too short or its required size overflows.
            fn mult(&self, scalars: &[u8], nbits: usize) -> $point {
                let npoints = self.len();
                if nbits > 256 {
                    panic!("scalar bit length exceeds 256 bits");
                }
                let (_, scalar_bytes) = scalar_layout(npoints, nbits)
                    .unwrap_or_else(|| panic!("scalars length mismatch"));

                if scalars.len() < scalar_bytes {
                    panic!("scalars length mismatch");
                }
                if npoints == 0 || nbits == 0 {
                    return <$point>::default();
                }

                let p: [*const $point_affine; 2] = [self.as_ptr(), ptr::null()];
                let s: [*const u8; 2] = [scalars.as_ptr(), ptr::null()];

                let mut ret = <$point>::default();
                unsafe {
                    let scratch_words = $scratch_sizeof(npoints) / 8;
                    let mut scratch = Vec::with_capacity(scratch_words);
                    scratch
                        .resize_with(scratch_words, MaybeUninit::<u64>::uninit);
                    $multi_scalar_mult(
                        &mut ret,
                        &p[0],
                        npoints,
                        &s[0],
                        nbits,
                        scratch.as_mut_ptr().cast::<u64>(),
                    );
                }
                ret
            }

            /// Adds these public subgroup points with variable-time work that
            /// may depend on point values. Empty input returns the identity.
            fn add(&self) -> $point {
                let npoints = self.len();

                if npoints == 0 {
                    return <$point>::default();
                }

                let p: [*const _; 2] = [self.as_ptr(), ptr::null()];
                let mut ret = <$point>::default();
                unsafe { $add(&mut ret, &p[0], npoints) };

                ret
            }

            fn validate(&self) -> Result<(), BLST_ERROR> {
                for i in 0..self.len() {
                    if unsafe { $is_inf(&self[i]) } {
                        return Err(BLST_ERROR::BLST_PK_IS_INFINITY);
                    }
                    if !unsafe { $in_group(&self[i]) } {
                        return Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP);
                    }
                }
                Ok(())
            }
        }

        #[cfg(test)]
        pippenger_test_mod!(
            $test_mod,
            $points,
            $point,
            $add_or_double,
            $generator,
            $mult,
        );
    };
}

#[cfg(test)]
include!("pippenger-test_mod.rs");

pippenger_mult_impl!(
    p1_affines,
    blst_p1,
    blst_p1_affine,
    blst_p1s_to_affine,
    blst_p1s_mult_pippenger_scratch_sizeof,
    blst_p1s_mult_pippenger,
    blst_p1_add_or_double,
    blst_p1_double,
    p1_multi_point,
    blst_p1_generator,
    blst_p1_mult,
    blst_p1s_add,
    blst_p1_affine_is_inf,
    blst_p1_affine_in_g1,
);

pippenger_mult_impl!(
    p2_affines,
    blst_p2,
    blst_p2_affine,
    blst_p2s_to_affine,
    blst_p2s_mult_pippenger_scratch_sizeof,
    blst_p2s_mult_pippenger,
    blst_p2_add_or_double,
    blst_p2_double,
    p2_multi_point,
    blst_p2_generator,
    blst_p2_mult,
    blst_p2s_add,
    blst_p2_affine_is_inf,
    blst_p2_affine_in_g2,
);

#[cfg(test)]
mod pippenger_layout_tests {
    use super::*;

    #[test]
    fn scalar_layout_rejects_overflow() {
        assert_eq!(scalar_layout(0, usize::MAX), Some((usize::MAX / 8 + 1, 0)));
        assert_eq!(
            scalar_layout(1, usize::MAX),
            Some((usize::MAX / 8 + 1, usize::MAX / 8 + 1))
        );
        assert_eq!(scalar_layout(usize::MAX, usize::MAX), None);
    }
}
