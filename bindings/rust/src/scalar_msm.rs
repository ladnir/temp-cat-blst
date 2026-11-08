//! Allocation-free evaluation of constant-time scalar-backend MSMs.

use alloc::vec::Vec;
use core::{
    mem::{size_of, MaybeUninit},
    ptr,
    sync::atomic::{compiler_fence, Ordering},
};

use super::{G1Affine, G1Projective, G2Affine, G2Projective, Scalar};
use crate::{
    blst_fr, blst_p1, blst_p1_affine, blst_p1s_mult_ct_affine_precompute_fr,
    blst_p1s_mult_ct_affine_precomputed_fr,
    blst_p1s_mult_ct_affine_precomputed_fr_table_sizeof, blst_p1s_mult_ct_fr,
    blst_p1s_mult_ct_fr_scratch_sizeof, blst_p2, blst_p2_affine,
    blst_p2s_mult_ct_fr, blst_p2s_mult_ct_fr_scratch_sizeof,
};

const _: () = {
    assert!(size_of::<u64>() == 8);
    assert!(size_of::<Scalar>() == size_of::<blst_fr>());
    assert!(size_of::<G1Affine>() == size_of::<blst_p1_affine>());
    assert!(size_of::<G2Affine>() == size_of::<blst_p2_affine>());
};

/// Errors from caller-owned constant-time MSM storage and evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CtMsmError {
    InvalidCapacity,
    AllocationFailed,
    LengthMismatch { points: usize, scalars: usize },
    CapacityExceeded { capacity: usize, actual: usize },
    BackendFailure,
}

struct AlignedWords {
    capacity: usize,
    bytes: usize,
    words: Vec<MaybeUninit<u64>>,
}

impl AlignedWords {
    fn new(capacity: usize, bytes: usize) -> Result<Self, CtMsmError> {
        if (capacity != 0 && bytes == 0)
            || bytes % size_of::<u64>() != 0
        {
            return Err(CtMsmError::InvalidCapacity);
        }
        let mut words = Vec::new();
        words
            .try_reserve_exact(bytes / size_of::<u64>())
            .map_err(|_| CtMsmError::AllocationFailed)?;
        // Native MSM initializes every scratch word before reading it.
        unsafe { words.set_len(bytes / size_of::<u64>()) };
        Ok(Self {
            capacity,
            bytes,
            words,
        })
    }

    #[inline]
    fn as_mut_ptr(&mut self) -> *mut u64 {
        self.words.as_mut_ptr().cast::<u64>()
    }

    fn check(&self, points: usize, scalars: usize) -> Result<(), CtMsmError> {
        if points != scalars {
            return Err(CtMsmError::LengthMismatch { points, scalars });
        }
        if points > self.capacity {
            return Err(CtMsmError::CapacityExceeded {
                capacity: self.capacity,
                actual: points,
            });
        }
        Ok(())
    }

    fn clear(&mut self) {
        for word in &mut self.words {
            // Volatile writes do not read potentially uninitialized scratch.
            unsafe { ptr::write_volatile(word.as_mut_ptr(), 0) };
        }
        compiler_fence(Ordering::SeqCst);
    }
}

/// Reusable, limb-aligned G1 constant-time MSM scratch.
pub struct G1CtMsmWorkspace(AlignedWords);

impl G1CtMsmWorkspace {
    pub fn required_bytes(max_points: usize) -> Result<usize, CtMsmError> {
        let bytes = unsafe { blst_p1s_mult_ct_fr_scratch_sizeof(max_points) };
        if (max_points != 0 && bytes == 0)
            || bytes % size_of::<u64>() != 0
        {
            Err(CtMsmError::InvalidCapacity)
        } else {
            Ok(bytes)
        }
    }

    pub fn new(max_points: usize) -> Result<Self, CtMsmError> {
        let bytes = Self::required_bytes(max_points)?;
        Ok(Self(AlignedWords::new(max_points, bytes)?))
    }

    pub fn capacity(&self) -> usize {
        self.0.capacity
    }

    /// Clears the owned scratch; also runs when this workspace is dropped.
    pub fn clear_secret(&mut self) {
        self.0.clear();
    }

    /// Computes an MSM with no allocation or scalar-dependent control flow.
    pub fn msm(
        &mut self,
        points: &[G1Affine],
        scalars: &[Scalar],
    ) -> Result<G1Projective, CtMsmError> {
        self.0.check(points.len(), scalars.len())?;
        if points.is_empty() {
            return Ok(G1Projective::identity());
        }
        let mut out = MaybeUninit::<blst_p1>::uninit();
        let scratch = self.0.as_mut_ptr();
        // Transparent wrappers are native-compatible; inputs, output, and
        // owned scratch are disjoint and capacity was checked above.
        let ok = unsafe {
            blst_p1s_mult_ct_fr(
                out.as_mut_ptr(),
                points.as_ptr().cast::<blst_p1_affine>(),
                scalars.as_ptr().cast::<blst_fr>(),
                points.len(),
                scratch,
                self.0.bytes,
            )
        };
        if !ok {
            return Err(CtMsmError::BackendFailure);
        }
        Ok(unsafe { G1Projective::from_raw_unchecked(out.assume_init()) })
    }
}

impl Drop for G1CtMsmWorkspace {
    fn drop(&mut self) {
        self.clear_secret();
    }
}

/// Reusable, limb-aligned G2 constant-time MSM scratch.
pub struct G2CtMsmWorkspace(AlignedWords);

impl G2CtMsmWorkspace {
    pub fn required_bytes(max_points: usize) -> Result<usize, CtMsmError> {
        let bytes = unsafe { blst_p2s_mult_ct_fr_scratch_sizeof(max_points) };
        if (max_points != 0 && bytes == 0)
            || bytes % size_of::<u64>() != 0
        {
            Err(CtMsmError::InvalidCapacity)
        } else {
            Ok(bytes)
        }
    }

    pub fn new(max_points: usize) -> Result<Self, CtMsmError> {
        let bytes = Self::required_bytes(max_points)?;
        Ok(Self(AlignedWords::new(max_points, bytes)?))
    }

    pub fn capacity(&self) -> usize {
        self.0.capacity
    }

    pub fn clear_secret(&mut self) {
        self.0.clear();
    }

    pub fn msm(
        &mut self,
        points: &[G2Affine],
        scalars: &[Scalar],
    ) -> Result<G2Projective, CtMsmError> {
        self.0.check(points.len(), scalars.len())?;
        if points.is_empty() {
            return Ok(G2Projective::identity());
        }
        let mut out = MaybeUninit::<blst_p2>::uninit();
        let scratch = self.0.as_mut_ptr();
        let ok = unsafe {
            blst_p2s_mult_ct_fr(
                out.as_mut_ptr(),
                points.as_ptr().cast::<blst_p2_affine>(),
                scalars.as_ptr().cast::<blst_fr>(),
                points.len(),
                scratch,
                self.0.bytes,
            )
        };
        if !ok {
            return Err(CtMsmError::BackendFailure);
        }
        Ok(unsafe { G2Projective::from_raw_unchecked(out.assume_init()) })
    }
}

impl Drop for G2CtMsmWorkspace {
    fn drop(&mut self) {
        self.clear_secret();
    }
}

/// Immutable precomputation for repeated G1 MSMs on public affine bases.
/// Any prefix of the original bases may be evaluated without rebuilding it.
pub struct G1FixedCtMsmTable(AlignedWords);

impl G1FixedCtMsmTable {
    pub fn new(points: &[G1Affine]) -> Result<Self, CtMsmError> {
        let bytes = unsafe {
            blst_p1s_mult_ct_affine_precomputed_fr_table_sizeof(points.len())
        };
        let mut table = AlignedWords::new(points.len(), bytes)?;
        if !points.is_empty() {
            let ok = unsafe {
                blst_p1s_mult_ct_affine_precompute_fr(
                    table.as_mut_ptr(),
                    bytes,
                    points.as_ptr().cast::<blst_p1_affine>(),
                    points.len(),
                )
            };
            if !ok {
                return Err(CtMsmError::BackendFailure);
            }
        }
        Ok(Self(table))
    }

    pub fn capacity(&self) -> usize {
        self.0.capacity
    }

    pub fn msm(
        &self,
        scalars: &[Scalar],
        scratch: &mut G1CtMsmWorkspace,
    ) -> Result<G1Projective, CtMsmError> {
        if scalars.len() > self.0.capacity {
            return Err(CtMsmError::CapacityExceeded {
                capacity: self.0.capacity,
                actual: scalars.len(),
            });
        }
        scratch.0.check(scalars.len(), scalars.len())?;
        if scalars.is_empty() {
            return Ok(G1Projective::identity());
        }
        let mut out = MaybeUninit::<blst_p1>::uninit();
        let scratch_ptr = scratch.0.as_mut_ptr();
        let ok = unsafe {
            blst_p1s_mult_ct_affine_precomputed_fr(
                out.as_mut_ptr(),
                self.0.words.as_ptr().cast::<u64>(),
                self.0.bytes,
                self.0.capacity,
                scalars.as_ptr().cast::<blst_fr>(),
                scalars.len(),
                scratch_ptr,
                scratch.0.bytes,
            )
        };
        if !ok {
            return Err(CtMsmError::BackendFailure);
        }
        Ok(unsafe { G1Projective::from_raw_unchecked(out.assume_init()) })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_ct_msm_matches_individual_products_and_checks_capacity() {
        let points = [
            G1Affine::identity(),
            G1Affine::generator(),
            (G1Projective::generator() * Scalar::from(3)).to_affine(),
        ];
        let scalars = [Scalar::from(17), Scalar::ZERO, Scalar::from(11)];
        let mut scratch = G1CtMsmWorkspace::new(3).unwrap();
        assert_eq!(scratch.msm(&[], &[]), Ok(G1Projective::identity()));
        let expected = G1Projective::from(points[2]) * scalars[2];
        assert_eq!(scratch.msm(&points, &scalars), Ok(expected));
        assert_eq!(
            scratch.msm(&points, &scalars[..2]),
            Err(CtMsmError::LengthMismatch {
                points: 3,
                scalars: 2,
            })
        );
        let mut small = G1CtMsmWorkspace::new(2).unwrap();
        assert_eq!(
            small.msm(&points, &scalars),
            Err(CtMsmError::CapacityExceeded {
                capacity: 2,
                actual: 3,
            })
        );
        scratch.clear_secret();
        assert!(scratch
            .0
            .words
            .iter()
            .all(|word| unsafe { *word.assume_init_ref() == 0 }));
    }

    #[test]
    fn fixed_g1_table_matches_direct_msm_for_every_prefix() {
        let generator = G1Projective::generator();
        let points = [
            generator.to_affine(),
            (generator * Scalar::from(2)).to_affine(),
            G1Affine::identity(),
            (generator * Scalar::from(7)).to_affine(),
        ];
        let scalars = [
            Scalar::from(13),
            Scalar::ZERO,
            Scalar::from(5),
            Scalar::from(19),
        ];
        let table = G1FixedCtMsmTable::new(&points).unwrap();
        let mut scratch = G1CtMsmWorkspace::new(points.len()).unwrap();
        assert_eq!(table.capacity(), 4);
        for len in 0..=points.len() {
            let expected =
                scratch.msm(&points[..len], &scalars[..len]).unwrap();
            assert_eq!(table.msm(&scalars[..len], &mut scratch), Ok(expected));
        }
        assert_eq!(G1FixedCtMsmTable::new(&[]).unwrap().capacity(), 0);
    }

    #[test]
    fn g2_ct_msm_matches_individual_products() {
        let points = [G2Affine::identity(), G2Affine::generator()];
        let scalars = [Scalar::from(9), Scalar::from(23)];
        let mut scratch = G2CtMsmWorkspace::new(2).unwrap();
        assert_eq!(
            scratch.msm(&points, &scalars),
            Ok(G2Projective::generator() * Scalar::from(23))
        );
        scratch.clear_secret();
        assert!(scratch
            .0
            .words
            .iter()
            .all(|word| unsafe { *word.assume_init_ref() == 0 }));
    }
}
