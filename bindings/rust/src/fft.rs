// Copyright Supranational LLC
// Copyright Category Labs
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

//! Wrappers for the owned BLST Fr and group FFTs.

use crate::*;

extern "C" {
    fn blst_fr_fft(values: *mut blst_fr, n: usize, inverse: bool) -> bool;
    fn blst_fr_fft_domain_len(n: usize) -> usize;
    fn blst_fr_fft_domain_init(
        domain: *mut blst_fr,
        domain_len: usize,
        n: usize,
    ) -> bool;
    fn blst_fr_fft_with_domain(
        values: *mut blst_fr,
        n: usize,
        domain: *const blst_fr,
        domain_len: usize,
        inverse: bool,
    ) -> bool;
    fn blst_p1_fft(values: *mut blst_p1, n: usize, inverse: bool) -> bool;
    fn blst_p2_fft(values: *mut blst_p2, n: usize, inverse: bool) -> bool;
    fn blst_p1_fft_with_domain(
        values: *mut blst_p1,
        n: usize,
        domain: *const blst_fr,
        domain_len: usize,
        inverse: bool,
    ) -> bool;
    fn blst_p2_fft_with_domain(
        values: *mut blst_p2,
        n: usize,
        domain: *const blst_fr,
        domain_len: usize,
        inverse: bool,
    ) -> bool;
    fn blst_p1_ifft_unscaled_with_domain(
        values: *mut blst_p1,
        n: usize,
        domain: *const blst_fr,
        domain_len: usize,
    ) -> bool;
    fn blst_p2_ifft_unscaled_with_domain(
        values: *mut blst_p2,
        n: usize,
        domain: *const blst_fr,
        domain_len: usize,
    ) -> bool;
    fn blst_fr_interpolate(
        out_coefficients: *mut blst_fr,
        points: *const blst_fr,
        values: *const blst_fr,
        n: usize,
    ) -> bool;
    fn blst_fr_interpolate_fft_domain_suffix_scratch_len(
        n: usize,
        fft_len: usize,
    ) -> usize;
    fn blst_fr_interpolate_fft_domain_suffix_with_scratch(
        out_suffix: *mut blst_fr,
        indices: *const usize,
        values: *const blst_fr,
        n: usize,
        ell: usize,
        fft_len: usize,
        scratch: *mut blst_fr,
        scratch_len: usize,
    ) -> bool;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FftError {
    EmptyInput,
    InvalidSize,
    BackendFailure,
}

/// Runs the BLST reference FFT in-place over Montgomery-form Fr values.
///
/// Every input must be reduced and Montgomery encoded. `values.len()` must be a
/// non-zero power of two no larger than `2^30`. Forward FFT uses powers of the
/// primitive root for the size; inverse FFT uses reciprocal powers and applies
/// the usual `1/n` scaling inside BLST. Root preparation may allocate. A
/// validation or allocation failure leaves values unchanged; an unexpected
/// native failure after the in-place transform begins may leave them partially
/// transformed. Work and addresses depend on public length/direction, not
/// coefficient values.
pub fn fr_fft(values: &mut [blst_fr], inverse: bool) -> Result<(), FftError> {
    check_fft_input(values)?;

    let ok = unsafe { blst_fr_fft(values.as_mut_ptr(), values.len(), inverse) };
    if ok {
        Ok(())
    } else {
        Err(FftError::BackendFailure)
    }
}

/// Returns the number of `blst_fr` elements required by a prepared domain.
/// Valid sizes are non-zero powers of two through `2^30`; the query allocates
/// nothing.
pub fn fr_fft_domain_len(n: usize) -> Result<usize, FftError> {
    check_fft_size(n)?;

    let len = unsafe { blst_fr_fft_domain_len(n) };
    if len == 0 {
        Err(FftError::BackendFailure)
    } else {
        Ok(len)
    }
}

/// Initializes exact-size caller-owned storage for a prepared FFT domain.
///
/// Domain construction fills forward root powers and the inverse-size scale in
/// Montgomery form. The returned storage may be reused for any number of
/// transforms of exactly `n` values. It allocates nothing. The domain must be
/// treated as unspecified if initialization returns an error and may be shared
/// immutably after success.
pub fn fr_fft_domain_init(
    domain: &mut [blst_fr],
    n: usize,
) -> Result<(), FftError> {
    let required = fr_fft_domain_len(n)?;
    if domain.len() != required {
        return Err(FftError::InvalidSize);
    }

    let ok = unsafe {
        blst_fr_fft_domain_init(domain.as_mut_ptr(), domain.len(), n)
    };
    if ok {
        Ok(())
    } else {
        Err(FftError::BackendFailure)
    }
}

/// Runs an allocation-free in-place FFT using a prepared caller-owned domain.
///
/// `domain` must have been initialized by `fr_fft_domain_init` for exactly the
/// current value count. Forward and inverse transforms both reuse its one
/// forward-root table; inverse execution includes the prepared `1/n` scale.
/// Values must be reduced Montgomery Fr elements and must not overlap domain.
/// The operation allocates nothing. Invalid shape/tag failures leave values
/// unchanged; an unexpected failure after in-place work begins may leave them
/// partially transformed. An immutable domain may be shared by concurrent calls
/// operating on disjoint mutable value slices.
pub fn fr_fft_with_domain(
    values: &mut [blst_fr],
    domain: &[blst_fr],
    inverse: bool,
) -> Result<(), FftError> {
    check_fft_input(values)?;
    if domain.len() != fr_fft_domain_len(values.len())? {
        return Err(FftError::InvalidSize);
    }

    let ok = unsafe {
        blst_fr_fft_with_domain(
            values.as_mut_ptr(),
            values.len(),
            domain.as_ptr(),
            domain.len(),
            inverse,
        )
    };
    if ok {
        Ok(())
    } else {
        Err(FftError::BackendFailure)
    }
}

/// Runs the in-place radix-two FFT on projective G1 points. Root preparation
/// may allocate. Group multiplication is variable-time.
///
/// # Safety
/// Every input must be a valid G1 subgroup point (infinity is allowed). Points
/// must be public; this operation is not constant-time in their values.
pub unsafe fn p1_fft(
    values: &mut [blst_p1],
    inverse: bool,
) -> Result<(), FftError> {
    check_fft_size(values.len())?;
    if unsafe { blst_p1_fft(values.as_mut_ptr(), values.len(), inverse) } {
        Ok(())
    } else {
        Err(FftError::BackendFailure)
    }
}

/// Runs the in-place radix-two FFT on projective G2 points. Root preparation
/// may allocate. Group multiplication is variable-time.
///
/// # Safety
/// Every input must be a valid G2 subgroup point (infinity is allowed). Points
/// must be public; this operation is not constant-time in their values.
pub unsafe fn p2_fft(
    values: &mut [blst_p2],
    inverse: bool,
) -> Result<(), FftError> {
    check_fft_size(values.len())?;
    if unsafe { blst_p2_fft(values.as_mut_ptr(), values.len(), inverse) } {
        Ok(())
    } else {
        Err(FftError::BackendFailure)
    }
}

/// Runs an allocation-free in-place G1 FFT with a prepared Fr domain.
///
/// # Safety
/// Every input must be a valid public G1 subgroup point. `domain` must have
/// been initialized for this size by `fr_fft_domain_init` and kept unchanged.
pub unsafe fn p1_fft_with_domain(
    values: &mut [blst_p1],
    domain: &[blst_fr],
    inverse: bool,
) -> Result<(), FftError> {
    check_fft_size(values.len())?;
    if domain.len() != fr_fft_domain_len(values.len())? {
        return Err(FftError::InvalidSize);
    }
    if unsafe {
        blst_p1_fft_with_domain(
            values.as_mut_ptr(),
            values.len(),
            domain.as_ptr(),
            domain.len(),
            inverse,
        )
    } {
        Ok(())
    } else {
        Err(FftError::BackendFailure)
    }
}

/// Runs an allocation-free in-place G2 FFT with a prepared Fr domain.
///
/// # Safety
/// Every input must be a valid public G2 subgroup point. `domain` must have
/// been initialized for this size by `fr_fft_domain_init` and kept unchanged.
pub unsafe fn p2_fft_with_domain(
    values: &mut [blst_p2],
    domain: &[blst_fr],
    inverse: bool,
) -> Result<(), FftError> {
    check_fft_size(values.len())?;
    if domain.len() != fr_fft_domain_len(values.len())? {
        return Err(FftError::InvalidSize);
    }
    if unsafe {
        blst_p2_fft_with_domain(
            values.as_mut_ptr(),
            values.len(),
            domain.as_ptr(),
            domain.len(),
            inverse,
        )
    } {
        Ok(())
    } else {
        Err(FftError::BackendFailure)
    }
}

/// Runs an allocation-free G1 inverse FFT without its `1/n` scale.
///
/// # Safety
/// Every input must be a valid public G1 subgroup point. `domain` must have
/// been initialized for this size by `fr_fft_domain_init` and kept unchanged.
pub unsafe fn p1_ifft_unscaled_with_domain(
    values: &mut [blst_p1],
    domain: &[blst_fr],
) -> Result<(), FftError> {
    check_fft_size(values.len())?;
    if domain.len() != fr_fft_domain_len(values.len())? {
        return Err(FftError::InvalidSize);
    }
    if unsafe {
        blst_p1_ifft_unscaled_with_domain(
            values.as_mut_ptr(),
            values.len(),
            domain.as_ptr(),
            domain.len(),
        )
    } {
        Ok(())
    } else {
        Err(FftError::BackendFailure)
    }
}

/// Runs an allocation-free G2 inverse FFT without its `1/n` scale.
///
/// # Safety
/// Every input must be a valid public G2 subgroup point. `domain` must have
/// been initialized for this size by `fr_fft_domain_init` and kept unchanged.
pub unsafe fn p2_ifft_unscaled_with_domain(
    values: &mut [blst_p2],
    domain: &[blst_fr],
) -> Result<(), FftError> {
    check_fft_size(values.len())?;
    if domain.len() != fr_fft_domain_len(values.len())? {
        return Err(FftError::InvalidSize);
    }
    if unsafe {
        blst_p2_ifft_unscaled_with_domain(
            values.as_mut_ptr(),
            values.len(),
            domain.as_ptr(),
            domain.len(),
        )
    } {
        Ok(())
    } else {
        Err(FftError::BackendFailure)
    }
}

/// Interpolates monomial coefficients from distinct Fr points and values.
///
/// `output.len()` must equal `points.len() == values.len()`. The output buffer
/// is caller-owned; no coefficient vector is returned across the FFI boundary.
/// Points and values must be reduced Montgomery-form Fr elements.
/// It is unchanged on error. The operation may allocate and is variable-time in
/// the points and public length; for fixed points, reconstruction shape does not
/// depend on values. Ordinary heap storage contains value-derived intermediates
/// and is not securely erased.
pub fn fr_interpolate(
    points: &[blst_fr],
    values: &[blst_fr],
    output: &mut [blst_fr],
) -> Result<(), FftError> {
    if points.is_empty() {
        return Err(FftError::EmptyInput);
    }
    if points.len() != values.len() || points.len() != output.len() {
        return Err(FftError::InvalidSize);
    }

    let ok = unsafe {
        blst_fr_interpolate(
            output.as_mut_ptr(),
            points.as_ptr(),
            values.as_ptr(),
            points.len(),
        )
    };
    if ok {
        Ok(())
    } else {
        Err(FftError::BackendFailure)
    }
}

/// Returns the caller-owned scratch length, in `blst_fr` elements, for
/// `fr_interpolate_fft_domain_suffix_with_scratch`.
pub fn fr_interpolate_fft_domain_suffix_scratch_len(
    n: usize,
    fft_len: usize,
) -> Result<usize, FftError> {
    if n == 0 {
        return Err(FftError::EmptyInput);
    }
    if n > fft_len || !fft_len.is_power_of_two() || fft_len > i32::MAX as usize
    {
        return Err(FftError::InvalidSize);
    }

    let len = unsafe {
        blst_fr_interpolate_fft_domain_suffix_scratch_len(n, fft_len)
    };
    if len == 0 {
        Err(FftError::InvalidSize)
    } else {
        Ok(len)
    }
}

/// Scratch-backed suffix interpolation over FFT-domain samples.
///
/// The operation allocates nothing. Its schedule depends on the public lengths
/// and indices, not on values. Output is unchanged on error, while scratch may
/// be modified after native shape validation succeeds. Scratch contains
/// value-derived intermediates and must be cleared by the caller if required.
/// Values must be reduced Montgomery-form Fr elements. Independent calls are
/// concurrency-safe because the wrapper exposes only caller-owned mutable
/// output and scratch.
pub fn fr_interpolate_fft_domain_suffix_with_scratch(
    indices: &[usize],
    values: &[blst_fr],
    ell: usize,
    fft_len: usize,
    output: &mut [blst_fr],
    scratch: &mut [blst_fr],
) -> Result<(), FftError> {
    if indices.is_empty() {
        return Err(FftError::EmptyInput);
    }
    if indices.len() != values.len()
        || ell > indices.len()
        || output.len() != indices.len() - ell
        || indices.len() > fft_len
        || !fft_len.is_power_of_two()
        || fft_len > i32::MAX as usize
    {
        return Err(FftError::InvalidSize);
    }

    let required =
        fr_interpolate_fft_domain_suffix_scratch_len(indices.len(), fft_len)?;
    if scratch.len() < required {
        return Err(FftError::InvalidSize);
    }

    let ok = unsafe {
        blst_fr_interpolate_fft_domain_suffix_with_scratch(
            output.as_mut_ptr(),
            indices.as_ptr(),
            values.as_ptr(),
            indices.len(),
            ell,
            fft_len,
            scratch.as_mut_ptr(),
            scratch.len(),
        )
    };
    if ok {
        Ok(())
    } else {
        Err(FftError::BackendFailure)
    }
}

fn check_fft_size(n: usize) -> Result<(), FftError> {
    if n == 0 {
        return Err(FftError::EmptyInput);
    }
    if !n.is_power_of_two() || n > i32::MAX as usize {
        return Err(FftError::InvalidSize);
    }

    Ok(())
}

fn check_fft_input(values: &[blst_fr]) -> Result<(), FftError> {
    check_fft_size(values.len())
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;
    use alloc::vec::Vec;

    fn fr_from_u64(value: u64) -> blst_fr {
        let words = [value, 0, 0, 0];
        let mut out = blst_fr::default();
        unsafe {
            blst_fr_from_uint64(&mut out, words.as_ptr());
        }
        out
    }

    fn fr_add(lhs: &blst_fr, rhs: &blst_fr) -> blst_fr {
        let mut out = blst_fr::default();
        unsafe {
            blst_fr_add(&mut out, lhs, rhs);
        }
        out
    }

    fn fr_sub(lhs: &blst_fr, rhs: &blst_fr) -> blst_fr {
        let mut out = blst_fr::default();
        unsafe {
            blst_fr_sub(&mut out, lhs, rhs);
        }
        out
    }

    fn fr_mul(lhs: &blst_fr, rhs: &blst_fr) -> blst_fr {
        let mut out = blst_fr::default();
        unsafe {
            blst_fr_mul(&mut out, lhs, rhs);
        }
        out
    }

    #[test]
    fn rejects_invalid_sizes() {
        assert_eq!(fr_fft(&mut [], false), Err(FftError::EmptyInput));

        let mut values = [fr_from_u64(1); 3];
        assert_eq!(fr_fft(&mut values, false), Err(FftError::InvalidSize));
    }

    #[test]
    fn round_trips_forward_then_inverse() {
        let input = [
            fr_from_u64(3),
            fr_from_u64(1),
            fr_from_u64(4),
            fr_from_u64(1),
            fr_from_u64(5),
            fr_from_u64(9),
            fr_from_u64(2),
            fr_from_u64(6),
        ];
        let mut values = input;

        fr_fft(&mut values, false).unwrap();
        assert_ne!(values, input);
        fr_fft(&mut values, true).unwrap();

        assert_eq!(values, input);
    }

    #[test]
    fn prepared_domain_rejects_invalid_shapes() {
        assert_eq!(fr_fft_domain_len(0), Err(FftError::EmptyInput));
        assert_eq!(fr_fft_domain_len(3), Err(FftError::InvalidSize));

        let mut short_domain = [fr_from_u64(0); 8];
        assert_eq!(
            fr_fft_domain_init(&mut short_domain, 8),
            Err(FftError::InvalidSize)
        );

        let mut values = [fr_from_u64(1); 8];
        assert_eq!(
            fr_fft_with_domain(&mut values, &short_domain, false),
            Err(FftError::InvalidSize)
        );
    }

    #[test]
    fn prepared_fft_matches_allocating_fft_and_round_trips() {
        for n in [1usize, 2, 4, 8, 32] {
            let input: Vec<_> = (0..n)
                .map(|i| fr_from_u64((i as u64 + 1) * 17 + 3))
                .collect();
            let mut domain =
                vec![fr_from_u64(0); fr_fft_domain_len(n).unwrap()];
            fr_fft_domain_init(&mut domain, n).unwrap();

            let mut expected = input.clone();
            fr_fft(&mut expected, false).unwrap();

            let mut actual = input.clone();
            fr_fft_with_domain(&mut actual, &domain, false).unwrap();
            assert_eq!(actual, expected, "forward mismatch for n={n}");

            fr_fft_with_domain(&mut actual, &domain, true).unwrap();
            assert_eq!(actual, input, "round-trip mismatch for n={n}");
        }
    }

    #[test]
    fn prepared_fft_matches_natural_order_dft() {
        const N: usize = 8;
        let input = [
            fr_from_u64(2),
            fr_from_u64(7),
            fr_from_u64(1),
            fr_from_u64(8),
            fr_from_u64(2),
            fr_from_u64(8),
            fr_from_u64(1),
            fr_from_u64(8),
        ];
        let mut domain = [fr_from_u64(0); N + 2];
        fr_fft_domain_init(&mut domain, N).unwrap();

        let mut expected = [fr_from_u64(0); N];
        for (evaluation_index, evaluation) in expected.iter_mut().enumerate() {
            for (coefficient_index, coefficient) in input.iter().enumerate() {
                let root = &domain[(evaluation_index * coefficient_index) % N];
                *evaluation = fr_add(evaluation, &fr_mul(coefficient, root));
            }
        }

        let mut actual = input;
        fr_fft_with_domain(&mut actual, &domain, false).unwrap();
        assert_eq!(actual, expected);
    }

    #[test]
    fn group_ffts_round_trip_with_prepared_domain() {
        const N: usize = 4;
        let mut domain = [blst_fr::default(); N + 2];
        fr_fft_domain_init(&mut domain, N).unwrap();

        let mut g1 = [unsafe { *blst_p1_generator() }; N];
        let mut g2 = [unsafe { *blst_p2_generator() }; N];
        for i in 1..N {
            let previous_g1 = g1[i - 1];
            let previous_g2 = g2[i - 1];
            unsafe {
                blst_p1_add_or_double(
                    &mut g1[i],
                    &previous_g1,
                    blst_p1_generator(),
                );
                blst_p2_add_or_double(
                    &mut g2[i],
                    &previous_g2,
                    blst_p2_generator(),
                );
            }
        }
        let input_g1 = g1;
        let input_g2 = g2;
        unsafe {
            p1_fft_with_domain(&mut g1, &domain, false).unwrap();
            p2_fft_with_domain(&mut g2, &domain, false).unwrap();
            let mut unscaled_g1 = g1;
            let mut unscaled_g2 = g2;
            p1_ifft_unscaled_with_domain(&mut unscaled_g1, &domain).unwrap();
            p2_ifft_unscaled_with_domain(&mut unscaled_g2, &domain).unwrap();
            for i in 0..N {
                let mut doubled_g1 = blst_p1::default();
                let mut quadrupled_g1 = blst_p1::default();
                blst_p1_add_or_double(
                    &mut doubled_g1,
                    &input_g1[i],
                    &input_g1[i],
                );
                blst_p1_add_or_double(
                    &mut quadrupled_g1,
                    &doubled_g1,
                    &doubled_g1,
                );
                assert!(blst_p1_is_equal(&unscaled_g1[i], &quadrupled_g1));
                let mut doubled_g2 = blst_p2::default();
                let mut quadrupled_g2 = blst_p2::default();
                blst_p2_add_or_double(
                    &mut doubled_g2,
                    &input_g2[i],
                    &input_g2[i],
                );
                blst_p2_add_or_double(
                    &mut quadrupled_g2,
                    &doubled_g2,
                    &doubled_g2,
                );
                assert!(blst_p2_is_equal(&unscaled_g2[i], &quadrupled_g2));
            }
            p1_fft_with_domain(&mut g1, &domain, true).unwrap();
            p2_fft_with_domain(&mut g2, &domain, true).unwrap();
            for i in 0..N {
                assert!(blst_p1_is_equal(&g1[i], &input_g1[i]));
                assert!(blst_p2_is_equal(&g2[i], &input_g2[i]));
            }
        }
    }

    #[test]
    fn size_two_matches_dft() {
        let a = fr_from_u64(11);
        let b = fr_from_u64(7);
        let mut values = [a, b];

        fr_fft(&mut values, false).unwrap();

        assert_eq!(values, [fr_add(&a, &b), fr_sub(&a, &b)]);
    }

    #[test]
    fn interpolate_recovers_monomial_coefficients() {
        let coefficients = [
            fr_from_u64(9),
            fr_from_u64(2),
            fr_from_u64(7),
            fr_from_u64(1),
            fr_from_u64(8),
        ];
        let points = [
            fr_from_u64(1),
            fr_from_u64(2),
            fr_from_u64(3),
            fr_from_u64(4),
            fr_from_u64(5),
        ];
        let mut values = [fr_from_u64(0); 5];
        for (out, point) in values.iter_mut().zip(&points) {
            let mut acc = fr_from_u64(0);
            for coefficient in coefficients.iter().rev() {
                acc = fr_add(&fr_mul(&acc, point), coefficient);
            }
            *out = acc;
        }

        let mut recovered = [fr_from_u64(0); 5];
        fr_interpolate(&points, &values, &mut recovered).unwrap();

        assert_eq!(recovered, coefficients);
    }

    #[test]
    fn suffix_interpolation_recovers_requested_coefficients() {
        let coefficients = [
            fr_from_u64(13),
            fr_from_u64(21),
            fr_from_u64(34),
            fr_from_u64(55),
            fr_from_u64(89),
        ];
        let fft_len = 8;
        let mut domain = [fr_from_u64(0); 8];
        domain[1] = fr_from_u64(1);
        fr_fft(&mut domain, false).unwrap();

        let indices = [0usize, 1, 3, 4, 6];
        let mut values = [fr_from_u64(0); 5];
        for (out, index) in values.iter_mut().zip(indices) {
            let point = &domain[index];
            let mut acc = fr_from_u64(0);
            for coefficient in coefficients.iter().rev() {
                acc = fr_add(&fr_mul(&acc, point), coefficient);
            }
            *out = acc;
        }

        let scratch_len = fr_interpolate_fft_domain_suffix_scratch_len(
            indices.len(),
            fft_len,
        )
        .unwrap();
        let mut scratch = vec![fr_from_u64(0); scratch_len];
        let mut scratch_suffix = [fr_from_u64(0); 3];
        fr_interpolate_fft_domain_suffix_with_scratch(
            &indices,
            &values,
            2,
            fft_len,
            &mut scratch_suffix,
            &mut scratch,
        )
        .unwrap();
        assert_eq!(scratch_suffix, coefficients[2..]);
    }

    #[test]
    fn suffix_interpolation_rejects_invalid_shapes() {
        let oversized_indices = [0usize; 9];
        let oversized_values = [fr_from_u64(1); 9];
        let mut oversized_output = [fr_from_u64(0); 9];
        assert_eq!(
            fr_interpolate_fft_domain_suffix_scratch_len(9, 8),
            Err(FftError::InvalidSize)
        );
        let mut oversized_scratch = [fr_from_u64(0); 1];
        assert_eq!(
            fr_interpolate_fft_domain_suffix_with_scratch(
                &oversized_indices,
                &oversized_values,
                0,
                8,
                &mut oversized_output,
                &mut oversized_scratch,
            ),
            Err(FftError::InvalidSize)
        );

        let indices = [0usize, 1];
        let values = [fr_from_u64(1), fr_from_u64(2)];
        let mut output = [fr_from_u64(0); 2];
        let required =
            fr_interpolate_fft_domain_suffix_scratch_len(2, 2).unwrap();
        let mut short_scratch = vec![fr_from_u64(0); required - 1];
        assert_eq!(
            fr_interpolate_fft_domain_suffix_with_scratch(
                &indices,
                &values,
                0,
                2,
                &mut output,
                &mut short_scratch,
            ),
            Err(FftError::InvalidSize)
        );
    }
}
