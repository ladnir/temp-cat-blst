// Copyright Supranational LLC
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

//! Small safe Rust wrapper for the owned BLST Fr FFT.

use crate::*;

extern "C" {
    fn blst_fr_fft(values: *mut blst_fr, n: usize, inverse: bool) -> bool;
    fn blst_fr_fft_transpose_forward(values: *mut blst_fr, n: usize) -> bool;
    fn blst_fr_interpolate(
        out_coefficients: *mut blst_fr,
        points: *const blst_fr,
        values: *const blst_fr,
        n: usize,
    ) -> bool;
    fn blst_fr_interpolate_fft_domain(
        out_coefficients: *mut blst_fr,
        indices: *const usize,
        values: *const blst_fr,
        n: usize,
        fft_len: usize,
    ) -> bool;
    fn blst_fr_interpolate_fft_domain_suffix(
        out_suffix: *mut blst_fr,
        indices: *const usize,
        values: *const blst_fr,
        n: usize,
        ell: usize,
        fft_len: usize,
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
/// `values.len()` must be a non-zero power of two. Forward FFT uses powers of
/// the primitive root for the size; inverse FFT uses reciprocal powers and
/// applies the usual `1/n` scaling inside BLST.
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
/// transforms of exactly `n` values.
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

/// Runs the transpose-equivalent of the forward BLST Fr FFT in-place.
///
/// This is the verifier kernel for BVE-style Vandermonde contractions. It uses
/// BLST's native full-size FFT implementation rather than materializing the
/// Vandermonde matrix in Rust.
pub fn fr_fft_transpose_forward(
    values: &mut [blst_fr],
) -> Result<(), FftError> {
    check_fft_input(values)?;

    let ok = unsafe {
        blst_fr_fft_transpose_forward(values.as_mut_ptr(), values.len())
    };
    if ok {
        Ok(())
    } else {
        Err(FftError::BackendFailure)
    }
}

/// Interpolates monomial coefficients from distinct Fr points and values.
///
/// `output.len()` must equal `points.len() == values.len()`. The output buffer
/// is caller-owned; no coefficient vector is returned across the FFI boundary.
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

/// Interpolates monomial coefficients from values at FFT-domain indices.
///
/// `indices[i]` selects the same root-of-unity point used by `fr_fft` at that
/// index in a domain of length `fft_len`.
pub fn fr_interpolate_fft_domain(
    indices: &[usize],
    values: &[blst_fr],
    fft_len: usize,
    output: &mut [blst_fr],
) -> Result<(), FftError> {
    if indices.is_empty() {
        return Err(FftError::EmptyInput);
    }
    if indices.len() != values.len()
        || indices.len() != output.len()
        || indices.len() > fft_len
        || !fft_len.is_power_of_two()
        || fft_len > i32::MAX as usize
    {
        return Err(FftError::InvalidSize);
    }

    let ok = unsafe {
        blst_fr_interpolate_fft_domain(
            output.as_mut_ptr(),
            indices.as_ptr(),
            values.as_ptr(),
            indices.len(),
            fft_len,
        )
    };
    if ok {
        Ok(())
    } else {
        Err(FftError::BackendFailure)
    }
}

/// Interpolates only coefficients `ell..n` from FFT-domain samples.
///
/// `output.len()` must equal `indices.len() - ell`.
pub fn fr_interpolate_fft_domain_suffix(
    indices: &[usize],
    values: &[blst_fr],
    ell: usize,
    fft_len: usize,
    output: &mut [blst_fr],
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

    let ok = unsafe {
        blst_fr_interpolate_fft_domain_suffix(
            output.as_mut_ptr(),
            indices.as_ptr(),
            values.as_ptr(),
            indices.len(),
            ell,
            fft_len,
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
    fn size_two_matches_dft() {
        let a = fr_from_u64(11);
        let b = fr_from_u64(7);
        let mut values = [a, b];

        fr_fft(&mut values, false).unwrap();

        assert_eq!(values, [fr_add(&a, &b), fr_sub(&a, &b)]);
    }

    #[test]
    fn transpose_forward_matches_forward_full_ntt() {
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
        let mut forward = input;
        let mut transposed = input;

        fr_fft(&mut forward, false).unwrap();
        fr_fft_transpose_forward(&mut transposed).unwrap();

        assert_eq!(transposed, forward);
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
    fn interpolate_fft_domain_indices_recovers_monomial_coefficients() {
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

        let mut recovered = [fr_from_u64(0); 5];
        fr_interpolate_fft_domain(&indices, &values, fft_len, &mut recovered)
            .unwrap();

        assert_eq!(recovered, coefficients);
    }

    #[test]
    fn suffix_interpolation_matches_full_interpolation() {
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

        let mut allocating_suffix = [fr_from_u64(0); 3];
        fr_interpolate_fft_domain_suffix(
            &indices,
            &values,
            2,
            fft_len,
            &mut allocating_suffix,
        )
        .unwrap();
        assert_eq!(allocating_suffix, coefficients[2..]);

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
            fr_interpolate_fft_domain(
                &oversized_indices,
                &oversized_values,
                8,
                &mut oversized_output,
            ),
            Err(FftError::InvalidSize)
        );
        assert_eq!(
            fr_interpolate_fft_domain_suffix(
                &oversized_indices,
                &oversized_values,
                0,
                8,
                &mut oversized_output,
            ),
            Err(FftError::InvalidSize)
        );
        assert_eq!(
            fr_interpolate_fft_domain_suffix_scratch_len(9, 8),
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
