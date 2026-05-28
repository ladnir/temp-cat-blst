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
        || !fft_len.is_power_of_two()
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

fn check_fft_input(values: &[blst_fr]) -> Result<(), FftError> {
    if values.is_empty() {
        return Err(FftError::EmptyInput);
    }
    if !values.len().is_power_of_two() {
        return Err(FftError::InvalidSize);
    }
    assert!(
        values.len() <= i32::MAX as usize,
        "blst Fr FFT length exceeds the C++ reference implementation limit"
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
