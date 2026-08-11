# Middle Product

The Rust middle-product surface is intentionally small:

- `MiddleProductPrecomputedB`
- `precompute_middle_product_b`
- `middle_product`
- `update`
- `combine`

Design points:

- `b` is provided as scalar-ABI `&[blst_p2]`
- `a` is provided as scalar-ABI `&[blst_p1]`
- `a0` is provided as one scalar-ABI `&blst_p1`
- output is written into caller-owned `&mut [blst_fp12]`
- precomputation owns an opaque native handle
- evaluation and update currently allocate native staging/work buffers on each
  call; allocation failure is returned as an error and leaves output unchanged

Mathematical and point-domain contract:

- for `a.len() == n` and `b.len() == 2*n`, output `k` is
  `product_i pairing(a[i], b[n + k - i])`
- every `blst_p1` and `blst_p2` input must already be reduced, on curve, and in
  its prime-order subgroup; infinity is accepted as the identity
- the native implementation performs no Release point or subgroup validation
- outputs are GT elements stored in `blst_fp12`; `combine` likewise requires GT
  inputs, not arbitrary Fp12 values

Length rules:

- `precompute_middle_product_b` expects `b.len() == 2 * n`
- `middle_product` expects `a.len() == n`
- `update` writes exactly `n` outputs for the precomputed handle
- `combine` expects `mid.len() == cnst.len()`
- `out.len()` must be at least `n`
- unsupported native FFT lengths fail precomputation rather than wrapping size
  arithmetic or attempting a transform

Backend behavior:

- the Rust ABI stays scalar/public
- the native implementation may internally use AVX and convert from scalar
  boundary types when the build enables that backend
- `combine` uses the GT group law; at the `blst_fp12` boundary this appears as
  elementwise multiplication because GT is represented inside Fp12
- all algorithms are variable-time with respect to their public points and
  lengths
- the opaque handle is immutable after successful construction; the safe Rust
  wrapper owns its lifetime and does not expose concurrent mutation
