# Middle Product

The Rust middle-product surface is intentionally small:

- `MiddleProductPrecomputedB`
- `precompute_middle_product_b`
- `middle_product`

Design points:

- `b` is provided as scalar-ABI `&[blst_p2]`
- `a` is provided as scalar-ABI `&[blst_p1]`
- output is written into caller-owned `&mut [blst_fp12]`
- the only native allocation is the opaque precomputed-`b` handle

Length rules:

- `precompute_middle_product_b` expects `b.len() == 2 * n`
- `middle_product` expects `a.len() == n`
- `out.len()` must be at least `n`

Backend behavior:

- the Rust ABI stays scalar/public
- the native implementation may internally use AVX and convert from scalar
  boundary types when the build enables that backend
