# src

Main C++ code (templates/headers, tests, benchmarks). Build with `make`; see repo root **README.md** for targets.

- **`*.hpp`** — RNS Montgomery, pairing, EC, scalar mult.
- **`test_*.cpp`** — unit tests.
- **`bench_*.cpp`** — Google Benchmark harnesses.
- **`gmp_wrapper.hpp`** — pulls in **`../cpu/precompute/gmp_wrapper.hpp`** and **`../test_data/test_values.hpp`**.
