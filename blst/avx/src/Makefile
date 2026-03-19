# Compiler selection: use CXX=g++ or CXX=clang++ (default: g++)
CXX ?= g++
INCLUDES = -I..

# Suppress warnings from BLST/vector code: static decls without defs, ignored attributes on __m512i, unknown pragmas (clang-only)
CXXFLAGS_WARN_SUPPRESS = -Wno-unused-function -Wno-ignored-attributes -Wno-unknown-pragmas
# Compiler-specific flags
ifneq ($(findstring clang++,$(CXX)),)
    CXXFLAGS = -std=c++20 -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) -DUSE_AVX512_FALLBACK
    CXXFLAGS_BENCH = -std=c++20 -O3 -march=native -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) -mavx512ifma
    CXXFLAGS_BENCH_FALLBACK = -std=c++20 -O3 -march=native -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) -DUSE_AVX512_FALLBACK
    CXXFLAGS_TEST_AVX = -std=c++20 -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) -mavx512ifma
    CXXFLAGS_TEST_FALLBACK = -std=c++20 -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) -DUSE_AVX512_FALLBACK
else
    CXXFLAGS = -std=c++2a -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) -DUSE_AVX512_FALLBACK
    CXXFLAGS_BENCH = -std=c++2a -O3 -march=native -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) -mavx512ifma
    CXXFLAGS_BENCH_FALLBACK = -std=c++2a -O3 -march=native -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) -DUSE_AVX512_FALLBACK
    CXXFLAGS_TEST_AVX = -std=c++2a -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) -mavx512ifma
    CXXFLAGS_TEST_FALLBACK = -std=c++2a -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) -DUSE_AVX512_FALLBACK
endif
TARGET_REDUCE = test_elementwise_reduce
SOURCE_REDUCE = test_elementwise_reduce.cpp
TARGET_INVERSION = test_inversion
SOURCE_INVERSION = test_inversion.cpp
TARGET_TEST_PAIRING = test_pairing
TARGET_TEST_PAIRING_AVX = test_pairing_avx
SOURCE_TEST_PAIRING = test_pairing.cpp
TARGET_TEST_EC = test_ec
SOURCE_TEST_EC = test_ec.cpp
TARGET_TEST_SCALAR_MULT = test_scalar_mult
SOURCE_TEST_SCALAR_MULT = test_scalar_mult.cpp
TARGET_BENCH_PAIRING = bench_pairing
TARGET_BENCH_PAIRING_FALLBACK = bench_pairing_fallback
SOURCE_BENCH_PAIRING = bench_pairing.cpp
TARGET_BENCH_PAIRING_50BIT = bench_pairing_50bit
TARGET_BENCH_PAIRING_50BIT_FALLBACK = bench_pairing_50bit_fallback
SOURCE_BENCH_PAIRING_50BIT = bench_pairing_50bit.cpp

# All header files in v4 directory
HEADERS = bounds.hpp constexpr_utils.hpp bounded_element.hpp ring_element.hpp \
          elementwise_reduce.hpp change_base.hpp bounded_ring.hpp \
          fp12.hpp final_exponentiation.hpp miller.hpp fp2.hpp inversion.hpp \
          ec.hpp conversion_inversion.hpp scalar_mult.hpp

# Google Benchmark library (adjust path as needed)
BENCHMARK_LIBS = -lbenchmark -lpthread
# BENCHMARK_LIBS = /home/ec2-user/benchmark/build2/src/libbenchmark.a -lpthread

.PHONY: all clean test test-reduce test-inversion test-pairing test-pairing-avx test-ec test-scalar-mult bench-pairing bench-pairing-fallback bench-pairing-50bit bench-pairing-50bit-fallback run

all: $(TARGET_REDUCE) $(TARGET_INVERSION) $(TARGET_TEST_PAIRING) $(TARGET_TEST_PAIRING_AVX) $(TARGET_TEST_EC) $(TARGET_TEST_SCALAR_MULT) $(TARGET_BENCH_PAIRING) $(TARGET_BENCH_PAIRING_FALLBACK) $(TARGET_BENCH_PAIRING_50BIT) $(TARGET_BENCH_PAIRING_50BIT_FALLBACK)

$(TARGET_REDUCE): $(SOURCE_REDUCE) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SOURCE_REDUCE) -o $(TARGET_REDUCE)

$(TARGET_INVERSION): $(SOURCE_INVERSION) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SOURCE_INVERSION) -o $(TARGET_INVERSION) -lgmp -lgmpxx

$(TARGET_TEST_PAIRING): $(SOURCE_TEST_PAIRING) $(HEADERS)
	$(CXX) $(CXXFLAGS) -D__ADX__ $(INCLUDES) $(SOURCE_TEST_PAIRING) -o $(TARGET_TEST_PAIRING) -lgmp -lgmpxx ../blst/libblst_small.a

$(TARGET_TEST_PAIRING_AVX): $(SOURCE_TEST_PAIRING) $(HEADERS)
	$(CXX) $(CXXFLAGS_TEST_AVX) -D__ADX__ $(INCLUDES) $(SOURCE_TEST_PAIRING) -o $(TARGET_TEST_PAIRING_AVX) -lgmp -lgmpxx ../blst/libblst_small.a

$(TARGET_TEST_EC): $(SOURCE_TEST_EC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -D__ADX__ $(INCLUDES) $(SOURCE_TEST_EC) -o $(TARGET_TEST_EC) -lgmp -lgmpxx ../blst/libblst_small.a

$(TARGET_TEST_SCALAR_MULT): $(SOURCE_TEST_SCALAR_MULT) $(HEADERS)
	$(CXX) $(CXXFLAGS_TEST_FALLBACK) -D__ADX__ $(INCLUDES) $(SOURCE_TEST_SCALAR_MULT) -o $(TARGET_TEST_SCALAR_MULT) -lgmp -lgmpxx ../blst/libblst_small.a

$(TARGET_BENCH_PAIRING): $(SOURCE_BENCH_PAIRING) $(HEADERS)
	$(CXX) $(CXXFLAGS_BENCH) -D__ADX__ $(INCLUDES) $(SOURCE_BENCH_PAIRING) -o $(TARGET_BENCH_PAIRING) -lgmp -lgmpxx ../blst/libblst_small.a $(BENCHMARK_LIBS)

$(TARGET_BENCH_PAIRING_FALLBACK): $(SOURCE_BENCH_PAIRING) $(HEADERS)
	$(CXX) $(CXXFLAGS_BENCH_FALLBACK) -D__ADX__ $(INCLUDES) $(SOURCE_BENCH_PAIRING) -o $(TARGET_BENCH_PAIRING_FALLBACK) -lgmp -lgmpxx ../blst/libblst_small.a $(BENCHMARK_LIBS)

$(TARGET_BENCH_PAIRING_50BIT): $(SOURCE_BENCH_PAIRING_50BIT) $(HEADERS)
	$(CXX) $(CXXFLAGS_BENCH) -D__ADX__ $(INCLUDES) $(SOURCE_BENCH_PAIRING_50BIT) -o $(TARGET_BENCH_PAIRING_50BIT) -lgmp -lgmpxx ../blst/libblst_small.a $(BENCHMARK_LIBS)

$(TARGET_BENCH_PAIRING_50BIT_FALLBACK): $(SOURCE_BENCH_PAIRING_50BIT) $(HEADERS)
	$(CXX) $(CXXFLAGS_BENCH_FALLBACK) -D__ADX__ $(INCLUDES) $(SOURCE_BENCH_PAIRING_50BIT) -o $(TARGET_BENCH_PAIRING_50BIT_FALLBACK) -lgmp -lgmpxx ../blst/libblst_small.a $(BENCHMARK_LIBS)

test: $(TARGET_REDUCE) $(TARGET_INVERSION)
	./$(TARGET_REDUCE)
	./$(TARGET_INVERSION)

test-reduce: $(TARGET_REDUCE)
	./$(TARGET_REDUCE)

test-inversion: $(TARGET_INVERSION)
	./$(TARGET_INVERSION)

test-pairing: $(TARGET_TEST_PAIRING)
	./$(TARGET_TEST_PAIRING)

test-pairing-avx: $(TARGET_TEST_PAIRING_AVX)
	./$(TARGET_TEST_PAIRING_AVX)

test-ec: $(TARGET_TEST_EC)
	./$(TARGET_TEST_EC)

test-scalar-mult: $(TARGET_TEST_SCALAR_MULT)
	./$(TARGET_TEST_SCALAR_MULT)

bench-pairing: $(TARGET_BENCH_PAIRING)
	./$(TARGET_BENCH_PAIRING)

bench-pairing-fallback: $(TARGET_BENCH_PAIRING_FALLBACK)
	./$(TARGET_BENCH_PAIRING_FALLBACK)

bench-pairing-50bit: $(TARGET_BENCH_PAIRING_50BIT)
	./$(TARGET_BENCH_PAIRING_50BIT)

bench-pairing-50bit-fallback: $(TARGET_BENCH_PAIRING_50BIT_FALLBACK)
	./$(TARGET_BENCH_PAIRING_50BIT_FALLBACK)

run: test

clean:
	rm -f $(TARGET_REDUCE) $(TARGET_INVERSION) $(TARGET_TEST_PAIRING) $(TARGET_TEST_PAIRING_AVX) $(TARGET_TEST_EC) $(TARGET_TEST_SCALAR_MULT) $(TARGET_BENCH_PAIRING) $(TARGET_BENCH_PAIRING_FALLBACK) $(TARGET_BENCH_PAIRING_50BIT) $(TARGET_BENCH_PAIRING_50BIT_FALLBACK)
