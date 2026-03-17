# ---- Benchmarking ----
include(FetchContent)

if(UNIX AND BLST_AVX512)
    enable_language(ASM)
endif()

# if cmake is new enough, set CMP0135 to NEW to avoid warnings from FetchContent
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

# Fetch Google Benchmark
FetchContent_Declare(
  googlebenchmark
  URL https://github.com/google/benchmark/archive/refs/tags/v1.8.3.zip
)

# Configure benchmark options
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googlebenchmark)

# Determine benchmark sources based on CUDA availability.
# CUDA-heavy benchmarks pull in thrust/cuda headers and should not be built
# when BLST_CUDA is disabled.
set(BENCHMARK_SOURCES)

if(BLST_CUDA)
    list(APPEND BENCHMARK_SOURCES
        "Fp12_bench.cpp"
        "main.cpp"
        "avxbls_bench.cpp"
        "Fp12Coop_bench.cpp"
        "FFT_bench.cpp"
    )
else()
    list(APPEND BENCHMARK_SOURCES
        "main.cpp"
        "Fp_avx_lufp_v2_bench.cpp"
        "Fp_avx_renorm_bench.cpp"
        "Fp_avx_k_bench.cpp"
        "Fp_avx_fused_baseline_bench.cpp"
        "Fr_bench.cpp"
        "G1_bench.cpp"
        "G2_bench.cpp"
        "pairing_bench.cpp"
        "FFT_bench.cpp"
    )
endif()

if(BLST_HAVE_VROOM_BASEFIELD)
    list(APPEND BENCHMARK_SOURCES "../blst/vroom/bench/Vroom_basefield_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "../blst/vroom/bench/Vroom_fp12_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "Vroom_fp12_scalar_baseline_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "../blst/vroom/bench/Vroom_fp12_cyclotomic_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "Vroom_fp12_scalar_cyclotomic_baseline_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "../blst/vroom/bench/Vroom_fp12_sparse_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "Vroom_fp12_sparse_scalar_baseline_bench.cpp")
endif()

if(BLST_PARALLEL_MSM)
    list(APPEND BENCHMARK_SOURCES "MSM_bench.cpp")
endif()

# Add CUDA benchmarks if CUDA is enabled
#if(BLST_CUDA)
#    # Enable CUDA language
#    enable_language(CUDA)
#    
#    # Find CUDA
#    find_package(CUDAToolkit REQUIRED)
#    
#    # Add CUDA benchmark sources
#    list(APPEND BENCHMARK_SOURCES 
#        "FpCuda_bench.cu"
#        "FrCuda_bench.cu"
#        "Fp2Cuda_bench.cu"
#        "Fp12Cuda_bench.cu"
#        "G1Cuda_bench.cu"
#        "G2Cuda_bench.cu"
#    )
#endif()

# Native AVX-512 benchmark-only fused k=1 mul+redc core ASM kernel.
if(UNIX AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64)$" AND BLST_ASM AND BLST_AVX512 AND NOT BLST_AVX512_PORTABLE)
    list(APPEND BENCHMARK_SOURCES "Fp_avx_fused_core_k1_asm.S")
    list(APPEND BENCHMARK_SOURCES "Fp_avx_mul_v2_k2_v6b_asm.S")
endif()

# Create benchmark executable
add_executable(blst_benchmarks ${BENCHMARK_SOURCES})

target_include_directories(blst_benchmarks PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(blst_benchmarks PRIVATE benchmark::benchmark blst_static blst_vroom_headers)

if(BLST_HAVE_VROOM_BASEFIELD AND MSVC)
    add_custom_command(TARGET blst_benchmarks POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<$<CONFIG:Debug>:${BLST_VROOM_MPIR_DEBUG_DLL}>
            $<$<NOT:$<CONFIG:Debug>>:${BLST_VROOM_MPIR_RELEASE_DLL}>
            $<TARGET_FILE_DIR:blst_benchmarks>
    )
endif()

# Set CUDA properties if CUDA is enabled
if(BLST_CUDA)
    set_source_files_properties(${BENCHMARK_SOURCES} PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(avxbls_bench.cpp PROPERTIES LANGUAGE CXX)

    set_target_properties(blst_benchmarks PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    target_compile_options(blst_benchmarks PRIVATE -linetable -Xptxas=-O3  -Xptxas=-v ) #-maxrregcount=64
    
endif()

# Add a custom target to run benchmarks
add_custom_target(run_benchmarks
    COMMAND blst_benchmarks
    DEPENDS blst_benchmarks
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running benchmarks..."
)
