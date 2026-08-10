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
        "Fp12Coop_bench.cpp"
        "FFT_bench.cpp"
        "MiddleProduct_multipair_bench.cpp"
    )
else()
    list(APPEND BENCHMARK_SOURCES
        "main.cpp"
        "Fr_bench.cpp"
        "G1_bench.cpp"
        "G2_bench.cpp"
        "ScalarBatchAffineMsm_bench.cpp"
        "ScalarPrimitive_bench.cpp"
        "pairing_bench.cpp"
        "ScalarMsm_bench.cpp"
        "FFT_bench.cpp"
        "MiddleProduct_multipair_bench.cpp"
    )
endif()

list(APPEND BENCHMARK_SOURCES
    "Trx_bench.cpp"
    "Btx_bench.cpp"
    "BtxBeat_head_to_head_bench.cpp"
    "Pfe_bench.cpp"
    "PfeBtx_head_to_head_bench.cpp")

if(BLST_HAVE_AVX_VALIDATION_GMP)
    list(APPEND BENCHMARK_SOURCES "avx/Avx_gt_glv4_scalar_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_basefield_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_basefield_inverse_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_basefield_reduce_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_basefield_50bit_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_basefield_inverse_50bit_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_basefield_reduce_50bit_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_curve_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_group_ops_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_curve_50bit_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_scalar_mult_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_scalar_mult_50bit_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_gt_glv4_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_fp12_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_fp12_50bit_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_fp12_scalar_baseline_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_fp12_cyclotomic_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_fp12_inverse_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_fp12_inverse_50bit_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_fp12_frobenius_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_fp12_finalexp_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_fp12_finalexp_wrapper_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_fp12_finalexp_wrapper_50bit_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_fp12_scalar_cyclotomic_baseline_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_fp12_sparse_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_fp12_sparse_scalar_baseline_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_miller_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_miller_50bit_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_pairing_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_pairing_50bit_bench.cpp")
    list(APPEND BENCHMARK_SOURCES "avx/Avx_pairing_scalar_baseline_bench.cpp")
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

# Create benchmark executable
add_executable(blst_benchmarks ${BENCHMARK_SOURCES})

target_include_directories(blst_benchmarks PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(blst_benchmarks PRIVATE ${CMAKE_CURRENT_LIST_DIR}/..)
target_link_libraries(blst_benchmarks PRIVATE benchmark::benchmark blst_static blst_avx_headers bte)

add_executable(group_scalar_mult_benchmarks
    "main.cpp"
    "Group_scalar_mult_bench.cpp"
)
target_include_directories(group_scalar_mult_benchmarks PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/..
)
target_link_libraries(group_scalar_mult_benchmarks PRIVATE
    benchmark::benchmark
    blst_static
    blst_avx_headers
)

if(BLST_HAVE_AVX_VALIDATION_GMP)
    if(BLST_AVX_FIELD_PERF_ONLY)
        add_executable(avx_group_ops_benchmarks
            "avx/Avx_matrixnok_head_to_head_bench.cpp"
        )
        if(BLST_VROOM_SOURCE_DIR)
            if(NOT EXISTS "${BLST_VROOM_SOURCE_DIR}/src/bls12_381_ring_params.hpp")
                message(FATAL_ERROR
                    "BLST_VROOM_SOURCE_DIR must name a VROOM source tree containing src/bls12_381_ring_params.hpp")
            endif()
            set(BLST_VROOM_REFERENCE_COMMIT "26bd9aee2168f167583cf096add355547a7b56d4")
            execute_process(
                COMMAND git -C "${BLST_VROOM_SOURCE_DIR}" rev-parse HEAD
                RESULT_VARIABLE BLST_VROOM_GIT_RESULT
                OUTPUT_VARIABLE BLST_VROOM_ACTUAL_COMMIT
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
            if(NOT BLST_VROOM_GIT_RESULT EQUAL 0 OR
                NOT BLST_VROOM_ACTUAL_COMMIT STREQUAL BLST_VROOM_REFERENCE_COMMIT)
                message(FATAL_ERROR
                    "Focused AVX comparison requires VROOM commit ${BLST_VROOM_REFERENCE_COMMIT}; "
                    "found '${BLST_VROOM_ACTUAL_COMMIT}'")
            endif()
            target_sources(avx_group_ops_benchmarks PRIVATE
                "avx/Vroom_matrixnok_fft_preview_bench.cpp"
            )
            target_include_directories(avx_group_ops_benchmarks PRIVATE
                "${BLST_VROOM_SOURCE_DIR}/src"
            )
        endif()
        target_link_libraries(avx_group_ops_benchmarks PRIVATE
            benchmark::benchmark_main
            blst_avx_headers
            blst_static
        )
    else()
        add_executable(avx_group_ops_benchmarks
            "main.cpp"
            "avx/Avx_group_ops_bench.cpp"
            "avx/Avx_basefield_50bit_bench.cpp"
        )
        target_link_libraries(avx_group_ops_benchmarks PRIVATE
            benchmark::benchmark
            blst_static
            blst_avx_headers
            bte
        )
    endif()
    target_include_directories(avx_group_ops_benchmarks PRIVATE ${CMAKE_CURRENT_LIST_DIR})
    target_include_directories(avx_group_ops_benchmarks PRIVATE ${CMAKE_CURRENT_LIST_DIR}/..)
    target_include_directories(avx_group_ops_benchmarks PRIVATE
        "${CMAKE_BINARY_DIR}"
        "${CMAKE_BINARY_DIR}/blst"
    )
    target_compile_features(avx_group_ops_benchmarks PRIVATE cxx_std_23)
    if(BLST_AVX_FIELD_PERF_ONLY AND BLST_AVX512)
        if(MSVC)
            target_compile_options(avx_group_ops_benchmarks PRIVATE /arch:AVX512)
        else()
            target_compile_options(avx_group_ops_benchmarks PRIVATE -march=native)
        endif()
    endif()
endif()

if(BLST_AVX512 OR BLST_AVX512_PORTABLE)
    add_executable(avx_pairing_ratio_benchmarks
        "main.cpp"
        "avx/Avx_pairing_ratio_bench.cpp"
    )
    target_include_directories(avx_pairing_ratio_benchmarks PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}
        ${CMAKE_CURRENT_LIST_DIR}/..
    )
    target_link_libraries(avx_pairing_ratio_benchmarks PRIVATE
        benchmark::benchmark
        blst_static
        blst_avx_headers
    )

    add_executable(avx_gt_msm_benchmarks
        "main.cpp"
        "avx/Avx_gt_msm_bench.cpp"
    )
    target_include_directories(avx_gt_msm_benchmarks PRIVATE ${CMAKE_CURRENT_LIST_DIR})
    target_include_directories(avx_gt_msm_benchmarks PRIVATE ${CMAKE_CURRENT_LIST_DIR}/..)
    target_link_libraries(avx_gt_msm_benchmarks PRIVATE benchmark::benchmark blst_static blst_avx_headers)
endif()

add_executable(bte_trx_benchmarks
    "main.cpp"
    "Trx_bench.cpp"
)
target_include_directories(bte_trx_benchmarks PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(bte_trx_benchmarks PRIVATE ${CMAKE_CURRENT_LIST_DIR}/..)
target_link_libraries(bte_trx_benchmarks PRIVATE benchmark::benchmark blst_static blst_avx_headers bte)

add_executable(bte_btx_benchmarks
    "main.cpp"
    "Btx_bench.cpp"
    "BtxBeat_head_to_head_bench.cpp"
)
target_include_directories(bte_btx_benchmarks PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(bte_btx_benchmarks PRIVATE ${CMAKE_CURRENT_LIST_DIR}/..)
target_link_libraries(bte_btx_benchmarks PRIVATE benchmark::benchmark blst_static blst_avx_headers bte)

add_executable(bte_pfe_benchmarks
    "main.cpp"
    "Pfe_bench.cpp"
    "PfeBtx_head_to_head_bench.cpp"
)
target_include_directories(bte_pfe_benchmarks PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(bte_pfe_benchmarks PRIVATE ${CMAKE_CURRENT_LIST_DIR}/..)
target_link_libraries(bte_pfe_benchmarks PRIVATE benchmark::benchmark blst_static blst_avx_headers bte)

if(BLST_HAVE_AVX_VALIDATION_GMP AND MSVC)
    add_custom_command(TARGET blst_benchmarks POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<$<CONFIG:Debug>:${BLST_AVX_MPIR_DEBUG_DLL}>
            $<$<NOT:$<CONFIG:Debug>>:${BLST_AVX_MPIR_RELEASE_DLL}>
            $<TARGET_FILE_DIR:blst_benchmarks>
    )
endif()

# Set CUDA properties if CUDA is enabled
if(BLST_CUDA)
    set_source_files_properties(${BENCHMARK_SOURCES} PROPERTIES LANGUAGE CUDA)

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

add_custom_target(run_bte_trx_benchmarks
    COMMAND bte_trx_benchmarks
    DEPENDS bte_trx_benchmarks
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running BTE TrX benchmarks..."
)

add_custom_target(run_bte_btx_benchmarks
    COMMAND bte_btx_benchmarks
    DEPENDS bte_btx_benchmarks
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running BTE BTX benchmarks..."
)

add_custom_target(run_bte_pfe_benchmarks
    COMMAND bte_pfe_benchmarks
    DEPENDS bte_pfe_benchmarks
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running BTE PFE benchmarks..."
)
