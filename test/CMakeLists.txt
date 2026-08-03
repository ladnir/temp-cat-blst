# ---- Testing ----
include(FetchContent)

# if cmake is new enough, set CMP0135 to NEW to avoid warnings from FetchContent
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


# Public top-level API tests should compile against whichever backend the
# default aliases resolve to.
set(PUBLIC_API_TEST_SOURCES
    "Backend_alias_tests.cpp"
    "BindingHeader_tests.cpp"
    "Btx_c_api_tests.cpp"
    "Btx_tests.cpp"
    "C_api_tests.cpp"
    "CompressedInfinity_tests.cpp"
    "Conformance_tests.cpp"
    "Fft_api_tests.cpp"
    "MiddleProduct_algorithm_tests.cpp"
    "MiddleProduct_c_api_tests.cpp"
    "Pfe_tests.cpp"
    "Trx_c_api_tests.cpp"
    "Trx_tests.cpp"
)

# Scalar backend tests are allowed to use scalar-only internals and helper APIs.
set(SCALAR_BACKEND_TEST_SOURCES
    "Fp_tests.cpp"
    "Fp2_tests.cpp"
    "Fp12_tests.cpp"
    "Fr_tests.cpp"
    "G1_tests.cpp"
    "G2_tests.cpp"
    "GT_tests.cpp"
    "Pairing_tests.cpp"
    "Interface_tests.cpp"
    "Fft_tests.cpp"
    "Scalar_header_include_tests.cpp"
    "scalar/Scalar_backend_tests.cpp"
    "scalar/Scalar_mult_window_tests.cpp"
)

set(AVX_VALIDATION_TEST_SOURCES)
set(AVX_TEST_SOURCES)

set(TEST_SOURCES
    ${PUBLIC_API_TEST_SOURCES}
    "Timer_tests.cpp"
)

if(NOT BLST_DEFAULT_AVX)
    list(APPEND TEST_SOURCES ${SCALAR_BACKEND_TEST_SOURCES})
endif()

if(BLST_HAVE_AVX_VALIDATION_GMP)
    list(APPEND AVX_VALIDATION_TEST_SOURCES
        "avx/runtime_constants_validation_support.cpp"
        "avx/Avx_runtime_constants_validation_tests.cpp"
    )
    list(APPEND TEST_SOURCES ${AVX_VALIDATION_TEST_SOURCES})
endif()

if(BLST_AVX512 OR BLST_AVX512_PORTABLE)
    list(APPEND AVX_TEST_SOURCES
        "avx/Avx_add_corner_cases_tests.cpp"
        "avx/Avx_gt_tests.cpp"
        "avx/Avx_msm_tests.cpp"
    )
    list(APPEND TEST_SOURCES ${AVX_TEST_SOURCES})
endif()

if(BLST_PARALLEL_MSM)
    list(APPEND TEST_SOURCES "Msm_tests.cpp")
endif()

# Add CUDA tests if CUDA is enabled
if(BLST_CUDA)
    # Enable CUDA language
    enable_language(CUDA)
    
    # Find CUDA
    find_package(CUDAToolkit REQUIRED)
    
    # Add CUDA test source
    list(APPEND TEST_SOURCES 
        "FpCuda_tests.cpp"
        "FrCuda_tests.cpp"
        "Fp2Cuda_tests.cpp"
        "Fp12Cuda_tests.cpp"
        "G1Cuda_tests.cpp"
        "G2Cuda_tests.cpp"
        "GTCuda_tests.cpp"
        "FftCuda_tests.cpp"
        "parallel_prefix_tests.cpp"
    )

endif()

# Create single test executable with both CPU and CUDA tests
add_executable(blst_tests ${TEST_SOURCES})

target_include_directories(blst_tests PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(blst_tests PRIVATE gtest_main blst_static blst_avx_headers bte)

if(BLST_HAVE_AVX_VALIDATION_GMP AND MSVC)
    add_custom_command(TARGET blst_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<$<CONFIG:Debug>:${BLST_AVX_MPIR_DEBUG_DLL}>
            $<$<NOT:$<CONFIG:Debug>>:${BLST_AVX_MPIR_RELEASE_DLL}>
            $<TARGET_FILE_DIR:blst_tests>
    )
endif()


# Add CUDA tests if CUDA is enabled
if(BLST_CUDA)
    set_source_files_properties(${TEST_SOURCES} PROPERTIES LANGUAGE CUDA)

    set_target_properties(blst_tests PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

include(GoogleTest)        # 4) glue between GTest and CTest
# 5) Auto-discover TEST() cases inside blst_tests and register them with CTest
gtest_discover_tests(blst_tests
  DISCOVERY_TIMEOUT 60
  DISCOVERY_MODE PRE_TEST
  # Optional:
  # WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  # TEST_PREFIX "unit."
  # PROPERTIES LABELS "unit"
)

if(BLST_ASAN OR BLST_UBSAN)
    # Keep the focused sanitizer change loop independent of the broad test
    # executable and the frozen BTx/BTE subsystem.
    add_executable(blst_sanitizer_boundary_tests
        "Backend_alias_tests.cpp"
        "BindingHeader_tests.cpp"
        "C_api_tests.cpp"
        "CompressedInfinity_tests.cpp"
        "Conformance_tests.cpp"
        "Fft_api_tests.cpp"
        "Interface_tests.cpp"
        "MiddleProduct_c_api_tests.cpp"
        "scalar/Scalar_mult_window_tests.cpp"
    )
    target_include_directories(blst_sanitizer_boundary_tests PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}
    )
    target_link_libraries(blst_sanitizer_boundary_tests PRIVATE
        gtest_main
        blst_static
        blst_avx_headers
    )
    gtest_discover_tests(blst_sanitizer_boundary_tests
        DISCOVERY_TIMEOUT 60
        DISCOVERY_MODE PRE_TEST
        TEST_PREFIX "sanitizer."
        PROPERTIES LABELS "sanitizer;boundary"
    )
endif()

if(BLST_AVX512 AND BLST_AVX512_PORTABLE)
    add_executable(blst_avx512_compat_tests
        "avx512_compat_test.cpp"
        "${PROJECT_SOURCE_DIR}/blst/backend/avx512/intrin_compat.cpp"
    )
    target_include_directories(blst_avx512_compat_tests PRIVATE
        "${PROJECT_SOURCE_DIR}"
        "${CMAKE_BINARY_DIR}/blst"
    )
    target_link_libraries(blst_avx512_compat_tests PRIVATE gtest_main)
    target_compile_features(blst_avx512_compat_tests PRIVATE cxx_std_23)

    if(MSVC)
        target_compile_options(blst_avx512_compat_tests PRIVATE /arch:AVX512)
    else()
        target_compile_options(blst_avx512_compat_tests PRIVATE -march=native)
    endif()

    gtest_discover_tests(blst_avx512_compat_tests
        DISCOVERY_TIMEOUT 60
        DISCOVERY_MODE PRE_TEST
        TEST_PREFIX "avx512_compat."
        PROPERTIES LABELS "avx512-compat"
    )
endif()
