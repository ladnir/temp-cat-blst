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
    "Avx_dense_transform_tests.cpp"
    "BindingHeader_tests.cpp"
    "C_api_tests.cpp"
    "MiddleProduct_c_api_tests.cpp"
)

# Scalar backend tests are allowed to use scalar-only internals and helper APIs.
set(SCALAR_BACKEND_TEST_SOURCES
    "Fp_tests.cpp"
    "Fp2_tests.cpp"
    "Fp12_tests.cpp"
    "Fr_tests.cpp"
    "BindingCpp_tests.cpp"
    "G1_tests.cpp"
    "G2_tests.cpp"
    "GT_tests.cpp"
    "Pairing_tests.cpp"
    "Interface_tests.cpp"
    "Fft_tests.cpp"
    "scalar/Scalar_backend_tests.cpp"
    "avx512_compat_test.cpp"
)

# AVX backend tests should exercise the explicit blst::avx surface.
set(AVX_BACKEND_TEST_SOURCES)
set(AVX_VALIDATION_TEST_SOURCES)

set(TEST_SOURCES
    ${PUBLIC_API_TEST_SOURCES}
)

if(NOT BLST_DEFAULT_AVX)
    list(APPEND TEST_SOURCES ${SCALAR_BACKEND_TEST_SOURCES})
endif()

if(BLST_HAVE_AVX_BASEFIELD AND NOT BLST_DEFAULT_AVX)
    list(APPEND TEST_SOURCES ${AVX_BACKEND_TEST_SOURCES})
endif()

if(BLST_HAVE_AVX_VALIDATION_GMP)
    list(APPEND AVX_VALIDATION_TEST_SOURCES
        "avx/runtime_constants_validation_support.cpp"
        "avx/Avx_runtime_constants_validation_tests.cpp"
    )
    list(APPEND TEST_SOURCES ${AVX_VALIDATION_TEST_SOURCES})
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
target_link_libraries(blst_tests PRIVATE gtest_main blst_static blst_avx_headers)

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
  # Optional:
  # WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  # TEST_PREFIX "unit."
  # PROPERTIES LABELS "unit"
)
