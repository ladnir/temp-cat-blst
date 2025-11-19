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


# Determine test sources based on CUDA availability
set(TEST_SOURCES
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
    "Msm_tests.cpp"
    "parallel_prefix_tests.cpp"
)

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
    )

endif()

# Create single test executable with both CPU and CUDA tests
add_executable(blst_tests ${TEST_SOURCES})

target_include_directories(blst_tests PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(blst_tests PRIVATE gtest_main blst_static)


# Add test with appropriate labels and properties
add_test(NAME blst_tests COMMAND blst_tests)


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