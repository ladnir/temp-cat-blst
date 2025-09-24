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

enable_testing()

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
)

# Add CUDA tests if CUDA is enabled
if(BLST_CUDA)
    # Enable CUDA language
    enable_language(CUDA)
    
    # Find CUDA
    find_package(CUDAToolkit REQUIRED)
    
    # Add CUDA test source
    list(APPEND TEST_SOURCES 
        "FpCuda_tests.cu"
        "FrCuda_tests.cu"

        "Fp2Cuda_tests.cu"
        "Fp12Cuda_tests.cu"
        "G1Cuda_tests.cu"
        "G2Cuda_tests.cu"
        "GTCuda_tests.cu"
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
    set_target_properties(blst_tests PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()