# ---- Testing ----
include(FetchContent)

# if cmake is new enough, set CMP0135 to NEW to avoid warnings from FetchContent
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

find_package(GTest 1.14 QUIET)
if(NOT TARGET GTest::gtest_main)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
        URL_HASH SHA256=1f357c27ca988c3f7c6b4bf68a9395005ac6761f034046e9dde0896e3aba00e4
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()


# Public top-level API tests should compile against whichever backend the
# default aliases resolve to.
set(PUBLIC_API_TEST_SOURCES
    "Backend_alias_tests.cpp"
    "BindingHeader_tests.cpp"
    "Btx_c_api_tests.cpp"
    "Btx_conformance_vector_tests.cpp"
    "Btx_tests.cpp"
    "C_api_tests.cpp"
    "CompressedInfinity_tests.cpp"
    "Conformance_tests.cpp"
    "Fft_api_tests.cpp"
    "MiddleProduct_algorithm_tests.cpp"
    "MiddleProduct_c_api_tests.cpp"
    "Pfe_tests.cpp"
    "Transparent_bte_tests.cpp"
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

# Audit regressions are ordinary unit tests. Keep them in the main test binary
# so each configuration has one discoverable runtime suite.
set(AUDIT_TEST_SOURCES
    "Addchain_audit_tests.cpp"
    "Aggregate_audit_tests.cpp"
    "Assembly_arithmetic_audit_tests.cpp"
    "Group_scalar_mult_tests.cpp"
    "Gt_deserialize_audit_tests.cpp"
    "Hash_to_field_audit_tests.cpp"
    "Keygen_audit_tests.cpp"
    "Map_to_group_audit_tests.cpp"
    "Scalar_point_audit_tests.cpp"
    "Sha256_audit_tests.cpp"
)

set(AVX_TEST_SOURCES)

set(TEST_SOURCES
    ${PUBLIC_API_TEST_SOURCES}
    ${AUDIT_TEST_SOURCES}
    "Timer_tests.cpp"
)

if(NOT BLST_DEFAULT_AVX)
    list(APPEND TEST_SOURCES ${SCALAR_BACKEND_TEST_SOURCES})
endif()

if(BLST_AVX512 OR BLST_AVX512_PORTABLE)
    list(APPEND AVX_TEST_SOURCES
        "avx/Avx_add_corner_cases_tests.cpp"
        "avx/Avx_gt_tests.cpp"
        "avx/Avx_matrixnok_primitive_tests.cpp"
        "avx/Avx_msm_tests.cpp"
    )
    list(APPEND TEST_SOURCES ${AVX_TEST_SOURCES})

    list(APPEND TEST_SOURCES
        "avx/Avx_c_api_contract_tests.cpp"
    )
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

# Sanitizer presets use the same test executable with a deliberately narrow
# source set. This avoids compiling frozen BTx/BTE and unrelated broad suites
# in the sanitizer change loop without creating a second runtime test binary.
if(BLST_ASAN OR BLST_UBSAN)
    set(TEST_SOURCES
        "Backend_alias_tests.cpp"
        "BindingHeader_tests.cpp"
        "C_api_tests.cpp"
        "CompressedInfinity_tests.cpp"
        "Conformance_tests.cpp"
        "Fft_api_tests.cpp"
        "Interface_tests.cpp"
        "MiddleProduct_algorithm_tests.cpp"
        "MiddleProduct_c_api_tests.cpp"
        "scalar/Scalar_mult_window_tests.cpp"
    )
    if(BLST_AVX512 OR BLST_AVX512_PORTABLE)
        list(APPEND TEST_SOURCES
            "avx/Avx_c_api_contract_tests.cpp"
            "avx/Avx_msm_tests.cpp"
        )
    endif()
endif()

# Create the single ordinary runtime test executable.
add_executable(blst_tests ${TEST_SOURCES})

target_include_directories(blst_tests PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(blst_tests PRIVATE GTest::gtest_main blst_static blst_avx_headers)
if(NOT BLST_ASAN AND NOT BLST_UBSAN)
    target_link_libraries(blst_tests PRIVATE bte)
endif()

# A standalone allocation-failure probe keeps global operator-new replacement
# out of the main GTest executable. Construction may allocate; native scalar
# and AVX evaluation must succeed after every allocation route is disabled.
add_executable(middle_product_workspace_allocation_test
    "middle_product_workspace_allocation_test.cpp"
)
target_include_directories(middle_product_workspace_allocation_test PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${CMAKE_BINARY_DIR}"
    "${CMAKE_BINARY_DIR}/blst"
)
target_link_libraries(middle_product_workspace_allocation_test PRIVATE
    blst_static
    blst_avx_headers
)
add_test(
    NAME unit.middle_product_workspace_allocation
    COMMAND middle_product_workspace_allocation_test
)
set_tests_properties(unit.middle_product_workspace_allocation PROPERTIES
    LABELS "unit;fft;allocation"
)

# Add CUDA tests if CUDA is enabled
if(BLST_CUDA)
    set_source_files_properties(${TEST_SOURCES} PROPERTIES LANGUAGE CUDA)

    set_target_properties(blst_tests PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

include(GoogleTest)        # 4) glue between GTest and CTest

# Fixed AVX parameters are security-sensitive generated inputs. Validate them
# independently without creating another test executable.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_test(NAME audit.matrixnok_constants
        COMMAND ${Python3_EXECUTABLE}
            "${PROJECT_SOURCE_DIR}/scripts/audit/validate_matrixnok_constants.py"
            --self-test
    )
    set_tests_properties(audit.matrixnok_constants PROPERTIES LABELS "audit;avx512")

    add_test(NAME audit.matrixnok_curve_constants
        COMMAND ${Python3_EXECUTABLE}
            "${PROJECT_SOURCE_DIR}/scripts/audit/validate_matrixnok_curve_constants.py"
            --self-test
    )
    set_tests_properties(audit.matrixnok_curve_constants PROPERTIES LABELS "audit;avx512")

    add_test(NAME audit.matrixnok_lean_constants
        COMMAND ${Python3_EXECUTABLE}
            "${PROJECT_SOURCE_DIR}/scripts/audit/generate_matrixnok_lean.py"
            --check
    )
    set_tests_properties(audit.matrixnok_lean_constants PROPERTIES LABELS "audit;avx512;formal_proof")

    add_test(NAME audit.matrixnok_dispatch
        COMMAND ${Python3_EXECUTABLE}
            "${PROJECT_SOURCE_DIR}/scripts/audit/validate_matrixnok_dispatch.py"
    )
    set_tests_properties(audit.matrixnok_dispatch PROPERTIES LABELS "audit;avx512;formal_proof")

    add_test(NAME audit.avx_proof_correspondence
        COMMAND ${Python3_EXECUTABLE}
            "${PROJECT_SOURCE_DIR}/scripts/audit/validate_avx_correspondence.py"
    )
    set_tests_properties(audit.avx_proof_correspondence PROPERTIES LABELS "audit;avx512;formal_proof")
else()
    message(WARNING "Python 3 not found; MatrixNoK constants, Lean binding, dispatch, and proof correspondence will not be validated by CTest")
endif()

# 5) Auto-discover TEST() cases inside blst_tests and register them with CTest
gtest_discover_tests(blst_tests
  DISCOVERY_TIMEOUT 60
  DISCOVERY_MODE PRE_TEST
  # Optional:
  # WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  # TEST_PREFIX "unit."
  # PROPERTIES LABELS "unit"
)

if(BLST_AVX512 AND BLST_AVX512_PORTABLE)
    set(AVX512_COMPAT_TEST_SOURCES
        "avx/Avx_matrixnok_primitive_tests.cpp"
        "avx/Avx_c_api_contract_tests.cpp"
        "${PROJECT_SOURCE_DIR}/blst/backend/backend_portable.cpp"
        "${PROJECT_SOURCE_DIR}/blst/scalar/fields/Fp.cpp"
        "${PROJECT_SOURCE_DIR}/blst/scalar/fields/Fp12.cpp"
        "${PROJECT_SOURCE_DIR}/blst/scalar/curve/GT.cpp"
        "${PROJECT_SOURCE_DIR}/blst/scalar/curve/point.cpp"
    )

    list(APPEND AVX512_COMPAT_TEST_SOURCES
        "${PROJECT_SOURCE_DIR}/blst/backend/avx512/intrin_compat.cpp"
        "avx512_compat_test.cpp"
    )

    add_executable(blst_avx512_compat_tests ${AVX512_COMPAT_TEST_SOURCES})
    target_include_directories(blst_avx512_compat_tests PRIVATE
        "${PROJECT_SOURCE_DIR}"
        "${CMAKE_BINARY_DIR}"
        "${CMAKE_BINARY_DIR}/blst"
    )
    # MatrixNoK differentials compile the scalar field, curve, and pairing
    # owners above as their independent oracle. The static archive supplies
    # only the unresolved low-level backend symbols; its duplicate scalar
    # objects are not extracted by the linker.
    target_link_libraries(blst_avx512_compat_tests PRIVATE
        GTest::gtest_main
        blst_static
        blst_avx_headers)
    target_compile_features(blst_avx512_compat_tests PRIVATE cxx_std_23)

    gtest_discover_tests(blst_avx512_compat_tests
        DISCOVERY_TIMEOUT 60
        DISCOVERY_MODE PRE_TEST
        TEST_PREFIX "avx512_compat."
        PROPERTIES LABELS "avx512-compat"
    )
endif()
