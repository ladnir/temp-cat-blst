# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

add_library(blst_fuzz_domain_compile OBJECT EXCLUDE_FROM_ALL
    tests/domain_compile_tests.cpp)
target_compile_features(blst_fuzz_domain_compile PRIVATE cxx_std_20)
target_include_directories(blst_fuzz_domain_compile PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(blst_fuzz_domain_compile PRIVATE blst_sanitizers)

add_executable(blst_fuzz_call_frame_tests EXCLUDE_FROM_ALL
    tests/call_frame_tests.cpp)
target_compile_features(blst_fuzz_call_frame_tests PRIVATE cxx_std_20)
target_include_directories(blst_fuzz_call_frame_tests PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(blst_fuzz_call_frame_tests PRIVATE blst_sanitizers)

add_executable(blst_fuzz_scalar_factory_tests EXCLUDE_FROM_ALL
    tests/scalar_factory_tests.cpp)
target_compile_features(blst_fuzz_scalar_factory_tests PRIVATE cxx_std_20)
target_include_directories(blst_fuzz_scalar_factory_tests PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(blst_fuzz_scalar_factory_tests PRIVATE blst_sanitizers)

add_custom_target(blst_fuzz_phase0
    DEPENDS
        blst_fuzz_domain_compile
        blst_fuzz_call_frame_tests
        blst_fuzz_scalar_factory_tests)

add_test(NAME blst_fuzz_call_frame_tests COMMAND blst_fuzz_call_frame_tests)
add_test(NAME blst_fuzz_scalar_factory_tests COMMAND blst_fuzz_scalar_factory_tests)
set_tests_properties(
    blst_fuzz_call_frame_tests
    blst_fuzz_scalar_factory_tests
    PROPERTIES LABELS "fuzz;phase0")
