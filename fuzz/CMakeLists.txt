# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

add_library(blst_fuzz_domain_compile OBJECT EXCLUDE_FROM_ALL
    tests/domain_compile_tests.cpp)
target_compile_features(blst_fuzz_domain_compile PRIVATE cxx_std_20)
target_include_directories(blst_fuzz_domain_compile PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(blst_fuzz_domain_compile PRIVATE blst_sanitizers)

add_library(blst_fuzz_adapter_compile OBJECT EXCLUDE_FROM_ALL
    internal/adapter.cpp)
target_compile_features(blst_fuzz_adapter_compile PRIVATE cxx_std_20)
target_include_directories(blst_fuzz_adapter_compile PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_BINARY_DIR}
    ${PROJECT_BINARY_DIR}/blst)
target_link_libraries(blst_fuzz_adapter_compile PRIVATE blst_sanitizers)

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

add_executable(blst_fuzz_adapter_factory_tests EXCLUDE_FROM_ALL
    tests/adapter_factory_tests.cpp
    internal/adapter.cpp)
target_compile_features(blst_fuzz_adapter_factory_tests PRIVATE cxx_std_20)
target_include_directories(blst_fuzz_adapter_factory_tests PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_BINARY_DIR}
    ${PROJECT_BINARY_DIR}/blst)
target_link_libraries(blst_fuzz_adapter_factory_tests PRIVATE
    blst_static
    blst_sanitizers)

add_executable(blst_fuzz_modulus_oracle_tests EXCLUDE_FROM_ALL
    tests/modulus_oracle_tests.cpp)
target_compile_features(blst_fuzz_modulus_oracle_tests PRIVATE cxx_std_20)
target_include_directories(blst_fuzz_modulus_oracle_tests PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(blst_fuzz_modulus_oracle_tests PRIVATE blst_sanitizers)

function(blst_add_fuzz_target target source max_len timeout_seconds)
    add_executable(${target} EXCLUDE_FROM_ALL
        support/entrypoint.cpp
        ${source}
        internal/adapter.cpp)
    target_compile_features(${target} PRIVATE cxx_std_20)
    target_include_directories(${target} PRIVATE
        ${PROJECT_SOURCE_DIR}
        ${PROJECT_BINARY_DIR}
        ${PROJECT_BINARY_DIR}/blst)
    target_compile_options(${target} PRIVATE -fsanitize=fuzzer-no-link)
    target_link_options(${target} PRIVATE -fsanitize=fuzzer)
    target_link_libraries(${target} PRIVATE blst_static blst_sanitizers)
    add_test(
        NAME ${target}_replay
        COMMAND ${target}
            -runs=0
            -max_len=${max_len}
            -timeout=${timeout_seconds}
            ${PROJECT_SOURCE_DIR}/fuzz/corpus/${target})
    set_tests_properties(${target}_replay PROPERTIES
        LABELS "fuzz;phase1;replay"
        RUN_SERIAL TRUE
        TIMEOUT 60)
endfunction()

blst_add_fuzz_target(blst_fuzz_decode targets/decode.cpp 192 2)
blst_add_fuzz_target(blst_fuzz_bytes targets/bytes.cpp 1024 2)

add_custom_target(blst_fuzz_phase0
    DEPENDS
        blst_fuzz_domain_compile
        blst_fuzz_adapter_compile
        blst_fuzz_call_frame_tests
        blst_fuzz_modulus_oracle_tests
        blst_fuzz_scalar_factory_tests)
add_dependencies(blst_fuzz_phase0 blst_fuzz_adapter_factory_tests)

add_test(NAME blst_fuzz_call_frame_tests COMMAND blst_fuzz_call_frame_tests)
add_test(NAME blst_fuzz_scalar_factory_tests COMMAND blst_fuzz_scalar_factory_tests)
add_test(NAME blst_fuzz_adapter_factory_tests COMMAND blst_fuzz_adapter_factory_tests)
add_test(NAME blst_fuzz_modulus_oracle_tests COMMAND blst_fuzz_modulus_oracle_tests)
set_tests_properties(
    blst_fuzz_call_frame_tests
    blst_fuzz_scalar_factory_tests
    blst_fuzz_adapter_factory_tests
    blst_fuzz_modulus_oracle_tests
    PROPERTIES LABELS "fuzz;phase0")

add_custom_target(blst_fuzz_phase1
    DEPENDS blst_fuzz_decode blst_fuzz_bytes)
