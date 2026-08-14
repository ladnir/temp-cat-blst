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

add_executable(blst_fuzz_seed_encoding_tests EXCLUDE_FROM_ALL
    tests/seed_encoding_tests.cpp)
target_compile_features(blst_fuzz_seed_encoding_tests PRIVATE cxx_std_20)
target_include_directories(blst_fuzz_seed_encoding_tests PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(blst_fuzz_seed_encoding_tests PRIVATE blst_sanitizers)

function(blst_add_fuzz_target target source max_len timeout_seconds)
    set(replay_timeout_seconds ${timeout_seconds})
    set(replay_test_timeout_seconds 60)
    # Sanitizer-instrumented Debug field and pairing arithmetic is
    # intentionally much slower than the RelWithDebInfo discovery profile.
    # Preserve the campaign's per-input discovery limit, but give deterministic
    # assertion replay enough time to complete the same retained seed.
    if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND timeout_seconds GREATER 2)
        set(replay_timeout_seconds 60)
        set(replay_test_timeout_seconds 120)
    endif()
    add_executable(${target} EXCLUDE_FROM_ALL
        support/entrypoint.cpp
        ${source}
        internal/adapter.cpp
        ${ARGN})
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
            -timeout=${replay_timeout_seconds}
            ${PROJECT_SOURCE_DIR}/fuzz/corpus/${target})
    set_tests_properties(${target}_replay PROPERTIES
        LABELS "fuzz;phase1;replay"
        RUN_SERIAL TRUE
        TIMEOUT ${replay_test_timeout_seconds})
endfunction()

blst_add_fuzz_target(blst_fuzz_decode targets/decode.cpp 192 2)
blst_add_fuzz_target(blst_fuzz_bytes targets/bytes.cpp 1024 2)
blst_add_fuzz_target(
    blst_fuzz_groups targets/groups.cpp 64 2 internal/adapter_groups.cpp)
blst_add_fuzz_target(
    blst_fuzz_fields targets/fields.cpp 256 2 internal/adapter_fields.cpp)
blst_add_fuzz_target(
    blst_fuzz_protocol targets/protocol.cpp 2048 10 internal/adapter_protocol.cpp)
blst_add_fuzz_target(
    blst_fuzz_uniq targets/uniq.cpp 4097 10 internal/adapter_uniq.cpp)
blst_add_fuzz_target(
    blst_fuzz_pairing_state targets/pairing_state.cpp 1024 10
    internal/adapter_pairing.cpp)
blst_add_fuzz_target(
    blst_fuzz_fft_middle targets/fft_middle.cpp 128 10
    internal/adapter_fft_middle.cpp)
blst_add_fuzz_target(
    blst_fuzz_msm targets/msm.cpp 64 10 internal/adapter_msm.cpp)
if(BLST_FUZZ_AVX)
    blst_add_fuzz_target(
        blst_fuzz_avx targets/avx.cpp 192 10 internal/adapter_avx.cpp)
    set_tests_properties(blst_fuzz_avx_replay PROPERTIES
        LABELS "fuzz;phase4;replay;manual-avx")
endif()
set_tests_properties(
    blst_fuzz_groups_replay
    blst_fuzz_fields_replay
    PROPERTIES LABELS "fuzz;phase2;replay")
set_tests_properties(
    blst_fuzz_protocol_replay
    blst_fuzz_pairing_state_replay
    blst_fuzz_fft_middle_replay
    blst_fuzz_msm_replay
    blst_fuzz_uniq_replay
    PROPERTIES LABELS "fuzz;phase3;replay")

add_custom_target(blst_fuzz_phase0
    DEPENDS
        blst_fuzz_domain_compile
        blst_fuzz_adapter_compile
        blst_fuzz_call_frame_tests
        blst_fuzz_modulus_oracle_tests
        blst_fuzz_seed_encoding_tests
        blst_fuzz_scalar_factory_tests)
add_dependencies(blst_fuzz_phase0 blst_fuzz_adapter_factory_tests)

add_test(NAME blst_fuzz_call_frame_tests COMMAND blst_fuzz_call_frame_tests)
add_test(NAME blst_fuzz_scalar_factory_tests COMMAND blst_fuzz_scalar_factory_tests)
add_test(NAME blst_fuzz_adapter_factory_tests COMMAND blst_fuzz_adapter_factory_tests)
add_test(NAME blst_fuzz_modulus_oracle_tests COMMAND blst_fuzz_modulus_oracle_tests)
add_test(NAME blst_fuzz_seed_encoding_tests COMMAND blst_fuzz_seed_encoding_tests)
set_tests_properties(
    blst_fuzz_call_frame_tests
    blst_fuzz_scalar_factory_tests
    blst_fuzz_adapter_factory_tests
    blst_fuzz_modulus_oracle_tests
    blst_fuzz_seed_encoding_tests
    PROPERTIES LABELS "fuzz;phase0")

add_custom_target(blst_fuzz_phase1
    DEPENDS blst_fuzz_decode blst_fuzz_bytes)

add_custom_target(blst_fuzz_phase2
    DEPENDS blst_fuzz_fields blst_fuzz_groups)

add_custom_target(blst_fuzz_phase3
    DEPENDS
        blst_fuzz_fft_middle
        blst_fuzz_msm
        blst_fuzz_pairing_state
        blst_fuzz_protocol
        blst_fuzz_uniq)

if(BLST_FUZZ_AVX)
    add_custom_target(blst_fuzz_phase4_avx DEPENDS blst_fuzz_avx)
endif()
