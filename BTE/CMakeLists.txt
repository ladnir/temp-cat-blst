include(FetchContent)

# The upstream AMD64 assembly source classification conflicts with this
# project's already-enabled MASM language under current CMake/MSVC. The
# upstream intrinsics backend retains runtime SIMD dispatch and avoids that
# generator-specific integration failure.
if(BLST_TARGET_IS_WASM)
    # Emscripten identifies its scalar target as x86 for compatibility. Do not
    # let upstream infer host SSE/AVX sources for a wasm32 object archive.
    set(BLAKE3_SIMD_TYPE "none" CACHE STRING "" FORCE)
elseif(MSVC AND BLST_TARGET_IS_X86_64)
    set(BLAKE3_SIMD_TYPE "x86-intrinsics" CACHE STRING "" FORCE)
endif()

FetchContent_Declare(
    blake3_upstream
    URL https://static.crates.io/crates/blake3/blake3-1.8.2.crate
    URL_HASH SHA256=3888aaa89e4b2a40fca9848e400f6a658a5a3978de7be858e209cafa8be9a4a0
    SOURCE_SUBDIR c
)
FetchContent_MakeAvailable(blake3_upstream)

set(BTE_SOURCES
    "building_blocks/blst_cpp.cpp"
    "building_blocks/bls_ots.cpp"
    "building_blocks/hkdf_sha256.cpp"
    "building_blocks/kzg.cpp"
    "building_blocks/lagrange.cpp"
    "building_blocks/polynomial.cpp"
    "building_blocks/serialization.cpp"
    "trx/core.cpp"
    "trx/c_api.cpp"
    "btx/core.cpp"
    "btx/admission.cpp"
    "btx/c_api.cpp"
    "btx/serialization.cpp"
    "pfe/core.cpp"
)

add_library(bte STATIC ${BTE_SOURCES})
target_include_directories(bte PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
target_link_libraries(bte PUBLIC blst_static PRIVATE BLAKE3::blake3)

target_compile_options(bte PUBLIC
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:MSVC>>:/std:c++latest>
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<NOT:$<CXX_COMPILER_ID:MSVC>>>:-std=c++23>
)

if(TARGET blst_sanitizers)
    target_link_libraries(bte PUBLIC blst_sanitizers)
endif()
