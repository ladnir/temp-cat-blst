# CMakeLists.txt
cmake_minimum_required(VERSION 3.27)
project(blst CXX)

option(BLST_BUILD_SHARED "Build shared library" OFF)
option(BLST_PORTABLE     "Build portable implementations and CPU-dispatched assembly fallbacks" ON)
option(BLST_ASM          "Enable assembly optimizations" ON)
option(BLST_ADX          "Enable ADX in the c code" ON)
option(BLST_DEFAULT_ASM  "Enable assembly optimizations by default" ON)
option(BLST_ASAN         "Enable Address Sanitizer" OFF)
option(BLST_AVX512       "Enable AVX512 code paths" OFF)
option(BLST_DEFAULT_AVX  "Prefer AVX implementations by default when available" ${BLST_AVX512})
option(BLST_AVX512_PORTABLE "Enable developer/testing AVX512 emulation" OFF)
option(BLST_ALLOW_UNSUPPORTED_AVX_COMPILER "Allow AVX builds with unsupported compilers" OFF)
option(BLST_ALLOW_UNTESTED_ARM64 "Allow experimental, unqualified ARM64 builds" OFF)
option(BLST_ALLOW_UNTESTED_CHERI "Allow experimental, unqualified CHERI pure-capability builds" OFF)
if(DEFINED BLAS_UBSAN AND NOT DEFINED BLST_UBSAN)
    set(BLST_UBSAN "${BLAS_UBSAN}")
endif()
option(BLST_UBSAN        "Enable Undefined Behavior Sanitizer" OFF)
option(BLST_TEST         "Build test executables" ON)
option(BLST_BENCHMARK    "Build benchmark executables" ON)
option(BLST_CUDA         "Enable CUDA support" OFF)
option(BLST_CUDA_MSM "Enable CUDA support" ${BLST_CUDA})
option(BLST_PARALLEL_MSM "Enable parallel MSM pipeline sources/tests/benchmarks" OFF)
option(BLST_CRNS_BATCH_VROOM_VENDOR "Build the vendored VROOM comparison backend" OFF)

set(BLST_TARGET_IS_X86_64 OFF)
set(BLST_TARGET_IS_ARM64 OFF)
set(BLST_TARGET_IS_APPLE_ARM64 OFF)
if(APPLE AND DEFINED CMAKE_OSX_ARCHITECTURES AND NOT CMAKE_OSX_ARCHITECTURES STREQUAL "")
    if(CMAKE_OSX_ARCHITECTURES STREQUAL "x86_64")
        set(BLST_TARGET_IS_X86_64 ON)
    elseif(CMAKE_OSX_ARCHITECTURES MATCHES "(^|;)arm64(;|$)")
        set(BLST_TARGET_IS_ARM64 ON)
    endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|X86_64|amd64|AMD64|x64|X64)$")
    set(BLST_TARGET_IS_X86_64 ON)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|AARCH64|arm64|ARM64|arm64ec|ARM64EC)$")
    set(BLST_TARGET_IS_ARM64 ON)
endif()
if(CMAKE_GENERATOR_PLATFORM MATCHES "^(ARM64|arm64|ARM64EC|arm64ec)$" OR
   CMAKE_CXX_COMPILER_ARCHITECTURE_ID MATCHES "^(ARM64|arm64|ARM64EC|arm64ec)$")
    set(BLST_TARGET_IS_X86_64 OFF)
    set(BLST_TARGET_IS_ARM64 ON)
endif()
if(BLST_TARGET_IS_ARM64 AND APPLE AND CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(BLST_TARGET_IS_APPLE_ARM64 ON)
endif()

include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_QUIET ON)
check_cxx_source_compiles([[
#ifndef __CHERI_PURE_CAPABILITY__
# error "not a CHERI pure-capability target"
#endif
int main() { return 0; }
]] BLST_TARGET_IS_CHERI_PURECAP)
unset(CMAKE_REQUIRED_QUIET)

if(BLST_TARGET_IS_CHERI_PURECAP)
    set(BLST_TARGET_IS_ARM64 ON)
    if(NOT (BLST_ALLOW_UNTESTED_ARM64 AND BLST_ALLOW_UNTESTED_CHERI))
        message(FATAL_ERROR
            "CHERI pure-capability builds are retained experimental and are not "
            "release-qualified. To attempt an unsupported build, explicitly set "
            "both BLST_ALLOW_UNTESTED_ARM64=ON and BLST_ALLOW_UNTESTED_CHERI=ON.")
    endif()
    message(WARNING
        "Building the unqualified CHERI pure-capability target. It may work, but "
        "it is not fully supported and must not be represented as production-ready.")
elseif(BLST_TARGET_IS_ARM64 AND NOT BLST_TARGET_IS_APPLE_ARM64)
    if(NOT BLST_ALLOW_UNTESTED_ARM64)
        message(FATAL_ERROR
            "ARM64 builds are retained experimental and are not release-qualified. "
            "To attempt an unsupported build, explicitly set "
            "BLST_ALLOW_UNTESTED_ARM64=ON.")
    endif()
    message(WARNING
        "Building the unqualified ARM64 target. It may work, but it is not fully "
        "supported and must not be represented as production-ready.")
endif()

if(BLST_TARGET_IS_APPLE_ARM64)
    set(BLST_APPLE_ARM64_SUPPORTED ON)
endif()

if(BLST_ADX AND NOT BLST_TARGET_IS_X86_64)
    message(STATUS
        "Disabling BLST_ADX for non-x86_64 target "
        "(CMAKE_SYSTEM_PROCESSOR='${CMAKE_SYSTEM_PROCESSOR}', "
        "CMAKE_OSX_ARCHITECTURES='${CMAKE_OSX_ARCHITECTURES}').")
    set(BLST_ADX OFF)
endif()

set(BLST_LIMB_SIZE "64" CACHE STRING "CatBLST limb size (64-bit only)")
set_property(CACHE BLST_LIMB_SIZE PROPERTY STRINGS "64")
if(NOT BLST_LIMB_SIZE STREQUAL "64")
    message(FATAL_ERROR
        "CatBLST supports only BLST_LIMB_SIZE=64; received "
        "'${BLST_LIMB_SIZE}'.")
endif()

if (BLST_ASM)
else()
    set(BLST_PORTABLE ON)
    set(BLST_DEFAULT_ASM OFF)
endif()

if(BLST_DEFAULT_AVX AND NOT (BLST_AVX512 OR BLST_AVX512_PORTABLE))
    message(FATAL_ERROR "BLST_DEFAULT_AVX requires BLST_AVX512 or BLST_AVX512_PORTABLE")
endif()

if(BLST_AVX512)
    set(BLST_AVX_COMPILER_SUPPORTED OFF)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang"
       AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 21.0.0)
        set(BLST_AVX_COMPILER_SUPPORTED ON)
    endif()
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
       AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0.0)
        set(BLST_AVX_COMPILER_SUPPORTED ON)
    endif()

    if(NOT BLST_AVX_COMPILER_SUPPORTED)
        if(BLST_ALLOW_UNSUPPORTED_AVX_COMPILER)
            message(WARNING
                "AVX support requires Clang 21 or newer, or GCC 15 or newer. "
                "Proceeding with unsupported compiler "
                "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} "
                "because BLST_ALLOW_UNSUPPORTED_AVX_COMPILER=ON.")
        else()
            message(FATAL_ERROR
                "AVX support requires Clang 21 or newer, or GCC 15 or newer. "
                "Detected ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}. "
                "Set BLST_ALLOW_UNSUPPORTED_AVX_COMPILER=ON to override this check.")
        endif()
    endif()
endif()

if(MSVC AND NOT DEFINED CMAKE_MSVC_RUNTIME_LIBRARY)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

if(DEFINED BLAS_UBSAN)
    message(WARNING "BLAS_UBSAN is deprecated; use BLST_UBSAN instead.")
endif()

add_library(blst_sanitizers INTERFACE)

if(BLST_ASAN OR BLST_UBSAN)
    if(MSVC)
        if(BLST_ASAN)
            target_compile_options(blst_sanitizers INTERFACE /fsanitize=address)
            target_link_options(blst_sanitizers INTERFACE /fsanitize=address)
        endif()

        if(BLST_UBSAN)
            if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
                target_compile_options(blst_sanitizers INTERFACE /clang:-fsanitize=undefined)
                target_link_options(blst_sanitizers INTERFACE /clang:-fsanitize=undefined)
            else()
                message(WARNING "UBSan is not supported with MSVC, only Clang on Windows")
            endif()
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        set(BLST_SANITIZER_COMPILE_FLAGS)
        set(BLST_SANITIZER_LINK_FLAGS)
        set(BLST_SANITIZER_RPATHS)
        if(BLST_ASAN)
            list(APPEND BLST_SANITIZER_COMPILE_FLAGS -fsanitize=address)
            list(APPEND BLST_SANITIZER_LINK_FLAGS -fsanitize=address)
        endif()
        if(BLST_UBSAN)
            list(APPEND BLST_SANITIZER_COMPILE_FLAGS -fsanitize=undefined)
            list(APPEND BLST_SANITIZER_LINK_FLAGS -fsanitize=undefined)
        endif()
        list(APPEND BLST_SANITIZER_COMPILE_FLAGS -fno-omit-frame-pointer)

        target_compile_options(blst_sanitizers INTERFACE
            $<$<COMPILE_LANGUAGE:C>:${BLST_SANITIZER_COMPILE_FLAGS}>
            $<$<COMPILE_LANGUAGE:CXX>:${BLST_SANITIZER_COMPILE_FLAGS}>)
        target_link_options(blst_sanitizers INTERFACE ${BLST_SANITIZER_LINK_FLAGS})

        if(UNIX)
            set(BLST_SANITIZER_RUNTIME_DIR "${CMAKE_BINARY_DIR}/sanitizer-runtime")

            if(BLST_ASAN)
                execute_process(
                    COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libasan.so
                    OUTPUT_VARIABLE BLST_ASAN_RUNTIME
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
                if(BLST_ASAN_RUNTIME AND EXISTS "${BLST_ASAN_RUNTIME}")
                    file(MAKE_DIRECTORY "${BLST_SANITIZER_RUNTIME_DIR}")
                    get_filename_component(BLST_ASAN_RUNTIME_REALPATH "${BLST_ASAN_RUNTIME}" REALPATH)
                    get_filename_component(BLST_ASAN_RUNTIME_DIR "${BLST_ASAN_RUNTIME_REALPATH}" DIRECTORY)
                    list(APPEND BLST_SANITIZER_RPATHS "${BLST_ASAN_RUNTIME_DIR}")
                    get_filename_component(BLST_ASAN_RUNTIME_NAME "${BLST_ASAN_RUNTIME_REALPATH}" NAME)
                    if(BLST_ASAN_RUNTIME_NAME MATCHES "^(.*\\.so\\.[0-9]+)\\..*$")
                        set(BLST_ASAN_SONAME "${CMAKE_MATCH_1}")
                        execute_process(
                            COMMAND ${CMAKE_COMMAND} -E create_symlink
                                "${BLST_ASAN_RUNTIME_REALPATH}"
                                "${BLST_SANITIZER_RUNTIME_DIR}/${BLST_ASAN_SONAME}")
                        list(APPEND BLST_SANITIZER_RPATHS "${BLST_SANITIZER_RUNTIME_DIR}")
                    endif()
                endif()
            endif()

            if(BLST_UBSAN)
                execute_process(
                    COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libubsan.so
                    OUTPUT_VARIABLE BLST_UBSAN_RUNTIME
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
                if(BLST_UBSAN_RUNTIME AND EXISTS "${BLST_UBSAN_RUNTIME}")
                    file(MAKE_DIRECTORY "${BLST_SANITIZER_RUNTIME_DIR}")
                    get_filename_component(BLST_UBSAN_RUNTIME_REALPATH "${BLST_UBSAN_RUNTIME}" REALPATH)
                    get_filename_component(BLST_UBSAN_RUNTIME_DIR "${BLST_UBSAN_RUNTIME_REALPATH}" DIRECTORY)
                    list(APPEND BLST_SANITIZER_RPATHS "${BLST_UBSAN_RUNTIME_DIR}")
                    get_filename_component(BLST_UBSAN_RUNTIME_NAME "${BLST_UBSAN_RUNTIME_REALPATH}" NAME)
                    if(BLST_UBSAN_RUNTIME_NAME MATCHES "^(.*\\.so\\.[0-9]+)\\..*$")
                        set(BLST_UBSAN_SONAME "${CMAKE_MATCH_1}")
                        execute_process(
                            COMMAND ${CMAKE_COMMAND} -E create_symlink
                                "${BLST_UBSAN_RUNTIME_REALPATH}"
                                "${BLST_SANITIZER_RUNTIME_DIR}/${BLST_UBSAN_SONAME}")
                        list(APPEND BLST_SANITIZER_RPATHS "${BLST_SANITIZER_RUNTIME_DIR}")
                    endif()
                endif()
            endif()

            list(REMOVE_DUPLICATES BLST_SANITIZER_RPATHS)
            foreach(BLST_SANITIZER_RPATH IN LISTS BLST_SANITIZER_RPATHS)
                target_link_options(blst_sanitizers INTERFACE "-Wl,-rpath,${BLST_SANITIZER_RPATH}")
            endforeach()
        endif()
    else()
        message(WARNING
            "Sanitizers requested, but no sanitizer flags are configured for "
            "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}.")
    endif()
endif()



message("Configuration:")
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    message("\n\nWARNING: Not in Release build mode, may be significantly slower.\n\n")
endif()

message("BLST_BUILD_SHARED    = ${BLST_BUILD_SHARED}")
message("BLST_PORTABLE        = ${BLST_PORTABLE}") # enable the portable impl
message("BLST_ASM             = ${BLST_ASM}")      # enable the asm impl (cpu)
message("BLST_ADX             = ${BLST_ADX}")      # enable the adx instruction set (cpu)
message("BLST_DEFAULT_ASM     = ${BLST_DEFAULT_ASM}") # use asm as the default impl.
message("BLST_AVX512          = ${BLST_AVX512}")  #
message("BLST_DEFAULT_AVX     = ${BLST_DEFAULT_AVX}")  #
message("BLST_AVX512_PORTABLE = ${BLST_AVX512_PORTABLE}")  #
message("BLST_ALLOW_UNSUPPORTED_AVX_COMPILER = ${BLST_ALLOW_UNSUPPORTED_AVX_COMPILER}")  #
message("BLST_ALLOW_UNTESTED_ARM64 = ${BLST_ALLOW_UNTESTED_ARM64}")
message("BLST_ALLOW_UNTESTED_CHERI = ${BLST_ALLOW_UNTESTED_CHERI}")
message("BLST_ASAN            = ${BLST_ASAN}")     # enable address sanitizer
message("BLST_UBSAN           = ${BLST_UBSAN}")    # enable undefined behavior sanitizer
message("BLST_CUDA            = ${BLST_CUDA}")     # enable the cude impl (gpu)
message("BLST_CUDA_MSM        = ${BLST_CUDA_MSM}")     # enable the cude impl (gpu)
message("BLST_PARALLEL_MSM    = ${BLST_PARALLEL_MSM}") # enable parallel MSM pipeline
message("BLST_CRNS_BATCH_VROOM_VENDOR = ${BLST_CRNS_BATCH_VROOM_VENDOR}") # build vendored VROOM comparison backend
message("BLST_FAST_BUILD      = ${BLST_FAST_BUILD}") # make the cude build faster but slower runtime.
message("BLST_LIMB_SIZE       = ${BLST_LIMB_SIZE}")
message("BLST_TEST            = ${BLST_TEST}")      # build tests
message("BLST_BENCHMARK       = ${BLST_BENCHMARK}") # build benchmarks


message(STATUS "C compiler: ${CMAKE_C_COMPILER}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ compiler ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "C++ compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

if(BLST_CRNS_BATCH_VROOM_VENDOR)
    message(WARNING "BLST_CRNS_BATCH_VROOM_VENDOR is temporarily unavailable during the VROOM backend reset.")
endif()
set(BLST_HAVE_CRNS_BATCH_VROOM_VENDOR OFF)


# Make sure these vars are copied into the try-compile's cache.
set(CMAKE_TRY_COMPILE_PLATFORM_VARIABLES
    CMAKE_CUDA_FLAGS
    CMAKE_CUDA_HOST_COMPILER
)

# Force the flag into cache so try-compile sees it.
set(CMAKE_CUDA_FLAGS "-allow-unsupported-compiler" CACHE STRING "nvcc flags" FORCE)

if(BLST_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        message(WARNING "CMAKE_CUDA_ARCHITECTURES not defined, defaulting to 89")
    endif()
        set(CMAKE_CUDA_ARCHITECTURES 89)
endif()

# applies to all targets
add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:--allow-unsupported-compiler>)
if(MSVC)
    add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/Zc:preprocessor>)
endif()

#add_subdirectory(check)
add_subdirectory(blst)
add_subdirectory(BTE)

enable_testing()

if(BLST_TEST)
    add_subdirectory(test)
endif()

# Add benchmark subdirectory if benchmarks are enabled
if(BLST_BENCHMARK)
    add_subdirectory(benchmark)
endif()

