# CMakeLists.txt
cmake_minimum_required(VERSION 3.27)
project(blst CXX)

option(BLST_BUILD_SHARED "Build shared library" OFF)
option(BLST_PORTABLE     "Force portable C implementation (no assembly)" ON)
option(BLST_ASM          "Enable assembly optimizations" ON)
option(BLST_ADX          "Enable ADX in the c code" ON)
option(BLST_DEFAULT_ASM  "Enable assembly optimizations by default" ON)
option(BLST_ASAN         "Enable Address Sanitizer" OFF)
option(BLST_AVX512       "Enable AVX512 code paths" OFF)
option(BLST_DEFAULT_AVX  "Prefer AVX implementations by default when available" ${BLST_AVX512})
option(BLST_AVX512_PORTABLE "Enable AVX512 emulation code paths" OFF)
option(BLAS_UBSAN        "Enable UB Sanitizer" OFF)
option(BLST_TEST         "Build test executables" ON)
option(BLST_BENCHMARK    "Build benchmark executables" ON)
option(BLST_CUDA         "Enable CUDA support" OFF)
option(BLST_CUDA_MSM "Enable CUDA support" ${BLST_CUDA})
option(BLST_PARALLEL_MSM "Enable parallel MSM pipeline sources/tests/benchmarks" OFF)
option(BLST_CRNS_BATCH_VROOM_VENDOR "Build the vendored VROOM comparison backend" OFF)

if( NOT DEFINED BLST_LIMB_SIZE)
	set(BLST_LIMB_SIZE 64)
endif()

if (BLST_ASM)
else()
    set(BLST_PORTABLE ON)
    set(BLST_DEFAULT_ASM OFF)
endif()

if(BLST_DEFAULT_AVX AND NOT (BLST_AVX512 OR BLST_AVX512_PORTABLE))
    message(FATAL_ERROR "BLST_DEFAULT_AVX requires BLST_AVX512 or BLST_AVX512_PORTABLE")
endif()

set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")



message("Configuration:")
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
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
message("BLST_CUDA            = ${BLST_CUDA}")     # enable the cude impl (gpu)
message("BLST_CUDA_MSM        = ${BLST_CUDA_MSM}")     # enable the cude impl (gpu)
message("BLST_PARALLEL_MSM    = ${BLST_PARALLEL_MSM}") # enable parallel MSM pipeline
message("BLST_CRNS_BATCH_VROOM_VENDOR = ${BLST_CRNS_BATCH_VROOM_VENDOR}") # build vendored VROOM comparison backend
message("BLST_FAST_BUILD      = ${BLST_FAST_BUILD}") # make the cude build faster but slower runtime.
message("BLST_LIMB_SIZE       = ${BLST_LIMB_SIZE}") # 32 or 64.
message("BLST_TEST            = ${BLST_TEST}")      # build tests
message("BLST_BENCHMARK       = ${BLST_BENCHMARK}") # build benchmarks


message(STATUS "C compiler: ${CMAKE_C_COMPILER}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ compiler ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "C++ compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

if(BLST_CRNS_BATCH_VROOM_VENDOR)
    find_path(BLST_GMP_INCLUDE_DIR gmpxx.h)
    find_library(BLST_GMP_LIBRARY gmp)
    find_library(BLST_GMPXX_LIBRARY gmpxx)

    if(BLST_GMP_INCLUDE_DIR AND BLST_GMP_LIBRARY AND BLST_GMPXX_LIBRARY)
        set(BLST_HAVE_CRNS_BATCH_VROOM_VENDOR ON)
        message(STATUS "Vendored VROOM backend enabled")
        message(STATUS "  GMP include: ${BLST_GMP_INCLUDE_DIR}")
        message(STATUS "  GMP libs: ${BLST_GMPXX_LIBRARY};${BLST_GMP_LIBRARY}")
    else()
        set(BLST_HAVE_CRNS_BATCH_VROOM_VENDOR OFF)
        message(WARNING "BLST_CRNS_BATCH_VROOM_VENDOR requested, but GMP/gmpxx was not found. Vendored VROOM targets will be skipped.")
    endif()
else()
    set(BLST_HAVE_CRNS_BATCH_VROOM_VENDOR OFF)
endif()


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

#add_subdirectory(check)
add_subdirectory(blst)

enable_testing()

if(BLST_TEST)
    add_subdirectory(test)
endif()

# Add benchmark subdirectory if benchmarks are enabled
if(BLST_BENCHMARK)
    add_subdirectory(benchmark)
endif()

