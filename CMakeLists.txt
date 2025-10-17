# CMakeLists.txt
cmake_minimum_required(VERSION 3.27)
project(blst CXX)

option(BLST_BUILD_SHARED "Build shared library" OFF)
option(BLST_PORTABLE     "Force portable C implementation (no assembly)" ON)
option(BLST_ASM          "Enable assembly optimizations" ON)
option(BLST_ADX          "Enable ADX in the c code" ON)
option(BLST_DEFAULT_ASM  "Enable assembly optimizations by default" ON)
option(BLST_ASAN         "Enable Address Sanitizer" OFF)
option(BLAS_UBSAN        "Enable UB Sanitizer" OFF)
option(BLST_BENCHMARK     "Build benchmark executables" ON)

if( NOT DEFINED BLST_LIMB_SIZE)
	set(BLST_LIMB_SIZE 64)
endif()

if (BLST_ASM)
else()
    set(BLST_PORTABLE ON)
    set(BLST_DEFAULT_ASM OFF)
endif()

set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")



option(BLST_CUDA "Enable CUDA support" OFF)

message("BLST_BUILD_SHARED = ${BLST_BUILD_SHARED}")
message("BLST_PORTABLE     = ${BLST_PORTABLE}") # enable the portable impl
message("BLST_ASM          = ${BLST_ASM}")      # enable the asm impl (cpu)
message("BLST_ADX          = ${BLST_ADX}")      # enable the adx instruction set (cpu)
message("BLST_DEFAULT_ASM  = ${BLST_DEFAULT_ASM}") # use asm as the default impl.
message("BLST_CUDA         = ${BLST_CUDA}")     # enable the cude impl (gpu)
message("BLST_FAST_BUILD   = ${BLST_FAST_BUILD}") # make the cude build faster but slower runtime.
message("BLST_LIMB_SIZE    = ${BLST_LIMB_SIZE}") # 32 or 64.
message("BLST_BENCHMARK    = ${BLST_BENCHMARK}") # build benchmarks


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

add_subdirectory(test)

# Add benchmark subdirectory if benchmarks are enabled
if(BLST_BENCHMARK)
    add_subdirectory(benchmark)
endif()

