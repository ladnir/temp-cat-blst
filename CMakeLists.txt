# CMakeLists.txt
cmake_minimum_required(VERSION 3.19)
project(blst C CXX)

option(BLST_BUILD_SHARED "Build shared library" OFF)
option(BLST_PORTABLE     "Force portable C implementation (no assembly)" ON)
option(BLST_ASM          "Enable assembly optimizations" ON)
option(BLST_ADX          "Enable ADX in the c code" ON)
option(BLST_DEFAULT_ASM  "Enable assembly optimizations by default" ON)
option(BLST_ASAN        "Enable Address Sanitizer" OFF)

if( NOT DEFINED BLST_LIMB_SIZE)
	set(BLST_LIMB_SIZE 64)
endif()


option(BLST_CUDA "Enable CUDA support" OFF)

message("BLST_BUILD_SHARED = ${BLST_BUILD_SHARED}")
message("BLST_PORTABLE     = ${BLST_PORTABLE}")
message("BLST_ASM          = ${BLST_ASM}")
message("BLST_ADX          = ${BLST_ADX}")
message("BLST_DEFAULT_ASM  = ${BLST_DEFAULT_ASM}")
message("BLST_CUDA         = ${BLST_CUDA}")
message("BLST_LIMB_SIZE    = ${BLST_LIMB_SIZE}")


# Make sure these vars are copied into the try-compile's cache.
set(CMAKE_TRY_COMPILE_PLATFORM_VARIABLES
    CMAKE_CUDA_FLAGS
    CMAKE_CUDA_HOST_COMPILER
)

# Force the flag into cache so try-compile sees it.
set(CMAKE_CUDA_FLAGS "-allow-unsupported-compiler" CACHE STRING "nvcc flags" FORCE)


# applies to all targets
add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:--allow-unsupported-compiler>)

#add_subdirectory(check)
add_subdirectory(src)
add_subdirectory(test)