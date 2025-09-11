# CMakeLists.txt
cmake_minimum_required(VERSION 3.19)
project(blst C CXX)

option(BLST_BUILD_SHARED "Build shared library" OFF)
option(BLST_PORTABLE     "Force portable C implementation (no assembly)" ON)
option(BLST_ASM          "Enable assembly optimizations" ON)
option(BLST_ADX          "Enable ADX in the c code" ON)
option(BLST_DEFAULT_ASM  "Enable assembly optimizations by default" ON)

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


add_subdirectory(src)

add_subdirectory(test)