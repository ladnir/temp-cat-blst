# CMakeLists.txt
cmake_minimum_required(VERSION 3.19)
project(blst C CXX)

option(BLST_BUILD_SHARED "Build shared library" OFF)
option(BLST_PORTABLE "Force portable C implementation (no assembly)" OFF)
option(BLST_CUDA "Enable CUDA support" OFF)

message("BLST_BUILD_SHARED = ${BLST_BUILD_SHARED}")
message("BLST_PORTABLE     = ${BLST_PORTABLE}")
message("BLST_CUDA         = ${BLST_CUDA}")


add_subdirectory(src)

add_subdirectory(test)