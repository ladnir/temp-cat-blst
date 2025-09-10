# CMakeLists.txt
cmake_minimum_required(VERSION 3.19)
project(blst C CXX)

option(BLST_BUILD_SHARED "Build shared library" OFF)
option(BLST_PORTABLE "Force portable C implementation (no assembly)" OFF)

message("BLST_BUILD_SHARED = ${BLST_BUILD_SHARED}")
message("BLST_PORTABLE     = ${BLST_PORTABLE}")

set(BLST_SRC src/server.c)

# ---- Static Library Configuration ----
set(BLST_STATIC_SRC ${BLST_SRC})
set(BLST_STATIC_OPTS "")

if(MSVC)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(ASM_PATTERN "${CMAKE_CURRENT_SOURCE_DIR}/build/win64/*-x86_64.asm")
    else()
        set(ASM_PATTERN "${CMAKE_CURRENT_SOURCE_DIR}/build/win64/*-armv8.asm")
    endif()
    file(GLOB ASM_SRC ${ASM_PATTERN})
    if(ASM_SRC AND NOT BLST_PORTABLE)
        list(APPEND BLST_STATIC_SRC ${ASM_SRC})
    endif()
    
    set(BLST_STATIC_OPTS 
        $<$<CONFIG:Debug>:/Zi /Od /MDd /DEBUG>
        $<$<CONFIG:Release>:/O2 /MD>
        $<$<CONFIG:RelWithDebInfo>:/Zi /O2 /MD /DEBUG>
    )
else()
    set(ASM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/build/assembly.S)
    if(EXISTS ${ASM_SRC} AND NOT BLST_PORTABLE)
        list(APPEND BLST_STATIC_SRC ${ASM_SRC})
    endif()
    
    set(BLST_STATIC_OPTS 
        $<$<CONFIG:Debug>:-g -O0 -fPIC>
        $<$<CONFIG:Release>:-O2 -Wall -Wextra -fPIC>
        $<$<CONFIG:RelWithDebInfo>:-g -O2 -Wall -Wextra -fPIC>
    )
endif()

add_library(blst_static ${BLST_STATIC_SRC})
set_target_properties(blst_static PROPERTIES OUTPUT_NAME blst)
target_compile_options(blst_static PRIVATE ${BLST_STATIC_OPTS})

# ---- Shared Library Configuration ----
if(BLST_BUILD_SHARED)
    set(BLST_SHARED_SRC ${BLST_SRC})
    set(BLST_SHARED_DEFS "")
    set(BLST_SHARED_OPTS "")
    set(BLST_SHARED_PROPS OUTPUT_NAME blst)

    if(MSVC)
        if(ASM_SRC AND NOT BLST_PORTABLE)
            list(APPEND BLST_SHARED_SRC ${ASM_SRC})
        endif()
        list(APPEND BLST_SHARED_DEFS __BLST_DLL_MAIN__)
        set(BLST_SHARED_OPTS /MD)
        set(BLST_SHARED_PROPS ${BLST_SHARED_PROPS} LINK_FLAGS "/DEF:${CMAKE_CURRENT_SOURCE_DIR}/build/win64/blst.def")
    else()
        if(EXISTS ${ASM_SRC} AND NOT BLST_PORTABLE)
            list(APPEND BLST_SHARED_SRC ${ASM_SRC})
        endif()
        set(BLST_SHARED_OPTS -fPIC)
        if(APPLE)
            set(BLST_SHARED_PROPS ${BLST_SHARED_PROPS} SUFFIX ".dylib")
        endif()
    endif()

    add_library(blst_shared SHARED ${BLST_SHARED_SRC})
    set_target_properties(blst_shared PROPERTIES ${BLST_SHARED_PROPS})
    target_compile_options(blst_shared PRIVATE ${BLST_SHARED_OPTS})
endif()

# ---- Configure config.h ----
if(BLST_PORTABLE)
    set(__BLST_PORTABLE__ 1)
    set(__BLST_NO_ASM__ 1)
else()
    unset(__BLST_PORTABLE__)
    unset(__BLST_NO_ASM__)
endif()

if(MSVC AND BLST_BUILD_SHARED)
    set(__BLST_DLL_MAIN__ 1)
else()
    unset(__BLST_DLL_MAIN__)
endif()

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

# ---- Testing ----
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()
add_executable(blst_tests 
"test/Fp_tests.cpp" 
"test/Fp2_tests.cpp" 
"test/Fp12_tests.cpp" 
"test/Fr_tests.cpp" 
"test/BindingCpp_tests.cpp"
"test/G1_tests.cpp"
"test/G2_tests.cpp"
"test/GT_tests.cpp"
"test/Pairing_tests.cpp"
)
target_include_directories(blst_tests PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(blst_tests PRIVATE gtest_main blst_static)

target_compile_definitions(blst_tests PRIVATE ${BLST_STATIC_DEFS})
target_compile_options(blst_tests PRIVATE ${BLST_STATIC_OPTS})
add_test(NAME blst_tests COMMAND blst_tests)