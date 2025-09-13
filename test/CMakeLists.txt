
# ---- Testing ----
include(FetchContent)
#cmake_policy(SET CMP0135 NEW)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()
add_executable(blst_tests 
"Fp_tests.cpp" 
"Fp2_tests.cpp" 
"Fp12_tests.cpp" 
"Fr_tests.cpp" 
"BindingCpp_tests.cpp"
"G1_tests.cpp"
"G2_tests.cpp"
"GT_tests.cpp"
"Pairing_tests.cpp"
 "Interface_tests.cpp")
target_include_directories(blst_tests PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(blst_tests PRIVATE gtest_main blst_static)

target_compile_definitions(blst_tests PRIVATE ${BLST_STATIC_DEFS})
target_compile_options(blst_tests PRIVATE ${BLST_STATIC_OPTS})
add_test(NAME blst_tests COMMAND blst_tests)



if(MSVC)
    set(BLST_DEBUG_OPT /Zi /Od /MDd /DEBUG /W4 /WX)
    set(BLST_RELEASE_OPT /O2 /MD)
    set(BLST_RELWITHDEBINFO_OPT /Zi /O2 /MD /DEBUG)
else()
    set(BLST_DEBUG_OPT -g -O0 -fPIC)
    set(BLST_RELEASE_OPT -O2 -Wall -Wextra -fPIC)
    set(BLST_RELWITHDEBINFO_OPT -g -O2 -Wall -Wextra -fPIC)
endif()

set(BLST_STATIC_OPTS 
    $<$<CONFIG:Debug>:${BLST_DEBUG_OPT}>
    $<$<CONFIG:Release>:${BLST_RELEASE_OPT}>
    $<$<CONFIG:RelWithDebInfo>:${BLST_RELWITHDEBINFO_OPT}>
)
target_compile_options(blst_tests PRIVATE ${BLST_STATIC_OPTS})
