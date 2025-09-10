
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
"Fp_tests.cpp" 
"Fp2_tests.cpp" 
"Fp12_tests.cpp" 
"Fr_tests.cpp" 
"BindingCpp_tests.cpp"
"G1_tests.cpp"
"G2_tests.cpp"
"GT_tests.cpp"
"Pairing_tests.cpp"
)
target_include_directories(blst_tests PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(blst_tests PRIVATE gtest_main blst_static)

target_compile_definitions(blst_tests PRIVATE ${BLST_STATIC_DEFS})
target_compile_options(blst_tests PRIVATE ${BLST_STATIC_OPTS})
add_test(NAME blst_tests COMMAND blst_tests)