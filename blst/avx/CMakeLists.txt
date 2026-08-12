# The selected AVX backend is header-driven. This target carries its include
# and compile usage requirements without external runtime dependencies.

add_library(blst_avx_headers INTERFACE)
target_compile_features(blst_avx_headers INTERFACE cxx_std_20)

if(BLST_AVX512)
    if(MSVC)
        target_compile_options(blst_avx_headers INTERFACE
            $<$<COMPILE_LANGUAGE:CXX>:/clang:-mavx2>
            $<$<COMPILE_LANGUAGE:CXX>:/clang:-mavx512f>
            $<$<COMPILE_LANGUAGE:CXX>:/clang:-mavx512dq>
            $<$<COMPILE_LANGUAGE:CXX>:/clang:-mavx512ifma>)
    else()
        target_compile_options(blst_avx_headers INTERFACE
            $<$<COMPILE_LANGUAGE:CXX>:-mavx2>
            $<$<COMPILE_LANGUAGE:CXX>:-mavx512f>
            $<$<COMPILE_LANGUAGE:CXX>:-mavx512dq>
            $<$<COMPILE_LANGUAGE:CXX>:-mavx512ifma>)
    endif()
endif()
