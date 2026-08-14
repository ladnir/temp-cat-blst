


# ---- Static Library Configuration ----
set(BLST_STATIC_SRC  
    "backend/backend_portable.cpp"
    "common/rb_tree.cpp"
    "sign/keygen.cpp"
    "scalar/fields/Fr.cpp"
    "scalar/fields/Fp.cpp"
    "scalar/fields/Fp12.cpp"
    "scalar/curve/GT.cpp"
    "scalar/curve/point.cpp"
    "scalar/curve/multi_scalar.cpp"
    "scalar/curve/multi_scalar_batch_affine.cpp"
    "sign/aggregate.cpp"
    "c_api.cpp"
    "avx_c_api.cpp"
    "fft/fr_interpolation.cpp"
    "fft_c_api.cpp"
    "middle_product_c_api.cpp"
    "exports.cpp"
)

if(BLST_PARALLEL_MSM)
    list(APPEND BLST_STATIC_SRC
        "cuda/curve/multi_scalar_par_G1.cpp"
        "cuda/curve/multi_scalar_par_G2.cpp"
        "cuda/curve/multi_scalar_par_GT.cpp"
    )
endif()

if(BLST_AVX512_PORTABLE)
    list(APPEND BLST_STATIC_SRC
        "backend/avx512/intrin_compat.cpp"
    )
endif()

if(MSVC)
    set(BLST_OPTS)
    set(BLST_DEBUG_OPT)
    set(BLST_RELEASE_OPT)
    set(BLST_RELWITHDEBINFO_OPT)
    if(BLST_TARGET_IS_X86_64)
        set(ASM_PATTERN "${CMAKE_CURRENT_SOURCE_DIR}/../build/win64/*-x86_64.asm")
    elseif(BLST_TARGET_IS_ARM64)
        set(ASM_PATTERN "${CMAKE_CURRENT_SOURCE_DIR}/../build/win64/*-armv8.asm")
    else()
        message(FATAL_ERROR "Unsupported MSVC assembly target")
    endif()
else()
    set(BLST_OPTS)
    set(BLST_DEBUG_OPT)
    set(BLST_RELEASE_OPT -fwrapv -Wall -Wextra)
    set(BLST_RELWITHDEBINFO_OPT -fwrapv -Wall -Wextra)

    set(ASM_PATTERN  ${CMAKE_CURRENT_SOURCE_DIR}/../build/assembly.S)
endif()


if(BLST_ADX)
    if (MSVC)
      set(BLST_OPTS  ${BLST_OPTS} /arch:AVX2)

      if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        set(BLST_OPTS  ${BLST_OPTS} /clang:-mbmi /clang:-mbmi2 /clang:-madx)
      endif()
    else()
      set(BLST_OPTS  ${BLST_OPTS} -mbmi2 -madx)
    endif()
endif()


if(BLST_ASM)
    list(APPEND BLST_STATIC_SRC "common/cpuid.cpp")
    file(GLOB ASM_SRC ${ASM_PATTERN})
    #list(APPEND BLST_STATIC_SRC ${ASM_SRC})
    message(STATUS "Assembly optimizations enabled, searching for ${ASM_PATTERN}")
    #message("found: ${ASM_SRC}")
    
    if(MSVC)
        enable_language(ASM_MASM)
        set_source_files_properties(${ASM_SRC} PROPERTIES LANGUAGE ASM_MASM)
        if(BLST_PORTABLE)
            # Match build/assembly.S: BLST_PORTABLE is the public build option,
            # while __BLST_PORTABLE__ is the assembly-internal dispatch guard.
            set_source_files_properties(${ASM_SRC} PROPERTIES
                COMPILE_DEFINITIONS "__BLST_PORTABLE__")
        endif()
    else()
        enable_language(ASM)
        set_source_files_properties(${ASM_SRC} PROPERTIES LANGUAGE ASM)
    endif()

endif()

if(BLST_CUDA)

    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    
    #set(BLST_STATIC_SRC ${BLST_STATIC_SRC} cuda_kernel.cpp)
    set(BLST_STATIC_SRC ${BLST_STATIC_SRC} 
        cuda/fields/Fr.cpp
        cuda/fields/Fp.cpp
        cuda/fields/Fp2.cpp
        cuda/fields/Fp12.cpp
        cuda/curve/G1.cpp
        cuda/curve/G2.cpp
        cuda/curve/GT.cpp
    )
    # Tell CMake to compile .c files with NVCC
    set_source_files_properties(${BLST_STATIC_SRC} PROPERTIES LANGUAGE CUDA)
    set(BLST_OPTS
        # CPU flags forwarded to NVCC's host compiler (only when list non-empty)
        $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<BOOL:${BLST_OPTS}>>:
        $<LIST:TRANSFORM,${BLST_OPTS},PREPEND,-Xcompiler=>>

        # Config-specific CPU flags forwarded to NVCC's host compiler
        $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CUDA>>:
        $<LIST:TRANSFORM,${BLST_DEBUG_OPT},PREPEND,-Xcompiler=>>
        $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CUDA>>:
        $<LIST:TRANSFORM,${BLST_RELEASE_OPT},PREPEND,-Xcompiler=>>
        $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANGUAGE:CUDA>>:
        $<LIST:TRANSFORM,${BLST_RELWITHDEBINFO_OPT},PREPEND,-Xcompiler=>>

        # CUDA device/linkage flags (CUDA-only)
        $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
        $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>
        #$<$<COMPILE_LANGUAGE:CUDA>:-arch=sm_89>
        $<$<COMPILE_LANGUAGE:CUDA>:-rdc=true>
        #$<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>           # cheap line tables
        $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CUDA>:-G -g -Xptxas -O0>> # expensive debug info
        $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CUDA>:-Xptxas -O3>> # 
        
        #verbose output for nvcc
        #$<$<COMPILE_LANGUAGE:CUDA>:
        #    -v                       # print every sub-command
        #    --keep --keep-dir=${CMAKE_BINARY_DIR}/nvcc_keep
        #    -Xptxas -v               # ptxas resource report (regs, smem, spills)
        #    -MMD
        #>
    )

    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    
        list(APPEND BLST_OPTS $<$<COMPILE_LANGUAGE:CUDA>:-lineinfo -Xptxas -O0>)

    endif()

    
else()

    
    set(BLST_OPTS  
        $<$<COMPILE_LANGUAGE:CXX>:${BLST_OPTS}>
        $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:${BLST_DEBUG_OPT}>
        $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:${BLST_RELEASE_OPT}>
        $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANGUAGE:CXX>>:${BLST_RELWITHDEBINFO_OPT}>
    )
endif()

add_library(blst_static ${BLST_STATIC_SRC} ${ASM_SRC})

if(BLST_FUZZ)
    target_compile_options(blst_static PRIVATE -fsanitize=fuzzer-no-link)
endif()
add_library(category_blst::blst ALIAS blst_static)
set_target_properties(blst_static PROPERTIES
    OUTPUT_NAME blst
    POSITION_INDEPENDENT_CODE ON)
#message(STATUS "BLST_OPTS: ${BLST_OPTS}")
target_compile_options(blst_static PRIVATE ${BLST_OPTS})
target_include_directories(blst_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
if(BLST_CUDA OR BLST_PARALLEL_MSM)
    target_compile_features(blst_static PUBLIC cxx_std_23)
else()
    target_compile_features(blst_static PUBLIC cxx_std_20)
endif()
if(BLST_CUDA)
    set_target_properties(blst_static PROPERTIES
        CUDA_STANDARD 23
        CUDA_STANDARD_REQUIRED ON
        CUDA_EXTENSIONS OFF)
endif()
# for in source builds, we need to include the build dir too (config.h)
target_include_directories(blst_static PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/..)


configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

target_include_directories(blst_static PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(avx)
if(BLST_DEFAULT_AVX)
    # The selected native C++ facade exposes AVX types and inline operations,
    # so its exact ISA compile contract is a legitimate consumer requirement.
    target_link_libraries(blst_static PUBLIC blst_avx_headers)
else()
    target_link_libraries(blst_static PRIVATE blst_avx_headers)
endif()

if(BLST_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    target_link_libraries(blst_static PRIVATE CUDA::cudart)
    target_include_directories(blst_static PUBLIC $<$<CXX_COMPILER_ID:MSVC>:${CUDAToolkit_INCLUDE_DIRS}>)
endif()

if(TARGET blst_sanitizers)
    target_link_libraries(blst_static PUBLIC blst_sanitizers)
endif()
