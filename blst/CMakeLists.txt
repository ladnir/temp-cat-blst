


# ---- Static Library Configuration ----
set(BLST_STATIC_SRC  mains/server.cpp)

if(MSVC)
    set(BLST_OPTS )
    set(BLST_DEBUG_OPT /Zi /Od /DEBUG ) #/W4 /WX
    set(BLST_RELEASE_OPT /O2 )
    set(BLST_RELWITHDEBINFO_OPT /Zi /O2 /DEBUG)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(ASM_PATTERN "${CMAKE_CURRENT_SOURCE_DIR}/../build/win64/*-x86_64.asm")
    else()
        set(ASM_PATTERN "${CMAKE_CURRENT_SOURCE_DIR}/../build/win64/*-armv8.asm")
    endif()
else()
    set(BLST_OPTS -fPIC)
    set(BLST_DEBUG_OPT -g -O0)
    set(BLST_RELEASE_OPT -O2 -Wall -Wextra)
    set(BLST_RELWITHDEBINFO_OPT -g -O2 -Wall -Wextra)

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
    file(GLOB ASM_SRC ${ASM_PATTERN})
    #list(APPEND BLST_STATIC_SRC ${ASM_SRC})
    message(STATUS "Assembly optimizations enabled, searching for ${ASM_PATTERN}")
    #message("found: ${ASM_SRC}")
    
    if(MSVC)
        enable_language(ASM_MASM)
        set_source_files_properties(${ASM_SRC} PROPERTIES LANGUAGE ASM_MASM)
    else()
        enable_language(ASM)
        set_source_files_properties(${ASM_SRC} PROPERTIES LANGUAGE ASM)
    endif()

endif()

if(BLST_CUDA)

    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    
    set(BLST_STATIC_SRC ${BLST_STATIC_SRC} cuda_kernel.cpp)

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
        #$<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>
        #$<$<COMPILE_LANGUAGE:CUDA>:-arch=sm_89>
        $<$<COMPILE_LANGUAGE:CUDA>:-rdc=true>
        #$<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>           # cheap line tables
        #$<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CUDA>:-G -g -Xptxas -O0>> # expensive debug info
        
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
set_target_properties(blst_static PROPERTIES OUTPUT_NAME blst)
#message(STATUS "BLST_OPTS: ${BLST_OPTS}")
target_compile_options(blst_static PUBLIC ${BLST_OPTS})
target_include_directories(blst_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

# for in source builds, we need to include the build dir too (config.h)
target_include_directories(blst_static PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/..)


configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

target_include_directories(blst_static PUBLIC ${CMAKE_CURRENT_BINARY_DIR})



if(BLST_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    target_link_libraries(blst_static PRIVATE CUDA::cudart)
    target_include_directories(blst_static PUBLIC $<$<CXX_COMPILER_ID:MSVC>:${CUDAToolkit_INCLUDE_DIRS}>)
endif()



# Function to enable ASan for a target
if(BLST_ASAN)
    if(MSVC)
        # Enable ASan for Visual Studio
        target_compile_options(blst_static PRIVATE /fsanitize=address)
            
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        # For GCC and Clang
        target_compile_options(blst_static PRIVATE -fsanitize=address -fno-omit-frame-pointer)
        target_link_options(blst_static PRIVATE -fsanitize=address)
    endif()
endif()

if(BLAS_UBSAN)
    if(MSVC)
        # Enable UBSan for Visual Studio
        if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        # UBSAN requires no iterator debugging in MSVC
            #target_compile_definitions(blst_static PUBLIC _ITERATOR_DEBUG_LEVEL=0;_HAS_ITERATOR_DEBUGGING=0)
            target_compile_options(blst_static PUBLIC /clang:-fsanitize=undefined)
        else()
            message(WARNING "UBSan is not supported with MSVC, only Clang on Windows")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        # For GCC and Clang
        target_compile_options(blst_static PUBLIC -fsanitize=undefined -fno-omit-frame-pointer)
        target_link_options(blst_static PUBLIC -fsanitize=undefined)
    endif()
endif()




#arch native
if (MSVC)
else()
    target_compile_options(blst_static PUBLIC -march=native)
endif()

