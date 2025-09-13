

set(BLST_SRC server.c)

# ---- Static Library Configuration ----
set(BLST_STATIC_SRC ${BLST_SRC})
add_compile_options(-save-temps=obj)

if(MSVC)
    set(BLST_DEBUG_OPT /Zi /Od /MDd /DEBUG /W4 /WX)
    set(BLST_RELEASE_OPT /O2 /MD)
    set(BLST_RELWITHDEBINFO_OPT /Zi /O2 /MD /DEBUG)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(ASM_PATTERN "${CMAKE_CURRENT_SOURCE_DIR}/../build/win64/*-x86_64.asm")
    else()
        set(ASM_PATTERN "${CMAKE_CURRENT_SOURCE_DIR}/../build/win64/*-armv8.asm")
    endif()
else()
    set(BLST_DEBUG_OPT -g -O0 -fPIC)
    set(BLST_RELEASE_OPT -O2 -Wall -Wextra -fPIC)
    set(BLST_RELWITHDEBINFO_OPT -g -O2 -Wall -Wextra -fPIC)

    set(ASM_PATTERN  ${CMAKE_CURRENT_SOURCE_DIR}/../build/assembly.S)
endif()


if(BLST_ASM)
    file(GLOB ASM_SRC ${ASM_PATTERN})
    #list(APPEND BLST_STATIC_SRC ${ASM_SRC})
    message(STATUS "Assembly optimizations enabled, searching for ${ASM_PATTERN}")
    message("found: ${ASM_SRC}")
    
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
    
    # Tell CMake to compile .c files with NVCC
    set_source_files_properties(${BLST_SRC} PROPERTIES LANGUAGE CUDA)
    
    set(BLST_STATIC_OPTS 
        $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CUDA>>:
          $<LIST:TRANSFORM,${BLST_DEBUG_OPT},PREPEND,-Xcompiler=>>
        $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CUDA>>:
          $<LIST:TRANSFORM,${BLST_RELEASE_OPT},PREPEND,-Xcompiler=>>
        $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANGUAGE:CUDA>>:
          $<LIST:TRANSFORM,${BLST_RELWITHDEBINFO_OPT},PREPEND,-Xcompiler=>>

        --expt-relaxed-constexpr
        --extended-lambda
        -arch=sm_70
        -rdc=true
        # Add verbose compilation flags
        #--verbose
        #--ptxas-options=-v
    )

else()

    
    #set(BLST_STATIC_OPTS 
    #    $<$<CONFIG:Debug>:${BLST_DEBUG_OPT}>
    #    $<$<CONFIG:Release>:${BLST_RELEASE_OPT}>
    #    $<$<CONFIG:RelWithDebInfo>:${BLST_RELWITHDEBINFO_OPT}>
    #)
    # Apply MSVC flags only to C files
    set(BLST_STATIC_OPTS 
        $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:${BLST_DEBUG_OPT}>
        $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:${BLST_RELEASE_OPT}>
        $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANGUAGE:C>>:${BLST_RELWITHDEBINFO_OPT}>
    )
endif()

add_library(blst_static ${BLST_STATIC_SRC} ${ASM_SRC})
set_target_properties(blst_static PROPERTIES OUTPUT_NAME blst)
target_compile_options(blst_static PRIVATE ${BLST_STATIC_OPTS})




configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

target_include_directories(blst_static PUBLIC ${CMAKE_CURRENT_BINARY_DIR})



if(BLST_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    target_link_libraries(blst_static PRIVATE CUDA::cudart)
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


if(BLST_ADX)
    if (MSVC)
      target_compile_options(blst_static PUBLIC /arch:AVX2)
    else()
      target_compile_options(blst_static PUBLIC -mbmi2 -madx)
    endif()
endif()


#arch native
if (MSVC)
else()
    target_compile_options(blst_static PUBLIC -march=native)
endif()

