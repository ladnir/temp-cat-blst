

set(BLST_SRC server.c)

# ---- Static Library Configuration ----
set(BLST_STATIC_SRC ${BLST_SRC})


if(MSVC)
    set(BLST_DEBUG_OPT /Zi /Od /MDd /DEBUG)
    set(BLST_RELEASE_OPT /O2 /MD)
    set(BLST_RELWITHDEBINFO_OPT /Zi /O2 /MD /DEBUG)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(ASM_PATTERN "${CMAKE_CURRENT_SOURCE_DIR}/build/win64/*-x86_64.asm")
    else()
        set(ASM_PATTERN "${CMAKE_CURRENT_SOURCE_DIR}/build/win64/*-armv8.asm")
    endif()
else()
    set(BLST_DEBUG_OPT -g -O0 -fPIC)
    set(BLST_RELEASE_OPT -O2 -Wall -Wextra -fPIC)
    set(BLST_RELWITHDEBINFO_OPT -g -O2 -Wall -Wextra -fPIC)

    set(ASM_PATTERN  ${CMAKE_CURRENT_SOURCE_DIR}/build/assembly.S)
endif()


file(GLOB ASM_SRC ${ASM_PATTERN})
if(ASM_SRC AND NOT BLST_PORTABLE)
    list(APPEND BLST_STATIC_SRC ${ASM_SRC})
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

    
    set(BLST_STATIC_OPTS 
        $<$<CONFIG:Debug>:${BLST_DEBUG_OPT}>
        $<$<CONFIG:Release>:${BLST_RELEASE_OPT}>
        $<$<CONFIG:RelWithDebInfo>:${BLST_RELWITHDEBINFO_OPT}>
    )
endif()

add_library(blst_static ${BLST_STATIC_SRC})
set_target_properties(blst_static PROPERTIES OUTPUT_NAME blst)
target_compile_options(blst_static PRIVATE ${BLST_STATIC_OPTS})

# ---- Shared Library Configuration ----
#if(BLST_BUILD_SHARED)
#    set(BLST_SHARED_SRC ${BLST_SRC})
#    set(BLST_SHARED_DEFS "")
#    set(BLST_SHARED_OPTS "")
#    set(BLST_SHARED_PROPS OUTPUT_NAME blst)
#
#    if(MSVC)
#        if(ASM_SRC AND NOT BLST_PORTABLE)
#            list(APPEND BLST_SHARED_SRC ${ASM_SRC})
#        endif()
#        list(APPEND BLST_SHARED_DEFS __BLST_DLL_MAIN__)
#        set(BLST_SHARED_OPTS /MD)
#        set(BLST_SHARED_PROPS ${BLST_SHARED_PROPS} LINK_FLAGS "/DEF:${CMAKE_CURRENT_SOURCE_DIR}/build/win64/#blst.def")
#    else()
#        if(EXISTS ${ASM_SRC} AND NOT BLST_PORTABLE)
#            list(APPEND BLST_SHARED_SRC ${ASM_SRC})
#        endif()
#        set(BLST_SHARED_OPTS -fPIC)
#        if(APPLE)
#            set(BLST_SHARED_PROPS ${BLST_SHARED_PROPS} SUFFIX ".dylib")
#        endif()
#    endif()
#
#    add_library(blst_shared SHARED ${BLST_SHARED_SRC})
#    set_target_properties(blst_shared PROPERTIES ${BLST_SHARED_PROPS})
#    target_compile_options(blst_shared PRIVATE ${BLST_SHARED_OPTS})
#endif()

# ---- Configure config.h ----
if(BLST_PORTABLE)
    set(__BLST_PORTABLE__ 1)
    set(__BLST_NO_ASM__ 1)
else()
    unset(__BLST_PORTABLE__)
    unset(__BLST_NO_ASM__)
endif()

#if(MSVC AND BLST_BUILD_SHARED)
#    set(__BLST_DLL_MAIN__ 1)
#else()
#    unset(__BLST_DLL_MAIN__)
#endif()

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
