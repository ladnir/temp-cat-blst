# Phase 4 base-field bring-up.
#
# The imported VROOM core is header-driven. This target only carries the
# dependency surface needed to compile the copied base-field layer.

set(BLST_VROOM_CORE_SOURCES)
set(BLST_VROOM_ADAPTER_SOURCES)
set(BLST_VROOM_TEST_SOURCES)
set(BLST_VROOM_BENCH_SOURCES)

add_library(blst_vroom_headers INTERFACE)

if(BLST_AVX512_PORTABLE)
    target_compile_definitions(blst_vroom_headers INTERFACE USE_AVX512_FALLBACK)
endif()

set(BLST_HAVE_VROOM_BASEFIELD OFF CACHE INTERNAL "Whether the copied VROOM base-field layer is enabled.")

if(MSVC)
    set(BLST_VROOM_VCPKG_ROOT "${CMAKE_SOURCE_DIR}/../vcpkg/installed/x64-windows")
    find_path(BLST_VROOM_GMPXX_INCLUDE_DIR
        NAMES gmpxx.h
        PATHS "${BLST_VROOM_VCPKG_ROOT}/include"
        NO_DEFAULT_PATH
    )
    find_library(BLST_VROOM_MPIR_RELEASE_LIB
        NAMES mpir
        PATHS "${BLST_VROOM_VCPKG_ROOT}/lib"
        NO_DEFAULT_PATH
    )
    find_library(BLST_VROOM_MPIR_DEBUG_LIB
        NAMES mpir
        PATHS "${BLST_VROOM_VCPKG_ROOT}/debug/lib"
        NO_DEFAULT_PATH
    )
    set(BLST_VROOM_MPIR_RELEASE_DLL "${BLST_VROOM_VCPKG_ROOT}/bin/mpir.dll")
    set(BLST_VROOM_MPIR_DEBUG_DLL "${BLST_VROOM_VCPKG_ROOT}/debug/bin/mpir.dll")

    if(BLST_VROOM_GMPXX_INCLUDE_DIR
        AND BLST_VROOM_MPIR_RELEASE_LIB
        AND BLST_VROOM_MPIR_DEBUG_LIB
        AND EXISTS "${BLST_VROOM_MPIR_RELEASE_DLL}"
        AND EXISTS "${BLST_VROOM_MPIR_DEBUG_DLL}")
        target_include_directories(blst_vroom_headers INTERFACE "${BLST_VROOM_GMPXX_INCLUDE_DIR}")
        target_link_libraries(blst_vroom_headers INTERFACE
            optimized "${BLST_VROOM_MPIR_RELEASE_LIB}"
            debug "${BLST_VROOM_MPIR_DEBUG_LIB}"
        )
        set(BLST_VROOM_MPIR_RELEASE_DLL "${BLST_VROOM_MPIR_RELEASE_DLL}" CACHE INTERNAL "Release MPIR runtime DLL for copied VROOM targets." FORCE)
        set(BLST_VROOM_MPIR_DEBUG_DLL "${BLST_VROOM_MPIR_DEBUG_DLL}" CACHE INTERNAL "Debug MPIR runtime DLL for copied VROOM targets." FORCE)
        set(BLST_HAVE_VROOM_BASEFIELD ON CACHE INTERNAL "Whether the copied VROOM base-field layer is enabled." FORCE)
    else()
        message(WARNING "VROOM base-field support disabled: MPIR headers/libs were not found under ${BLST_VROOM_VCPKG_ROOT}.")
    endif()
else()
    find_path(BLST_VROOM_GMPXX_INCLUDE_DIR NAMES gmpxx.h)
    find_library(BLST_VROOM_GMPXX_LIBRARY NAMES gmpxx mpirxx)
    find_library(BLST_VROOM_GMP_LIBRARY NAMES gmp mpir)

    if(BLST_VROOM_GMPXX_INCLUDE_DIR AND BLST_VROOM_GMPXX_LIBRARY AND BLST_VROOM_GMP_LIBRARY)
        target_include_directories(blst_vroom_headers INTERFACE "${BLST_VROOM_GMPXX_INCLUDE_DIR}")
        target_link_libraries(blst_vroom_headers INTERFACE
            "${BLST_VROOM_GMPXX_LIBRARY}"
            "${BLST_VROOM_GMP_LIBRARY}"
        )
        set(BLST_HAVE_VROOM_BASEFIELD ON CACHE INTERNAL "Whether the copied VROOM base-field layer is enabled." FORCE)
    else()
        message(WARNING "VROOM base-field support disabled: gmpxx/gmp were not found.")
    endif()
endif()
