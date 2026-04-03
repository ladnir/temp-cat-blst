# Phase 4 base-field bring-up.
#
# The promoted AVX backend is header-driven. This target only carries the
# compile definitions and optional validation-only GMP linkage.

add_library(blst_avx_headers INTERFACE)

set(BLST_HAVE_AVX_BASEFIELD OFF CACHE INTERNAL "Whether the copied AVX base-field layer is enabled.")
set(BLST_HAVE_AVX_VALIDATION_GMP OFF CACHE INTERNAL "Whether GMP/MPIR is available for AVX constant validation.")

if(BLST_AVX_VALIDATE_CONSTANTS)
    target_compile_definitions(blst_avx_headers INTERFACE BLST_AVX_VALIDATE_CONSTANTS=1)
endif()

if(BLST_AVX512 OR BLST_AVX512_PORTABLE)
    set(BLST_HAVE_AVX_BASEFIELD ON CACHE INTERNAL "Whether the copied AVX base-field layer is enabled." FORCE)

    if(BLST_AVX_VALIDATE_CONSTANTS)
        if(MSVC)
            set(BLST_AVX_VCPKG_ROOT "${CMAKE_SOURCE_DIR}/../vcpkg/installed/x64-windows")
            find_path(BLST_AVX_GMPXX_INCLUDE_DIR
                NAMES gmpxx.h
                PATHS "${BLST_AVX_VCPKG_ROOT}/include"
                NO_DEFAULT_PATH
            )
            find_library(BLST_AVX_MPIR_RELEASE_LIB
                NAMES mpir
                PATHS "${BLST_AVX_VCPKG_ROOT}/lib"
                NO_DEFAULT_PATH
            )
            find_library(BLST_AVX_MPIR_DEBUG_LIB
                NAMES mpir
                PATHS "${BLST_AVX_VCPKG_ROOT}/debug/lib"
                NO_DEFAULT_PATH
            )
            set(BLST_AVX_MPIR_RELEASE_DLL "${BLST_AVX_VCPKG_ROOT}/bin/mpir.dll")
            set(BLST_AVX_MPIR_DEBUG_DLL "${BLST_AVX_VCPKG_ROOT}/debug/bin/mpir.dll")

            if(BLST_AVX_GMPXX_INCLUDE_DIR
                AND BLST_AVX_MPIR_RELEASE_LIB
                AND BLST_AVX_MPIR_DEBUG_LIB
                AND EXISTS "${BLST_AVX_MPIR_RELEASE_DLL}"
                AND EXISTS "${BLST_AVX_MPIR_DEBUG_DLL}")
                target_include_directories(blst_avx_headers INTERFACE "${BLST_AVX_GMPXX_INCLUDE_DIR}")
                target_link_libraries(blst_avx_headers INTERFACE
                    optimized "${BLST_AVX_MPIR_RELEASE_LIB}"
                    debug "${BLST_AVX_MPIR_DEBUG_LIB}"
                )
                set(BLST_AVX_MPIR_RELEASE_DLL "${BLST_AVX_MPIR_RELEASE_DLL}" CACHE INTERNAL "Release MPIR runtime DLL for copied AVX targets." FORCE)
                set(BLST_AVX_MPIR_DEBUG_DLL "${BLST_AVX_MPIR_DEBUG_DLL}" CACHE INTERNAL "Debug MPIR runtime DLL for copied AVX targets." FORCE)
                set(BLST_HAVE_AVX_VALIDATION_GMP ON CACHE INTERNAL "Whether GMP/MPIR is available for AVX constant validation." FORCE)
            else()
                message(FATAL_ERROR "BLST_AVX_VALIDATE_CONSTANTS=ON requires MPIR headers/libs under ${BLST_AVX_VCPKG_ROOT}.")
            endif()
        else()
            find_path(BLST_AVX_GMPXX_INCLUDE_DIR NAMES gmpxx.h)
            find_library(BLST_AVX_GMPXX_LIBRARY NAMES gmpxx mpirxx)
            find_library(BLST_AVX_GMP_LIBRARY NAMES gmp mpir)

            if(BLST_AVX_GMPXX_INCLUDE_DIR AND BLST_AVX_GMPXX_LIBRARY AND BLST_AVX_GMP_LIBRARY)
                target_include_directories(blst_avx_headers INTERFACE "${BLST_AVX_GMPXX_INCLUDE_DIR}")
                target_link_libraries(blst_avx_headers INTERFACE
                    "${BLST_AVX_GMPXX_LIBRARY}"
                    "${BLST_AVX_GMP_LIBRARY}"
                )
                set(BLST_HAVE_AVX_VALIDATION_GMP ON CACHE INTERNAL "Whether GMP/MPIR is available for AVX constant validation." FORCE)
            else()
                message(FATAL_ERROR "BLST_AVX_VALIDATE_CONSTANTS=ON requires gmpxx/gmp to be installed.")
            endif()
        endif()
    endif()
endif()
