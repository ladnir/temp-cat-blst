set(BTE_TRX_SOURCES
    "building_blocks/blst_cpp.cpp"
    "building_blocks/bls_ots.cpp"
    "building_blocks/hkdf_sha256.cpp"
    "building_blocks/kzg.cpp"
    "building_blocks/lagrange.cpp"
    "building_blocks/polynomial.cpp"
    "building_blocks/serialization.cpp"
    "trx/core.cpp"
)

add_library(bte_trx STATIC ${BTE_TRX_SOURCES})
target_include_directories(bte_trx PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
target_link_libraries(bte_trx PUBLIC blst_static)

target_compile_options(bte_trx PUBLIC
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:MSVC>>:/std:c++latest>
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<NOT:$<CXX_COMPILER_ID:MSVC>>>:-std=c++23>
)

if(TARGET blst_sanitizers)
    target_link_libraries(bte_trx PUBLIC blst_sanitizers)
endif()
