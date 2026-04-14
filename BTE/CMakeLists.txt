set(BTE_BUILDING_BLOCK_SOURCES
    "building_blocks/blst_cpp.cpp"
    "building_blocks/bls_ots.cpp"
    "building_blocks/hkdf_sha256.cpp"
    "building_blocks/kzg.cpp"
    "building_blocks/lagrange.cpp"
    "building_blocks/polynomial.cpp"
    "building_blocks/serialization.cpp"
)

add_library(bte_building_blocks STATIC ${BTE_BUILDING_BLOCK_SOURCES})
target_include_directories(bte_building_blocks PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
target_link_libraries(bte_building_blocks PUBLIC blst_static)

target_compile_options(bte_building_blocks PUBLIC
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:MSVC>>:/std:c++latest>
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<NOT:$<CXX_COMPILER_ID:MSVC>>>:-std=c++23>
)

if(TARGET blst_sanitizers)
    target_link_libraries(bte_building_blocks PUBLIC blst_sanitizers)
endif()

add_library(bte_trx STATIC "trx/core.cpp")
target_include_directories(bte_trx PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
target_link_libraries(bte_trx PUBLIC bte_building_blocks)
target_compile_options(bte_trx PUBLIC
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:MSVC>>:/std:c++latest>
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<NOT:$<CXX_COMPILER_ID:MSVC>>>:-std=c++23>
)

if(TARGET blst_sanitizers)
    target_link_libraries(bte_trx PUBLIC blst_sanitizers)
endif()

add_library(bte_btx STATIC
    "btx/core.cpp"
    "btx/serialization.cpp"
)
target_include_directories(bte_btx PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
target_link_libraries(bte_btx PUBLIC bte_building_blocks)
target_compile_options(bte_btx PUBLIC
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:MSVC>>:/std:c++latest>
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<NOT:$<CXX_COMPILER_ID:MSVC>>>:-std=c++23>
)

if(TARGET blst_sanitizers)
    target_link_libraries(bte_btx PUBLIC blst_sanitizers)
endif()
