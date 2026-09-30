# DuckVEP v2 host (preview): a self-contained build of the stable C API v2
# extension. It is deliberately a separate CMake project so the v1 build
# (../../CMakeLists.txt, the distributed extension) is not touched by it.
#
#   cmake -S host_v2 -B build/cmake_v2 -DCMAKE_BUILD_TYPE=Release
#   cmake --build build/cmake_v2
#
# `make release_v2` from the repository root runs these and appends the footer.
cmake_minimum_required(VERSION 3.5...3.29)
project(duckvep_v2 LANGUAGES C)

get_filename_component(DUCKVEP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
set(DUCKDB_CAPI_V2_DIR "${DUCKVEP_ROOT}/duckdb_capi_v2" CACHE PATH "Pinned DuckDB C API v2 headers")
foreach(header duckdb_v2.h duckdb_extension_v2.h)
    if(NOT EXISTS "${DUCKDB_CAPI_V2_DIR}/${header}")
        message(FATAL_ERROR "Missing ${header}; run scripts/fetch-v2-sdk.py --refresh")
    endif()
endforeach()

file(GLOB KERNEL_SOURCES CONFIGURE_DEPENDS "${DUCKVEP_ROOT}/src/kernel/src/*.c")
set(DUCKVEP_SOURCE_ROOT "${DUCKVEP_ROOT}")
file(GLOB CORE_SOURCES CONFIGURE_DEPENDS "${DUCKVEP_ROOT}/src/core/*.c")
add_library(duckvep SHARED
    host_v2.c
    host_v2_nested.c
    host_v2_builders.c
    host_v2_model.c
    host_v2_annotate.c
    host_v2_haplotypes.c
    host_v2_coding_calls.c
    host_v2_budget.c
    ${CORE_SOURCES}
    ${KERNEL_SOURCES}
    "${DUCKVEP_ROOT}/src/duckvep_reference.c"
    "${DUCKVEP_ROOT}/third_party/cgranges/cgranges.c")
set_source_files_properties("${DUCKVEP_ROOT}/third_party/cgranges/cgranges.c" PROPERTIES COMPILE_OPTIONS "-w")
# The v1 build names its artifact duckvep too; this one lives in its own tree.
set_target_properties(duckvep PROPERTIES PREFIX "lib" C_VISIBILITY_PRESET hidden)
target_compile_definitions(duckvep PRIVATE DUCKDB_EXTENSION_NAME=duckvep)
target_include_directories(duckvep PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}" "${DUCKVEP_ROOT}/src" "${DUCKVEP_ROOT}/src/kernel/include" "${DUCKVEP_ROOT}/src/kernel/src"
    "${DUCKVEP_ROOT}/third_party/cgranges")
target_include_directories(duckvep SYSTEM PRIVATE "${DUCKDB_CAPI_V2_DIR}")
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(duckvep PRIVATE -Wall -Wextra -Wpedantic -Wformat=2 -Wstrict-prototypes
        "-ffile-prefix-map=${DUCKVEP_ROOT}=." -ffunction-sections -fdata-sections -Werror)
endif()
# Indexed-FASTA validation of the model needs HTSlib (local faidx and zlib only), as in v1.
include("${DUCKVEP_ROOT}/cmake/htslib.cmake")
target_link_libraries(duckvep PRIVATE ${HTSLIB_LINK_LIBS} m)
if(APPLE)
    target_link_options(duckvep PRIVATE "-Wl,-dead_strip"
        "-Wl,-exported_symbols_list,${CMAKE_CURRENT_SOURCE_DIR}/duckvep_v2.exports")
elseif(UNIX)
    target_link_options(duckvep PRIVATE "-Wl,--gc-sections" "-Wl,--exclude-libs,ALL"
        "-Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/duckvep_v2.map")
endif()
