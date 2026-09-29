cmake_minimum_required(VERSION 3.5...3.29)
option(DUCKDB_WASM_EXTENSION "Build for DuckDB-Wasm" OFF)
option(DUCKVEP_WEBR_EXTENSION "Build a webR Emscripten side module" OFF)

if(NOT DEFINED EXTENSION_NAME)
    message(FATAL_ERROR "DuckDB extension name is required")
endif()
add_definitions(-DDUCKDB_EXTENSION_NAME=${EXTENSION_NAME})

if(APPLE)
    if(NOT DEFINED OSX_BUILD_ARCH OR "${OSX_BUILD_ARCH}" STREQUAL "")
        if(DEFINED VCPKG_TARGET_TRIPLET)
            if(VCPKG_TARGET_TRIPLET MATCHES "^x64-osx")
                set(OSX_BUILD_ARCH "x86_64")
            elseif(VCPKG_TARGET_TRIPLET MATCHES "^arm64-osx")
                set(OSX_BUILD_ARCH "arm64")
            endif()
        endif()
        if((NOT DEFINED OSX_BUILD_ARCH OR "${OSX_BUILD_ARCH}" STREQUAL "") AND
           DEFINED CMAKE_OSX_ARCHITECTURES AND NOT "${CMAKE_OSX_ARCHITECTURES}" STREQUAL "")
            list(GET CMAKE_OSX_ARCHITECTURES 0 OSX_BUILD_ARCH)
        endif()
    endif()
    if(DEFINED OSX_BUILD_ARCH AND NOT "${OSX_BUILD_ARCH}" STREQUAL "")
        set(CMAKE_OSX_ARCHITECTURES "${OSX_BUILD_ARCH}" CACHE STRING "" FORCE)
    endif()
endif()

if(DEFINED DUCKDB_EXTENSION_API_VERSION_UNSTABLE)
    add_definitions(-DDUCKDB_EXTENSION_API_VERSION_UNSTABLE=${DUCKDB_EXTENSION_API_VERSION_UNSTABLE})
    if(DEFINED TARGET_DUCKDB_VERSION_MAJOR)
        add_definitions(-DDUCKDB_EXTENSION_API_VERSION_MAJOR=${TARGET_DUCKDB_VERSION_MAJOR})
        add_definitions(-DDUCKDB_EXTENSION_API_VERSION_MINOR=${TARGET_DUCKDB_VERSION_MINOR})
        add_definitions(-DDUCKDB_EXTENSION_API_VERSION_PATCH=${TARGET_DUCKDB_VERSION_PATCH})
    endif()
endif()

project(${EXTENSION_NAME} LANGUAGES C)
if(DUCKVEP_WEBR_EXTENSION)
    # Emscripten side modules are loadable shared libraries in webR.
    set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)
endif()
file(GLOB KERNEL_SOURCES CONFIGURE_DEPENDS src/kernel/src/*.c)
set(EXTENSION_SOURCES
    src/duckvep.c
    src/duckvep_registration.c
    src/duckvep_builder.c
    src/duckvep_prepare_sql.c
    src/duckvep_structural_sql.c
    src/duckvep_repeat.c
    src/duckvep_annotate.c
    src/duckvep_ensembl.c
    src/duckvep_haplotype_sql.c
    src/duckvep_model.c
    src/duckvep_phase_sql.c
    src/duckvep_reference.c
    src/duckvep_sql.c
    third_party/cgranges/cgranges.c
    ${KERNEL_SOURCES}
)
if(DUCKDB_WASM_EXTENSION AND NOT DUCKVEP_WEBR_EXTENSION)
    add_library(${EXTENSION_NAME} STATIC ${EXTENSION_SOURCES})
else()
    add_library(${EXTENSION_NAME} SHARED ${EXTENSION_SOURCES})
endif()
if(DUCKDB_WASM_EXTENSION)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE DUCKVEP_WASM_DUCKDB_RUNTIME=1)
    target_compile_options(${EXTENSION_NAME} PRIVATE -fwasm-exceptions -include ${CMAKE_SOURCE_DIR}/src/include/wasm_socket_compat.h)
endif()
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(VISIBILITY_INLINES_HIDDEN ON)
if(NOT APPLE AND NOT DUCKDB_WASM_EXTENSION)
    set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -s")
endif()
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(${EXTENSION_NAME} PRIVATE -Wall -Wextra -Wpedantic -Wformat=2 -Wstrict-prototypes "-ffile-prefix-map=${CMAKE_SOURCE_DIR}=." -ffunction-sections -fdata-sections)
    if(NOT EMSCRIPTEN)
        target_compile_options(${EXTENSION_NAME} PRIVATE -Werror)
    endif()
    set_source_files_properties(third_party/cgranges/cgranges.c PROPERTIES COMPILE_OPTIONS "-w")
endif()
target_include_directories(${EXTENSION_NAME} PRIVATE src/include src src/kernel/include src/kernel/src third_party/cgranges)
target_include_directories(${EXTENSION_NAME} SYSTEM PRIVATE duckdb_capi)
include(cmake/htslib.cmake)
target_link_libraries(${EXTENSION_NAME} PRIVATE ${HTSLIB_LINK_LIBS})
if(DUCKVEP_WEBR_EXTENSION)
    target_link_options(${EXTENSION_NAME} PRIVATE -fwasm-exceptions -sSIDE_MODULE=1 -Wl,--export=duckvep_init_c_api)
elseif(DUCKDB_WASM_EXTENSION)
    include(cmake/wasm_archive.cmake)
elseif(APPLE)
    target_link_options(${EXTENSION_NAME} PRIVATE
        "-Wl,-dead_strip"
        "-Wl,-exported_symbols_list,${CMAKE_SOURCE_DIR}/cmake/duckvep.exports")
elseif(UNIX)
    target_link_options(${EXTENSION_NAME} PRIVATE
        "-Wl,--gc-sections"
        "-Wl,--exclude-libs,ALL"
        "-Wl,--version-script=${CMAKE_SOURCE_DIR}/cmake/duckvep.map")
elseif(WIN32)
    target_link_options(${EXTENSION_NAME} PRIVATE "-Wl,--gc-sections")
endif()
