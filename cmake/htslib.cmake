# Indexed FASTA needs HTSlib's local faidx/BGZF path and zlib only.
include(ExternalProject)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
find_program(MAKE_COMMAND NAMES gmake make REQUIRED)
find_program(SH_COMMAND NAMES sh bash REQUIRED)

set(HTSLIB_SRC_DIR "${CMAKE_SOURCE_DIR}/third_party/htslib")
set(HTSLIB_CFLAGS "${CMAKE_C_FLAGS} -O2 -fPIC -ffunction-sections -fdata-sections")
set(HTSLIB_CPPFLAGS "")
set(HTSLIB_LDFLAGS "${CMAKE_EXE_LINKER_FLAGS}")
if(APPLE AND DEFINED OSX_BUILD_ARCH AND NOT "${OSX_BUILD_ARCH}" STREQUAL "")
    string(APPEND HTSLIB_CFLAGS " -arch ${OSX_BUILD_ARCH}")
    string(APPEND HTSLIB_LDFLAGS " -arch ${OSX_BUILD_ARCH}")
endif()
if(DEFINED VCPKG_INSTALLED_DIR AND DEFINED VCPKG_TARGET_TRIPLET)
    set(_vcpkg_prefix "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}")
    if(EXISTS "${_vcpkg_prefix}")
        string(APPEND HTSLIB_CPPFLAGS " -I${_vcpkg_prefix}/include")
        string(APPEND HTSLIB_LDFLAGS " -L${_vcpkg_prefix}/lib")
    endif()
endif()
if(DUCKDB_WASM_EXTENSION)
    string(APPEND HTSLIB_CFLAGS " -fwasm-exceptions -s USE_ZLIB=1")
    string(APPEND HTSLIB_LDFLAGS " -s USE_ZLIB=1")
    string(APPEND HTSLIB_CPPFLAGS " -DDUCKVEP_WASM_DUCKDB_RUNTIME=1 -include ${CMAKE_SOURCE_DIR}/src/include/wasm_socket_compat.h")
endif()

set(HTSLIB_CONFIGURE_FLAGS
    --disable-libcurl --disable-plugins --disable-s3 --disable-gcs
    --disable-bz2 --disable-lzma --without-libdeflate)
set(HTSLIB_CONFIGURE_ENV_VARS "")
if(DUCKDB_WASM_EXTENSION)
    list(APPEND HTSLIB_CONFIGURE_ENV_VARS
        "ac_cv_search_recv=none required"
        "ac_cv_lib_z_inflate=yes"
        "ac_cv_func_fork=no" "ac_cv_func_vfork=no"
        "ac_cv_func_getrandom=no")
    list(APPEND HTSLIB_CONFIGURE_FLAGS --host=wasm32-unknown-emscripten)
    if(DEFINED HTSLIB_BUILD_TRIPLET AND NOT "${HTSLIB_BUILD_TRIPLET}" STREQUAL "")
        list(APPEND HTSLIB_CONFIGURE_FLAGS "--build=${HTSLIB_BUILD_TRIPLET}")
    endif()
    list(APPEND HTSLIB_CONFIGURE_ENV_VARS ${HTSLIB_AUTOCONF_CACHE})
endif()
set(HTSLIB_MAKE_ARGS
    "CC=${CMAKE_C_COMPILER}"
    "CFLAGS=${HTSLIB_CFLAGS}"
    "CPPFLAGS=${HTSLIB_CPPFLAGS}"
    "LDFLAGS=${HTSLIB_LDFLAGS}"
    "LIBS=-lz")

ExternalProject_Add(htslib_build
    SOURCE_DIR "${HTSLIB_SRC_DIR}"
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env
        "CC=${CMAKE_C_COMPILER}"
        "CFLAGS=${HTSLIB_CFLAGS}"
        "CPPFLAGS=${HTSLIB_CPPFLAGS}"
        "LDFLAGS=${HTSLIB_LDFLAGS}"
        "LIBS=-lz"
        ${HTSLIB_CONFIGURE_ENV_VARS}
        ${SH_COMMAND} ./configure ${HTSLIB_CONFIGURE_FLAGS}
    BUILD_COMMAND ${MAKE_COMMAND} -j lib-static ${HTSLIB_MAKE_ARGS}
    BUILD_ALWAYS TRUE
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS "${HTSLIB_SRC_DIR}/libhts.a"
    LOG_BUILD TRUE
)

# HTSlib's Windows objects export every API symbol; the extension only exports
# its DuckDB entrypoint.
if(WIN32)
    if(NOT CMAKE_OBJCOPY)
        message(FATAL_ERROR "objcopy is required to strip htslib's DLL exports on Windows")
    endif()
    ExternalProject_Add_Step(htslib_build strip_dllexport
        COMMAND "${CMAKE_OBJCOPY}" --remove-section=.drectve "${HTSLIB_SRC_DIR}/libhts.a"
        DEPENDEES build
        ALWAYS TRUE
    )
endif()
add_library(hts STATIC IMPORTED GLOBAL)
set_target_properties(hts PROPERTIES
    IMPORTED_LOCATION "${HTSLIB_SRC_DIR}/libhts.a"
    INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_SRC_DIR}")
add_dependencies(hts htslib_build)
set(HTSLIB_LINK_LIBS hts ZLIB::ZLIB Threads::Threads)
if(WIN32 AND MINGW)
    find_library(GNUREGEX_LIBRARY NAMES gnurx regex)
    if(GNUREGEX_LIBRARY)
        list(APPEND HTSLIB_LINK_LIBS ${GNUREGEX_LIBRARY})
    endif()
    list(APPEND HTSLIB_LINK_LIBS ws2_32)
elseif(NOT APPLE AND NOT DUCKDB_WASM_EXTENSION)
    list(APPEND HTSLIB_LINK_LIBS m)
endif()
message(STATUS "htslib link libs: ${HTSLIB_LINK_LIBS}")
