if(DUCKDB_WASM_EXTENSION)
    # DuckDB-Wasm packages the extension archive. Merge htslib and its
    # static dependencies into that archive.
    set(DUCKVEP_WASM_FAT_ARCHIVE_INPUTS "${HTSLIB_BUILD_DIR}/libhts.a")

    if(TARGET ZLIB::ZLIB)
        list(APPEND DUCKVEP_WASM_FAT_ARCHIVE_INPUTS "$<TARGET_FILE:ZLIB::ZLIB>")
    elseif(ZLIB_LIBRARIES)
        list(APPEND DUCKVEP_WASM_FAT_ARCHIVE_INPUTS ${ZLIB_LIBRARIES})
    endif()

    if(BZIP2_FOUND)
        if(TARGET BZip2::BZip2)
            list(APPEND DUCKVEP_WASM_FAT_ARCHIVE_INPUTS "$<TARGET_FILE:BZip2::BZip2>")
        elseif(BZIP2_LIBRARIES)
            list(APPEND DUCKVEP_WASM_FAT_ARCHIVE_INPUTS ${BZIP2_LIBRARIES})
        endif()
    endif()

    if(LIBLZMA_FOUND)
        if(TARGET LibLZMA::LibLZMA)
            list(APPEND DUCKVEP_WASM_FAT_ARCHIVE_INPUTS "$<TARGET_FILE:LibLZMA::LibLZMA>")
        elseif(LIBLZMA_LIBRARIES)
            list(APPEND DUCKVEP_WASM_FAT_ARCHIVE_INPUTS ${LIBLZMA_LIBRARIES})
        endif()
    endif()

    if(LIBDEFLATE_FOUND)
        list(APPEND DUCKVEP_WASM_FAT_ARCHIVE_INPUTS "${LIBDEFLATE_LIBRARY}")
    endif()

    if(CURL_FOUND AND NOT DUCKDB_WASM_EXTENSION)
        if(WIN32 AND MINGW)
            list(APPEND DUCKVEP_WASM_FAT_ARCHIVE_INPUTS ${PC_LIBCURL_LDFLAGS_LIST})
        elseif(TARGET CURL::libcurl)
            list(APPEND DUCKVEP_WASM_FAT_ARCHIVE_INPUTS "$<TARGET_FILE:CURL::libcurl>")
        endif()
    endif()

    if(OPENSSL_FOUND AND NOT DUCKDB_WASM_EXTENSION AND TARGET OpenSSL::Crypto)
        list(APPEND DUCKVEP_WASM_FAT_ARCHIVE_INPUTS "$<TARGET_FILE:OpenSSL::Crypto>")
    endif()

    if(WIN32 AND MINGW AND GNUREGEX_LIBRARY)
        list(APPEND DUCKVEP_WASM_FAT_ARCHIVE_INPUTS "${GNUREGEX_LIBRARY}")
    endif()

    list(REMOVE_DUPLICATES DUCKVEP_WASM_FAT_ARCHIVE_INPUTS)
    string(JOIN "\n" DUCKVEP_WASM_FAT_ARCHIVE_INPUTS_CONTENT ${DUCKVEP_WASM_FAT_ARCHIVE_INPUTS})
    file(GENERATE OUTPUT "${CMAKE_BINARY_DIR}/wasm_fat_archive_inputs.txt"
        CONTENT "${DUCKVEP_WASM_FAT_ARCHIVE_INPUTS_CONTENT}\n")

    add_dependencies(${EXTENSION_NAME} htslib_build)
    add_custom_command(TARGET ${EXTENSION_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND}
            -DARCHIVE=$<TARGET_FILE:${EXTENSION_NAME}>
            -DINPUT_LIST=${CMAKE_BINARY_DIR}/wasm_fat_archive_inputs.txt
            -DARCHIVER=${CMAKE_AR}
            -DRANLIB=${CMAKE_RANLIB}
            -P ${CMAKE_SOURCE_DIR}/cmake/merge_static_archives.cmake
        VERBATIM)
endif()
