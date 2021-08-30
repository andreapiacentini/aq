# Link scripts
list( APPEND tools_files
    aq_data_downloader.sh
)

foreach( FILENAME ${tools_files} )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                             ${PROJECT_SOURCE_DIR}/aq/tools/${FILENAME}
                             ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
endforeach()
