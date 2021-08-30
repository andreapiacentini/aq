add_subdirectory( aq )
add_subdirectory( mains )

ecbuild_add_test( TARGET aq_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )
