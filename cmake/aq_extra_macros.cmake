# (C) Copyright 2018-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Extra macros to eliminate repetition

# Macro to link list of files from source to destination
macro( LINK_FILES filelist src_dir dst_dir )
  foreach(FILENAME ${filelist})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${src_dir}/${FILENAME}
      ${dst_dir}/${FILENAME}
    )
  endforeach(FILENAME)
endmacro()

macro( LINK_FILES_DIR filelist dst_dir )
  foreach(FILENAME ${filelist})
    execute_process( COMMAND ln -sf ${FILENAME} ${dst_dir} )
  endforeach(FILENAME)
endmacro()
