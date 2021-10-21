#!/usr/bin/env bash

# Loop over input files
for item in "$@"; do
  if [ -e Data/${item} ]; then
    echo "Item ${item} already present"
  else
    echo "Linking ${item}: ln -sf ${AQ_JEDI_LOCAL_PATH_TESTFILES}/aq-tests/${item} Data/${item}"
    ln -sf ${AQ_JEDI_LOCAL_PATH_TESTFILES}/aq-tests/${item} Data/${item}
  fi
done

exit 0
