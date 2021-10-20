#!/usr/bin/env bash

# Loop over input files
for item in "$@"; do
  if [ -e Data/${item} ]; then
    echo "Item ${item} already present"
  else
    echo "Linking ${item}: ln -sf ${LOCAL_PATH_AQ_JEDI_TESTFILES}/aq-tests-tier2/${item} Data/${item}"
    ln -sf ${LOCAL_PATH_AQ_JEDI_TESTFILES}/aq-tests-tier2/${item} Data/${item}
  fi
done

exit 0
