#!/usr/bin/env bash

# Check for the curl command
if ! command -v curl &> /dev/null; then
    echo "Command curl could not be found"
    exit 1
fi

# Loop over input files
for file in "$@"; do
  mkdir -p `dirname Data/${file}`
  if [ -f Data/${file} ]; then
    echo "File ${file} already present"
  else
    echo "Downloading ${file}: curl -u anonymous: \"ftp://ftp.cerfacs.fr/globc/andrea/oops-testdata/aq-tests-tier2/${file}\" -o Data/${file}"
    curl -u anonymous: "ftp://ftp.cerfacs.fr/globc/andrea/oops-testdata/aq-tests-tier2/${file}" -o Data/${file}
  fi
done

exit 0
