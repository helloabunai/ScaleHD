#!/bin/bash

#run_scalehd_test.sh

set -euo pipefail

echo "Running ScaleHD demo tests..."

echo "" > ContinuousTesting/ScaleHDLog.txt
# Run ScaleHD on example demo BAM and reference FASTA
# Adjust file names/paths according to the actual demo data in ScaleHD repo
ScaleHD -v -b -c "ContinuousTesting/ConfigurationFile_MiSeq_runX.xml" -o ContinuousTesting -j MiSeq_runX_ScaleHDresults

echo "Validating checksums..."

# Generate actual checksums
find ContinuousTesting/MiSeq_runX_ScaleHDresults -maxdepth 1 -type f -name "*.csv" \
  -exec sha256sum {} \; \
  | awk '{split($2, a, "/"); print a[length(a)] "," $1}' > ContinuousTesting/actual_checksums.csv

# Compare against expected
diff_output=$(diff <(sort ContinuousTesting/scalehd_test_expected_checksums.csv) <(sort ContinuousTesting/actual_checksums.csv) || true)

if [[ -z "$diff_output" ]]; then
    echo "All checksums match. Demo test passed."
else
    echo "Checksum mismatch found:"
    echo "$diff_output"
    exit 1
fi

