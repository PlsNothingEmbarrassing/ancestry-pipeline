#!/bin/bash

# Define variables for file paths
TEST=/home/c.c24082291/PLINK_Bioinfo/database/genotypes/unknown-qc-v9/unknown-qc-v9
OUTDIR=/home/c.c24082291/PLINK_Bioinfo/project/PCA_ASSESSMENT/tmp
RESULTS=/home/c.c24082291/PLINK_Bioinfo/project/PCA_ASSESSMENT/results

# Load the plink module
module load plink/1.9

# Check if the user provided a keep file
if [ -z "$1" ]; then
    echo "Error: No .keep file specified."
    exit 1
fi

# Assign the first argument to KEEP_FILE
KEEP_FILE=$1

# Check if the keep file exists
if [ ! -f "$KEEP_FILE" ]; then
    echo "Error: The specified .keep file does not exist."
    exit 1
fi

# Check if the keep file is empty
if [ ! -s "$KEEP_FILE" ]; then
    echo "Error: The specified .keep file is empty."
    exit 1
fi

# Check if the keep file has the correct extension
if [[ "$KEEP_FILE" != *.keep ]]; then
    echo "Error: The specified file is not in the correct format."
    exit 1
fi

# Run plink to filter the dataset based on the keep file and create a new binary file
plink --bfile $TEST \
--keep "$KEEP_FILE" \
--make-bed \
--out "${OUTDIR}/test_gbr"

# Split the X chromosome pseudo-autosomal region
plink --bfile "${OUTDIR}/test_gbr" \
--split-x b37 no-fail \
--make-bed \
--out "${OUTDIR}/test_gbr_split"

# Run plink to check the sex of individuals in the filtered dataset
plink --bfile "${OUTDIR}/test_gbr_split" \
--check-sex \
--out "${OUTDIR}/test_gbr_sexcheck"

# Extract individuals identified as female (sex code 2) from the sex check results
awk '$4 == 2' "${OUTDIR}/test_gbr_sexcheck.sexcheck" > "${OUTDIR}/test_gbr_female.sexcheck"