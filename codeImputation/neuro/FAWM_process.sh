#!/bin/bash
# set -e
# FA_LIST="/media/disk2/beijing_dti/Neuroimaging/listBrainFiles.txt"  # Input file list
# CSV_LIST="/media/disk2/beijing_dti/Neuroimaging/listBraincsv.txt"  # Output file list
# ATLAS="/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_WMmask2.nii.gz"  # Atlas file
# 
# # Loop through each FA image in the list
# while IFS= read -r FA_IMAGE; do
#   FA_IMAGE=$(echo "$FA_IMAGE" | tr -d '"')
#   OUT_FILE=$(echo "$FA_IMAGE" | tr -d '"')
#   # Swap 'imputation_brains' with 'imputation_csv' for the output path
#   OUTPUT_IMAGE=${FA_IMAGE/imputation_brains/imputation_csv}
#   
#   # Extract the output directory path
#   OUTPUT_DIR=$(dirname "$OUTPUT_IMAGE")
# 
#   # Create the output directory if it doesn't exist
#   mkdir -p "$OUTPUT_DIR"
# 
#   # Get the base file name (without extension)
#   BASENAME=$(basename "$FA_IMAGE" .nii.gz)
# 
#   # Run c3d and save the result as a CSV in the new path
#   c3d "$FA_IMAGE" -spacing 2x2x2mm -o "$FA_IMAGE"
#   c3d "$FA_IMAGE" "$ATLAS" -lstat > "$OUTPUT_DIR/${BASENAME}.csv"
# 
#   echo "Processed: $FA_IMAGE -> $OUTPUT_DIR/${BASENAME}.csv"
# 
# done < "$FA_LIST"

# Input files
FA_IMAGE_FILE="/media/disk2/beijing_dti/Neuroimaging/listCalibFiles.txt"
OUTPUT_FILE_FILE="/media/disk2/beijing_dti/Neuroimaging/listCalibcsv.txt"

# Atlas path
ATLAS="/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_FAmask2.nii.gz" 

# Run jobs in parallel
paste "$FA_IMAGE_FILE" "$OUTPUT_FILE_FILE" | tr -d '"' | xargs -P 8 -n 2 bash -c '
    FA_IMAGE="$0"
    OUTPUT_FILE="$1"
    OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
    mkdir -p "$OUTPUT_DIR"
    c3d "$FA_IMAGE" "'"$ATLAS"'" -lstat > "$OUTPUT_FILE"
    echo "Processed $FA_IMAGE to $OUTPUT_FILE"
'