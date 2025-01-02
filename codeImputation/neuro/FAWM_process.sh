#!/bin/bash

# Define your FA images directory and the output directory
FA_DIR=~/ImageImputation/Neuroimaging/train_T1_FA/synthetic_1mm
OUTPUT_DIR=~/ImageImputation/Neuroimaging/train_T1_FA/output/stats
ATLAS="/home/xiongj3/fsl/data/atlases/JHU/JHU-ICBM-tracts-maxprob-thr0-1mm.nii.gz"  # The JHU EVE atlas

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop over all FA images matching the pattern subject_*_reassembled_brain.nii.gz
for FA_IMAGE in "$FA_DIR"/subject_*_reassembled_brain.nii.gz; do
  # Ensure the file exists
  if [ ! -e "$FA_IMAGE" ]; then
    echo "No matching files found for pattern: $FA_IMAGE"
    continue
  fi

  # Get the base name of the FA image (e.g., subject_0)
  SUBJECT_NAME=$(basename "$FA_IMAGE" "_reassembled_brain.nii.gz")

  # Extract FA statistics for each subject
  c3d "$FA_IMAGE" "$ATLAS" -lstat > "$OUTPUT_DIR/${SUBJECT_NAME}_FA_stats_1mm.csv"
done


