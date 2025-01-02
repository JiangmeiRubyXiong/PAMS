#!/bin/bash

# Base directory
base_dir="/media/disk2/beijing_dti/anatNew"

# Standard templates
ref_img="/home/xiongj3/fsl/data/standard/MNI152_T1_1mm.nii.gz"
ref_brain_img="/home/xiongj3/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz"

# Function to process one subject
process_subject() {
    subject_dir="$1"
    subject_id=$(basename "$subject_dir")
    input_image="$subject_dir/anat/NIfTI/defaced_MPRAGE.nii.gz"

    # Check if the input image exists
    if [[ -f "$input_image" ]]; then
        echo "Processing subject: $subject_id"
        
        # Create a directory for the subject
        output_dir="./${subject_id}"
        mkdir -p "$output_dir"

        # Pre-alignment
        flirt -in "$input_image" \
              -ref "$ref_img" \
              -out "$output_dir/prealigned.nii.gz"
        
        # Bias correction
        fast -B --type=1 "$output_dir/prealigned.nii.gz"
        
        # Skull stripping
        bet "$output_dir/prealigned_restore.nii.gz" \
            "$output_dir/brain_extracted.nii.gz" -f 0.4 -R
        
        # Linear registration
        flirt -in "$output_dir/brain_extracted.nii.gz" \
              -ref "$ref_brain_img" \
              -out "$output_dir/linear_transform.nii.gz" \
              -omat "$output_dir/linear_transform.mat"
        
        # Nonlinear transformation
        fnirt --in="$output_dir/brain_extracted.nii.gz" \
              --ref="$ref_brain_img" \
              --aff="$output_dir/linear_transform.mat" \
              --iout="$output_dir/final_output.nii.gz" \
              --cout="$output_dir/nonlinear_warp_coefficients.nii.gz"

        echo "Finished processing subject: $subject_id"
    else
        echo "Input image not found for subject: $subject_id"
    fi
}

export -f process_subject

# Find all subject directories and process them in parallel
for subject_dir in "$base_dir"/*; do
    process_subject "$subject_dir" &
    while [[ $(jobs -r -p | wc -l) -ge 16 ]]; do
        wait -n
    done
done
wait