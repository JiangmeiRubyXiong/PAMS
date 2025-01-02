#!/bin/bash

# Define paths
root_dir="/media/disk2/beijing_dti"
csv_file="${root_dir}/enhanced/unpack/BeijingEnhancedDTIProcessed.csv"
reference_file="/home/xiongj3/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz"

# Extract the NewID column (assuming it's the second column)
ids=$(seq 0 88)

# Loop through each ID
for id in $ids; do
    
    # Define the input FA file path
    input_file="${root_dir}/Neuroimaging/train_T1_FA/synthetic_brain/subject_${id}_reassembled_brain.nii.gz"
    
    # Define the output FA file path
    output_file="${root_dir}/Neuroimaging/train_T1_FA/synthetic_brain_2mm/subject_${id}_reassembled_brain.nii.gz"
    
    # Check if the input file exists
    if [ -f "$input_file" ]; then
        echo "Processing ID $id..."
        # Run the flirt command to reslice the FA file to 2mm isotropic voxels
        flirt -in "$input_file" -ref "$reference_file" -applyisoxfm 2 -out "$output_file"
    else
        echo "FA file for ID $id not found: $input_file"
    fi
done

echo "Reslicing complete!"