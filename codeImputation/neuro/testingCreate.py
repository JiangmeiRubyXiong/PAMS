import pandas as pd
import nibabel as nib
import numpy as np
import os

# File paths
csv_path = "/media/disk2/beijing_dti/enhanced/unpack/BeijingEnhancedDTIProcessed1.csv"
t1_folder = "/media/disk2/beijing_dti/enhanced/T1_registered_1mm"
fa_folder = "/media/disk2/beijing_dti/enhanced/unpack"

# Read the subject IDs from the CSV file and ensure they are strings with leading zeros
subject_data = pd.read_csv(csv_path)
subject_data['NewID'] = subject_data['NewID'].astype(str).str.zfill(3)

# Randomly sample 89 subjects
n_samples = 89
sampled_subject_data = subject_data.sample(n=n_samples, random_state=42)  # Set random_state for reproducibility

unsampled = pd.concat([sampled_subject_data,subject_data]).drop_duplicates(keep=False)

sampled_subject_data.to_csv("/media/disk2/beijing_dti/processed/train_df.csv", index=False)
unsampled.to_csv('/media/disk2/beijing_dti/processed/test_df.csv', index=False)



for slicei in range(91, 182):
  # Initialize lists to hold the T1 and FA slices
  t1_slices = []
  fa_slices = []
  for index, row in unsampled.iterrows():
      # Paths for T1 and FA files
      subject_id = row.iat[0]
      new_id = row[1]
      t1_path = os.path.join(t1_folder, f"{subject_id}", "final_output.nii.gz")
      fa_path = os.path.join(fa_folder, f"{new_id}", "scalars_standard", "DTI_FA.nii.gz")
      
      # Check if both files exist
      if not os.path.exists(t1_path):
          print(f"T1 file not found for subject {subject_id}: {t1_path}")
          continue
      
         # Load T1 and FA NIfTI files
      t1_img = nib.load(t1_path)
      fa_img = nib.load(fa_path)
      
      # Get data arrays
      t1_data = t1_img.get_fdata()
      fa_data = fa_img.get_fdata()
      
      # Extract one slice (e.g., the middle slice along the third axis)
      # slice_idx_t1 = t1_data.shape[2] // 2  # middle slice
      # slice_idx_fa = fa_data.shape[2] // 2  # middle slice (assuming same orientation)
      t1_slice = t1_data[:, :, slicei]  # Extract T1 slice
      fa_slice = fa_data[:, :, slicei]  # Extract FA slice
      
      # Expand dimensions to [X, Y, 1]
      t1_slice = np.expand_dims(t1_slice, axis=2)
      fa_slice = np.expand_dims(fa_slice, axis=2)
      
      # Append slices to their respective lists
      t1_slices.append(t1_slice)
      fa_slices.append(fa_slice)
  
  # Stack all slices along the fourth dimension
  t1_training_data = np.stack(t1_slices, axis=3)
  fa_training_data = np.stack(fa_slices, axis=3)
  
  print(f"T1 training data shape: {t1_training_data.shape} (should be [X, Y, 1, {n_samples}])")
  
  # Optional: Save the stacked data as NIfTI files
  t1_output_path = "/media/disk2/beijing_dti/processed/testingData/t1_testing_data"+str(slicei)+".nii.gz"
  fa_output_path = "/media/disk2/beijing_dti/processed/testingData/fa_testing_data"+str(slicei)+".nii.gz"
  
  t1_stacked_img = nib.Nifti1Image(t1_training_data, affine=np.eye(4))
  fa_stacked_img = nib.Nifti1Image(fa_training_data, affine=np.eye(4))
  
  nib.save(t1_stacked_img, t1_output_path)
  nib.save(fa_stacked_img, fa_output_path)
  
  print(f"T1 training data saved to {t1_output_path}")
  print(f"FA training data saved to {fa_output_path}")
