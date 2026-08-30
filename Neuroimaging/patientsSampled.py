import pandas as pd
import nibabel as nib
import numpy as np
import random
import os

# File paths
csv_path = "/media/disk2/beijing_dti/enhanced/unpack/BeijingEnhancedDTIProcessed1.csv"

# Read the subject IDs from the CSV file and ensure they are strings with leading zeros
subject_data = pd.read_csv(csv_path)
subject_data['NewID'] = subject_data['NewID'].astype(str).str.zfill(3)
# Add a column with the original row numbers
print(subject_data.index)
subject_data['original_row'] = subject_data.index

# Randomly sample 89 subjects
n_samples = 89
sampled_subject_data = subject_data.sample(n=n_samples, random_state=42)  # Set random_state for reproducibility
# Get the rows that were NOT sampled
test_subject_data = subject_data.drop(sampled_subject_data.index)

# Save the sampled data to a CSV file
sampled_subject_data.to_csv('/media/disk2/beijing_dti/Neuroimaging/train_subject_data.csv', index=False)

# Save the non-sampled data to another CSV file
test_subject_data.to_csv('/media/disk2/beijing_dti/Neuroimaging/testt_subject_data.csv', index=False)
