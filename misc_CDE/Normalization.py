#!/usr/bin/env python
# coding: utf-8

# In[2]:


from PIL import Image
import numpy as np
import os
import matplotlib.pyplot as plt

# In[11]:


# Normalize and save .tif image with compression
def normalize_image_with_compression(image_path, output_path, compression='tiff_lzw'):
    
    img = Image.open(image_path)
    img_array = np.array(img, dtype=np.float32)

    # Calculate the mean intensity
    mean_intensity = np.mean(img_array[img_array>0])
    # Apply the normalization: log10(1 + raw_intensity / image_mean)
    normalized_array = np.log10(1 + (img_array / mean_intensity))

    # Save with LZW compression (or another compression method)
    normalized_img = Image.fromarray(normalized_array)
    normalized_img.save(output_path, compression=compression)

Image.MAX_IMAGE_PIXELS = None


# In[40]:


# Function to process all images in the folders listed in the .txt file
def process_images_in_folders(txt_file_path, output_base_folder):
    # Step 1: Read the .txt file with the folder names
    with open(txt_file_path, 'r') as f:
        folder_paths = [line.strip() for line in f.readlines()]

    # Step 2: Loop through each folder and each image file
    for folder_path in folder_paths:
        folder_path0 = folder_path[1:]
        folder_path = "/media/disk2/GCA" + folder_path0
        # Ensure the folder exists
        if os.path.isdir(folder_path):
            # Get all the image files in the folder
            for file_name in os.listdir(folder_path):
                if file_name.lower().endswith('.tif'):  # Process only .tif/.tiff files
                    input_image_path = os.path.join(folder_path, file_name)
                    
                    # Create an output path (in the same folder or different)
                    output_folder = output_base_folder + folder_path0
                    os.makedirs(output_folder, exist_ok=True)  # Create output folder if it doesn't exist
                    
                    output_image_path = os.path.join(output_folder, file_name)
                    
                    # Apply the normalization function
                    normalize_image_with_compression(input_image_path, output_image_path)
                    print(f"Processed {file_name} in {folder_path}")
        else:
            print(f"Folder not found: {folder_path}")


# In[41]:


process_images_in_folders('/media/disk2/GCA/create_dir.txt', "/media/disk2/GCA_norm")

