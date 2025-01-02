import CycleGAN
from CycleGAN import *

# Create a CycleGAN on GPU 0
myCycleGAN = CycleGAN(0)

import os


for slicez in range(72, 182):
        
    input_dir="/home/xiongj3/ImageImputation/Neuroimaging/train_T1_FA/trainingSlices"
    trainA_dir = input_dir+"/t1_training_data_"+str(slicez)+".nii.gz"
    trainB_dir = input_dir+"/fa_training_data_"+str(slicez)+".nii.gz"
    models_dir = "/home/xiongj3/ImageImputation/Neuroimaging/train_T1_FA/models/slice"+str(slicez)
    try:
        os.mkdir(models_dir)
        print(f"Directory '{models_dir}' created successfully.")
    except FileExistsError:
        print(f"Directory '{models_dir}' already exists.")
    
    output_sample_dir = models_dir+'/outputsample.png'
    batch_size = 10
    epochs = 200
    normalization_factor_A = 1000
    normalization_factor_B = 1
    import time

    # Start the timer
    start_time = time.time()
    myCycleGAN.train(trainA_dir, normalization_factor_A, trainB_dir, normalization_factor_B, models_dir, batch_size, epochs, output_sample_dir=output_sample_dir, output_sample_channels=1)
    end_time = time.time()

    # Compute the elapsed time
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")
