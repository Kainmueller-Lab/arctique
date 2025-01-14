import numpy as np
import shutil 
import os 
from pathlib import Path


full_dataset_folder = Path("C:/Users/cwinklm/Documents/Data/v_review_sample") #Path("/home/claudia/synthetic_uncertainty/data/v_review_full/render_review_full/")
target_folder =  Path("C:/Users/cwinklm/Documents/Data/v_review_sample10")

N_SAMPLES = 10
seed = 123

train_images = full_dataset_folder.joinpath("train/images")
test_images = full_dataset_folder.joinpath("test/images")
images_w_variations = full_dataset_folder.joinpath("variations/nuclei_intensity/0_00/images")


test_idx = [int(''.join(filter(str.isdigit, filename))) for filename in os.listdir(test_images)]
test_idx.sort()
train_idx = [int(''.join(filter(str.isdigit, filename))) for filename in os.listdir(train_images)]
train_idx.sort()
images_w_variations_idx = [int(''.join(filter(str.isdigit, filename))) for filename in os.listdir(images_w_variations)]
images_w_variations_idx.sort()


# check if all varied images are in test set 
print("IDX for varaitions", images_w_variations_idx)
var_in_test = True 
for idx in images_w_variations_idx: 
    if not idx in test_idx: 
        var_in_test = False
print("All varied images are in testset:", var_in_test)

np.random.seed(seed)
trainset_sample = np.random.choice(train_idx, N_SAMPLES, replace=False)
print(trainset_sample)

IDXs = {"train": trainset_sample, "test": images_w_variations_idx}

# copy variations: 
shutil.copytree(full_dataset_folder.joinpath("variations"), target_folder.joinpath("variations"))


# move subset
for mode, sample_idx in IDXs.items(): 
    print(mode, sample_idx)
    os.makedirs(target_folder.joinpath(f"{mode}/images"))
    for idx in sample_idx:
        shutil.copy(full_dataset_folder.joinpath(f"{mode}/images/img_{idx}.png"), target_folder.joinpath(f"{mode}/images/img_{idx}.png"))

    for mask_folder in os.listdir(full_dataset_folder.joinpath(f"{mode}/masks/")): 
        os.makedirs(target_folder.joinpath(f"{mode}/masks/{mask_folder}"))

        if mask_folder == "instance_3d_indexing": 
            for idx in sample_idx:
                shutil.copytree(full_dataset_folder.joinpath(f"{mode}/masks/{mask_folder}/{idx}"), target_folder.joinpath(f"{mode}/masks/{mask_folder}/{idx}"))
        else: 
            for idx in sample_idx:
                shutil.copy(full_dataset_folder.joinpath(f"{mode}/masks/{mask_folder}/{idx}.tif"), target_folder.joinpath(f"{mode}/masks/{mask_folder}/{idx}.tif"))