import yaml
import os
from pathlib import Path
import numpy as np
from sanity_check_functions import sanity_check_single_masks, sanity_check_semantic, sanity_check_FG_BG, sanity_check_instance_3class
from mask_helpers import make_val_split, read_val_idx_from_folder, copy_dataset_to_target_folders,  make_noisefree_variations, make_noisy_train_val_semantic, make_noisy_train_val_single, make_noisy_train_val_FG_BG, make_noisy_train_val_instance_3class, make_noisefree_test_set




def main(): 
    config_file_name = "/home/claudia/synthetic_uncertainty/src/data_prep/noisy_data_config.yaml"
    with open(config_file_name, 'r') as file:
        noisy_data_config = yaml.safe_load(file)

    DATA_PATH = Path(noisy_data_config["data_source"])
    print(DATA_PATH)

    TARGET_PATH = Path(noisy_data_config["target_location"])
    print(TARGET_PATH)
    if not os.path.isdir(TARGET_PATH): os.makedirs(TARGET_PATH)

    OVERWRITE = noisy_data_config["overwrite_exitsting"]
    SANITY_CHECK = noisy_data_config["sanity_check"]["execute"]
    FILE_EXTENSION = noisy_data_config["file_extension"]

    if noisy_data_config["create_new_train_val_set"]:
        if noisy_data_config["val_split"]["make_val_split"]: 
            val_split_config_dict = noisy_data_config["val_split"]
            original_image_folder = DATA_PATH.joinpath("train/images")
            val_idx = make_val_split(original_image_folder, val_split_config_dict)
            print("val idx:", val_idx)
        else:
            val_idx=[]
        copy_dataset_to_target_folders(DATA_PATH, TARGET_PATH, val_idx = val_idx, overwrite=OVERWRITE, file_extension=FILE_EXTENSION)
    else: 
        val_idx = read_val_idx_from_folder(TARGET_PATH)

    if SANITY_CHECK: 
        os.makedirs(TARGET_PATH.joinpath("sanity_check"), exist_ok=True)
        sanity_check_idx = [int(name.strip("."+'tif')) for name in os.listdir(TARGET_PATH.joinpath("train/masks/semantic_indexing"))]
        sanity_check_idx = np.random.choice(sanity_check_idx, size=noisy_data_config["sanity_check"]["samples"], replace=False)
    
    if noisy_data_config["make_noisefree_testset"]: 
        make_noisefree_test_set(DATA_PATH, TARGET_PATH, overwrite=OVERWRITE)

    if noisy_data_config["make_single_masks"]:
        make_noisy_train_val_single(TARGET_PATH, val_idx=val_idx, overwrite=OVERWRITE, file_extension=FILE_EXTENSION)
        if SANITY_CHECK:  sanity_check_single_masks(TARGET_PATH.joinpath("train/masks"), TARGET_PATH.joinpath("sanity_check"), sanity_check_idx)

    noisy_semantic_params = noisy_data_config["semantic"]
    if noisy_semantic_params["make_noisy_labels"]: 
        make_noisy_train_val_semantic(TARGET_PATH, noisy_semantic_params, val_idx=val_idx, overwrite=OVERWRITE, file_extension=FILE_EXTENSION)
        if SANITY_CHECK: sanity_check_semantic(TARGET_PATH.joinpath("train/masks"), TARGET_PATH.joinpath("sanity_check"), sanity_check_idx)

    noisy_FG_BG_params = noisy_data_config["FG_BG"]
    if noisy_FG_BG_params["make_noisy_labels"]:
        make_noisy_train_val_FG_BG(TARGET_PATH, noisy_FG_BG_params, val_idx=val_idx, overwrite=OVERWRITE, file_extension=FILE_EXTENSION)
        if SANITY_CHECK: sanity_check_FG_BG(TARGET_PATH.joinpath("train/masks"), TARGET_PATH.joinpath("sanity_check"), sanity_check_idx)
    
    noisy_instance_3class_params =  noisy_data_config["instance_3class"]
    if noisy_instance_3class_params["make_noisy_labels"]:
        make_noisy_train_val_instance_3class(TARGET_PATH, noisy_instance_3class_params, val_idx=val_idx, overwrite=OVERWRITE, file_extension=FILE_EXTENSION)
        if SANITY_CHECK: sanity_check_instance_3class(TARGET_PATH.joinpath("train/masks"), TARGET_PATH.joinpath("sanity_check"), sanity_check_idx)


    if noisy_data_config["variations"]["process_noisefree"]: 
        make_noisefree_variations(DATA_PATH, TARGET_PATH, overwrite=OVERWRITE, file_extension=FILE_EXTENSION)
    if noisy_data_config["variations"]["make_noisy_labels"]:
        raise NotImplementedError("Not yet implemented")

if __name__ == "__main__": 
    main()








