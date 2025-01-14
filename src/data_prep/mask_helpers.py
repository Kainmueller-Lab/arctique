import numpy as np
from pathlib import Path
import os 
from PIL import Image
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib
import shutil
import torch
from scipy import ndimage
import time



#######################################################
#  TODO
# double-check creation of instance masks
# add comments to functions 
#######################################################


def setup_folder(location, overwrite, makedir=False): 
    '''
    Function to setup a folder. If the folder already exists, it can be overwritten.
    '''

    if os.path.isdir(location): 
        if overwrite==False: 
            raise ValueError(f"Cannot create {location} because it already exists. Specify a different folder or set 'overwrite=True'")

        else: 
             print(f"Overwritig esisting folder {location}")
             shutil.rmtree(location)
    
    if makedir: 
        os.makedirs(location)


def make_val_split(train_image_folder, val_split_config_dict):
    '''
    Function to create a validation split.
    '''
    if val_split_config_dict["make_val_split"] == True: 
        
        val_frac = val_split_config_dict["val_frac"]

        all_filenames = os.listdir(train_image_folder)
        all_idx = [int(''.join(filter(str.isdigit, filename))) for filename in all_filenames]
        
        np.random.seed(val_split_config_dict["val_frac_seed"])
        val_idx = np.random.choice(all_idx, size=int(val_frac*len(all_idx)), replace=False)

        return val_idx
    else: 
        return []
    

def read_val_idx_from_folder(target_folder):
    '''
    Function to read the validation indices from a folder. Used when the target folder already contains a validation split.
    '''
    if not os.path.isdir(target_folder.joinpath("val")): 
        return []
    else: 
        val_idx = [int(name.strip(".png").strip("img_")) for name in os.listdir(target_folder.joinpath("val/images"))]
        return val_idx


def copy_dataset_to_target_folders(source_folder, target_folder, val_idx = [], overwrite=False, file_extension="png"): 
    '''
    Function to copy the dataset to the target folder. All images and masks will be converted to the specified file extension.
    '''

    source_folder = Path(source_folder)
    target_folder = Path(target_folder)
    os.makedirs(target_folder, exist_ok=True)
    
    # copy test folder as-is
    print(f"copy test set to target location: {target_folder}")
    setup_folder(target_folder.joinpath("test"), overwrite=overwrite)
    shutil.copytree(source_folder.joinpath("test"), target_folder.joinpath("test"))


    print(f"copy train set to target location: {target_folder}")
    if len(val_idx) >0: 

        if os.path.isdir(target_folder.joinpath("val/images")): 
            val_idx.sort()  
            all_filenames = os.listdir(target_folder.joinpath("val/images"))
            all_idx = [int(''.join(filter(str.isdigit, filename))) for filename in all_filenames]
            all_idx.sort()

            if overwrite == True: 
                if list(val_idx) == list(all_idx): 
                    print("removing exitsing image folders")
                    shutil.rmtree(target_folder.joinpath("val/images"))
                    shutil.rmtree(target_folder.joinpath("train/images"))
                else: 
                    raise ValueError(f"Cannot overwrite {target_folder.joinpath('val/images')} because it contains different indices than currently specified.")

        os.makedirs(target_folder.joinpath("train/images"), exist_ok=True)
        os.makedirs(target_folder.joinpath("val/images"), exist_ok=True)

        subdirs = os.listdir(source_folder.joinpath("train/masks"))

        for subdir in subdirs: 

            setup_folder(target_folder.joinpath(f"train/masks/{subdir}"), overwrite=overwrite, makedir=True)
            if len(val_idx) >0: 
                setup_folder(target_folder.joinpath(f"val/masks/{subdir}"), overwrite=overwrite, makedir=True)


        all_filenames = os.listdir(source_folder.joinpath("train/images"))
        all_idx = [int(''.join(filter(str.isdigit, filename))) for filename in all_filenames]

        for idx in all_idx: 
            location = "val" if idx in val_idx else "train"
            shutil.copy(source_folder.joinpath(f"train/images/img_{idx}.{'png'}"), target_folder.joinpath(f"{location}/images/img_{idx}.{file_extension}"))
            for subdir in subdirs:
                if subdir in ['instance_3d_indexing']: 
                    # for mask_slice_file in os.listdir(source_folder.joinpath(f"train/masks/{subdir}/{idx}")):
                    #     print(mask_slice_file)
                    #     save_name = mask_slice_file.split('.')
                    #     print(save_name) ### 
                    #     shutil.copy(source_folder.joinpath(f"train/masks/{subdir}/{idx}/{mask_slice_file}"), target_folder.joinpath(f"{location}/masks/{subdir}/{idx}/{save_name[0]}.{'png'}"))
                    shutil.copytree(source_folder.joinpath(f"train/masks/{subdir}/{idx}"), target_folder.joinpath(f"{location}/masks/{subdir}/{idx}"))
                else:
                    shutil.copy(source_folder.joinpath(f"train/masks/{subdir}/{idx}.{'tif'}"), target_folder.joinpath(f"{location}/masks/{subdir}/{idx}.{file_extension}"))
    else: 
        setup_folder(target_folder.joinpath("train"), overwrite=overwrite)
        shutil.copytree(source_folder.joinpath("train"), target_folder.joinpath("train"))


def make_noisefree_test_set(data_location, target_folder, overwrite=False, file_extension="png"): 

    noiseless_FG_BG_test_mask_path = target_folder.joinpath(f"test/masks/FG_BG_noise_0")
    setup_folder(noiseless_FG_BG_test_mask_path, overwrite=overwrite)
    shutil.copytree(data_location.joinpath(f"test/masks/semantic_indexing"), noiseless_FG_BG_test_mask_path)

    for mask_name in os.listdir(noiseless_FG_BG_test_mask_path):
        mask = np.array(Image.open(noiseless_FG_BG_test_mask_path.joinpath(mask_name)))
        mask[mask>0] = 1
        mask_png = Image.fromarray(mask.astype(np.uint8))
        mask_label, _ = mask_name.split(".")
        mask_png.save(noiseless_FG_BG_test_mask_path.joinpath(f"{mask_label}.{file_extension}"))
        os.remove(noiseless_FG_BG_test_mask_path.joinpath(f"{mask_label}.tif"))


    noiseless_semantic_test_mask_path = target_folder.joinpath(f"test/masks/semantic_noise_0")
    setup_folder(noiseless_semantic_test_mask_path, overwrite=overwrite)
    shutil.copytree(data_location.joinpath(f"test/masks/semantic_indexing"), noiseless_semantic_test_mask_path)

    for mask_name in os.listdir(noiseless_semantic_test_mask_path):
        mask = np.array(Image.open(noiseless_semantic_test_mask_path.joinpath(mask_name)))
        mask_png = Image.fromarray(mask)
        mask_label, _ = mask_name.split(".")
        mask_png.save(noiseless_semantic_test_mask_path.joinpath(f"{mask_label}.{file_extension}"))
        os.remove(noiseless_semantic_test_mask_path.joinpath(f"{mask_label}.tif"))


    noiseless_instance_3class_test_mask_path = target_folder.joinpath(f"test/masks/instance_3class_noise_0")
    setup_folder(noiseless_instance_3class_test_mask_path, overwrite=overwrite)

    shutil.copytree(data_location.joinpath(f"test/masks/instance_indexing"), noiseless_instance_3class_test_mask_path)
    for mask_name in os.listdir(noiseless_instance_3class_test_mask_path):
        instance_mask = np.array(Image.open(noiseless_instance_3class_test_mask_path.joinpath(mask_name)))
        instance_3class_mask = build_instance_3class_mask(instance_mask)
        instance_3class_mask = Image.fromarray(instance_3class_mask.astype(np.uint8))
        mask_label, _ = mask_name.split(".")
        instance_3class_mask.save(noiseless_instance_3class_test_mask_path.joinpath(f"{mask_label}.{file_extension}"))
        os.remove(noiseless_instance_3class_test_mask_path.joinpath(f"{mask_label}.tif"))


    noiseless_instance_3class_test_mask_path = target_folder.joinpath(f"test/masks/large_cell_borders")
    setup_folder(noiseless_instance_3class_test_mask_path, overwrite=overwrite)

    shutil.copytree(data_location.joinpath(f"test/masks/instance_indexing"), noiseless_instance_3class_test_mask_path)
    for mask_name in os.listdir(noiseless_instance_3class_test_mask_path):
        instance_mask = np.array(Image.open(noiseless_instance_3class_test_mask_path.joinpath(mask_name)))
        instance_3class_mask = build_instance_3class_mask(instance_mask, border_width=10)
        instance_3class_mask[instance_3class_mask<2]=0
        instance_3class_mask[instance_3class_mask==2]=1
        instance_3class_mask = Image.fromarray(instance_3class_mask.astype(np.uint8))
        mask_label, _ = mask_name.split(".")
        instance_3class_mask.save(noiseless_instance_3class_test_mask_path.joinpath(f"{mask_label}.{file_extension}"))
        os.remove(noiseless_instance_3class_test_mask_path.joinpath(f"{mask_label}.tif"))


    noiseless_instance_3class_test_mask_path = target_folder.joinpath(f"test/masks/large_cell_borders_epithel")
    setup_folder(noiseless_instance_3class_test_mask_path, overwrite=overwrite)

    shutil.copytree(data_location.joinpath(f"test/masks/instance_indexing"), noiseless_instance_3class_test_mask_path)
    for mask_name in os.listdir(noiseless_instance_3class_test_mask_path):

        mask_label, _ = mask_name.split(".")
        semantic_mask = np.array(Image.open(noiseless_semantic_test_mask_path.joinpath(f"{mask_label}.{file_extension}")))
        semantic_mask[semantic_mask !=1] = 0 

        instance_mask = np.array(Image.open(noiseless_instance_3class_test_mask_path.joinpath(mask_name)))
        instance_mask = np.multiply(semantic_mask, instance_mask)

        instance_3class_mask = build_instance_3class_mask(instance_mask, border_width=10)
        instance_3class_mask[instance_3class_mask<2]=0
        instance_3class_mask[instance_3class_mask==2]=1
        instance_3class_mask = Image.fromarray(instance_3class_mask.astype(np.uint8))
        #mask_label, _ = mask_name.split(".")
        instance_3class_mask.save(noiseless_instance_3class_test_mask_path.joinpath(f"{mask_label}.{file_extension}"))
        os.remove(noiseless_instance_3class_test_mask_path.joinpath(f"{mask_label}.tif"))


def make_noisefree_variations(data_location, target_folder, overwrite=True, file_extension="png"):
    var_location  = data_location.joinpath("variations")
    setup_folder(target_folder.joinpath("variations"), overwrite=overwrite, makedir=True)
    print("setup", target_folder.joinpath("variations"))

    for variation in os.listdir(var_location):
        for noise_level in os.listdir(var_location.joinpath(variation)):
            print(variation, noise_level)

            variation_image_path = target_folder.joinpath(f"variations/{variation}/{noise_level}/images")
            setup_folder(variation_image_path, overwrite=overwrite, makedir=True)
            print("setup", variation_image_path)
            IDXSs = [int(name.strip("img_").strip(".png")) for name in os.listdir(var_location.joinpath(f"{variation}/{noise_level}/images"))]
            for idx in IDXSs:
                shutil.copy(var_location.joinpath(f"{variation}/{noise_level}/images/img_{idx}.{'png'}"), variation_image_path.joinpath(f"img_{idx}.{file_extension}"))


            noiseless_FG_BG_variation_mask_path = target_folder.joinpath(f"variations/{variation}/{noise_level}/masks/FG_BG_noise_0")
            setup_folder(noiseless_FG_BG_variation_mask_path, overwrite=overwrite)
            shutil.copytree(var_location.joinpath(f"{variation}/{noise_level}/masks/semantic_indexing"), noiseless_FG_BG_variation_mask_path)

            for mask_name in os.listdir(noiseless_FG_BG_variation_mask_path):
                mask = np.array(Image.open(noiseless_FG_BG_variation_mask_path.joinpath(mask_name)))
                mask[mask>0] = 1
                mask_png = Image.fromarray(mask.astype(np.uint8))
                mask_label, _ = mask_name.split(".")
                mask_png.save(noiseless_FG_BG_variation_mask_path.joinpath(f"{mask_label}.{file_extension}"))
                os.remove(noiseless_FG_BG_variation_mask_path.joinpath(f"{mask_label}.tif"))


            noiseless_semantic_variation_mask_path = target_folder.joinpath(f"variations/{variation}/{noise_level}/masks/semantic_noise_0")
            setup_folder(noiseless_semantic_variation_mask_path, overwrite=overwrite)
            shutil.copytree(var_location.joinpath(f"{variation}/{noise_level}/masks/semantic_indexing"), noiseless_semantic_variation_mask_path)

            for mask_name in os.listdir(noiseless_semantic_variation_mask_path):
                mask = np.array(Image.open(noiseless_semantic_variation_mask_path.joinpath(mask_name)))
                mask_png = Image.fromarray(mask)
                mask_label, _ = mask_name.split(".")
                mask_png.save(noiseless_semantic_variation_mask_path.joinpath(f"{mask_label}.{file_extension}"))
                os.remove(noiseless_semantic_variation_mask_path.joinpath(f"{mask_label}.tif"))


            noiseless_instance_3class_variation_mask_path = target_folder.joinpath(f"variations/{variation}/{noise_level}/masks/instance_3class_noise_0")
            setup_folder(noiseless_instance_3class_variation_mask_path, overwrite=overwrite)

            shutil.copytree(data_location.joinpath(f"test/masks/instance_indexing"), noiseless_instance_3class_variation_mask_path)
            for mask_name in os.listdir(noiseless_instance_3class_variation_mask_path):
                instance_mask = np.array(Image.open(noiseless_instance_3class_variation_mask_path.joinpath(mask_name)))
                instance_3class_mask = build_instance_3class_mask(instance_mask)
                instance_3class_mask = Image.fromarray(instance_3class_mask.astype(np.uint8))
                mask_label, _ = mask_name.split(".")
                instance_3class_mask.save(noiseless_instance_3class_variation_mask_path.joinpath(f"{mask_label}.{file_extension}"))
                os.remove(noiseless_instance_3class_variation_mask_path.joinpath(f"{mask_label}.tif"))




def make_noisy_train_val_semantic(data_location, params, val_idx=[], overwrite=False, file_extension="png"): 

    print("generating noiseless semantic masks with corrected labels")

    for split in ["train", "val"]: 
        if split=="val" and len(val_idx)==0: 
            continue

        noiseless_semantic_mask_path = data_location.joinpath(f"{split}/masks/semantic_noise_0")
        setup_folder(noiseless_semantic_mask_path, overwrite=overwrite)

        print("generating noise free dataset")
        start = time.time()

        shutil.copytree(data_location.joinpath(f"{split}/masks/semantic_indexing"), noiseless_semantic_mask_path)

        for mask_name in os.listdir(noiseless_semantic_mask_path):
            mask = np.array(Image.open(noiseless_semantic_mask_path.joinpath(mask_name)))
            mask_png = Image.fromarray(mask)
            mask_png.save(noiseless_semantic_mask_path.joinpath(f"{mask_name.split('.')[0]}.{file_extension}"))

        stop = time.time()
        print(f"this took {np.round(stop-start, 4)} seconds")


    for nl in params["noise_levels"]: 

        noise_params = params["noise_params"]
        noise_params["global_flip_prob"] = nl/100

        print(f"generating noise level {nl}")
        start = time.time()
        for split in ["train", "val"]: 
            if split=="val" and len(val_idx)==0: 
                continue
    
            source_folder = data_location.joinpath(f"{split}/masks/semantic_noise_0")
            target_folder = data_location.joinpath(f"{split}/masks/semantic_noise_{nl}")
            setup_folder(target_folder, overwrite=overwrite, makedir=True)


            all_mask_names = os.listdir(source_folder)
            all_mask_idx = [int(''.join(filter(str.isdigit, filename))) for filename in all_mask_names]

            for mask_idx in all_mask_idx: 
                noisy_mask_array = build_noisy_semantic_mask(source_folder, 
                                        mask_idx, 
                                        clean_semantic_mask_folder = data_location.joinpath(f"{split}/masks/semantic_noise_0"),
                                        clean_instance_mask_folder = data_location.joinpath(f"{split}/masks/instance_indexing"),
                                        **params["noise_params"])
                
                noisy_mask = Image.fromarray(noisy_mask_array)
                noisy_mask.save(target_folder.joinpath(f"{mask_idx}.{file_extension}"))
        stop = time.time()
        print(f"this took {np.round(stop-start, 4)} seconds")


def make_noisy_train_val_instance_3class(data_location, params, val_idx=[], overwrite=False, file_extension="png"):
    neighbor_dict = {}

    for split in ["train", "val"]: 

        if split=="val" and len(val_idx)==0: 
            continue
        print("generating noiseless insta masks")        
        noiseless_instance_3class_mask_path = data_location.joinpath(f"{split}/masks/instance_3class_noise_0")
        setup_folder(noiseless_instance_3class_mask_path, overwrite)

        shutil.copytree(data_location.joinpath(f"{split}/masks/instance_indexing"), noiseless_instance_3class_mask_path)
        for mask_name in os.listdir(noiseless_instance_3class_mask_path):
            instance_mask = np.array(Image.open(noiseless_instance_3class_mask_path.joinpath(mask_name)))
            instance_3class_mask = build_instance_3class_mask(instance_mask, border_width=params["border_width"])
            instance_3class_mask = Image.fromarray(instance_3class_mask.astype(np.uint8))
            instance_3class_mask.save(noiseless_instance_3class_mask_path.joinpath(f"{mask_name.split('.')[0]}.{file_extension}"))
            
    for nl_idx, nl in enumerate(params["noise_levels"]): 
        start = time.time()
        
        print(f"generating noise level {nl}: neighboring cells are merged with probability: {nl/100}")

        for split in ["train", "val"]: 
            if split=="val" and len(val_idx)==0:  
                continue

            instance_label_mask_folder = data_location.joinpath(f"{split}/masks/instance_indexing")
            target_folder = data_location.joinpath(f"{split}/masks/instance_3class_noise_{nl}")
            setup_folder(target_folder, overwrite=overwrite, makedir=True)

            all_mask_names = os.listdir(instance_label_mask_folder)
            all_mask_idx = [int(''.join(filter(str.isdigit, filename))) for filename in all_mask_names]

          
            for idx, mask_idx in enumerate(all_mask_idx): 
                print(f"{(idx+1):03d}/{len(all_mask_idx)}", end="\r")
                if nl_idx ==0:
                    instance_mask = np.array(Image.open(instance_label_mask_folder.joinpath(f"{mask_idx}.{file_extension}")), dtype=int)
                    nieghboring_cells, _ = get_list_of_neighboring_cells(instance_mask)
                    neighbor_dict[mask_idx] = nieghboring_cells


                #merged_instance_3class_mask = build_merged_instance_mask(single_mask_folder.joinpath(f"{mask_idx}"), three_class=False, merge_percent=(nl/100), border_width=params["border_width"])
                merged_instance_3class_mask = build_merged_instance_mask(instance_label_mask_folder.joinpath(f"{mask_idx}.{file_extension}"), three_class=True, merge_percent=(nl/100), 
                                                                         border_width=params["border_width"], neighbors=neighbor_dict[mask_idx])

                noisy_mask = Image.fromarray(merged_instance_3class_mask.astype(np.uint8))
                noisy_mask.save(target_folder.joinpath(f"{mask_idx}.{file_extension}"))
        stop = time.time()
        print(f"this took {np.round(stop-start, 4)} seconds")


def make_noisy_train_val_single(data_location, val_idx=[], overwrite=False, file_extension=".png"):
    for split in ["train", "val"]: 
        if split=="val" and len(val_idx)==0: 
            continue

        start = time.time()
        print(f"creating individual cells for {split}")
        mask_folder_3d = data_location.joinpath(f'{split}/masks/instance_3d_indexing/')
        target_folder = data_location.joinpath(f'{split}/masks/single_masks/')

        setup_folder(target_folder, overwrite=overwrite)


        make_individual_masks(mask_folder_3d)
        stop = time.time()
        print(f"this took {np.round(stop-start, 4)} seconds")




def make_noisy_train_val_FG_BG(data_location, params, val_idx=[], overwrite=False, file_extension="png"):

    for split in ["train", "val"]: 

        if split=="val" and len(val_idx)==0: 
            continue
        print("generating noiseless FG-BG masks")        
        noiseless_FG_BG_mask_path = data_location.joinpath(f"{split}/masks/FG_BG_noise_0")
        setup_folder(noiseless_FG_BG_mask_path, overwrite=overwrite)

        shutil.copytree(data_location.joinpath(f"{split}/masks/semantic_indexing"), noiseless_FG_BG_mask_path)
        for mask_name in os.listdir(noiseless_FG_BG_mask_path):
            mask = np.array(Image.open(noiseless_FG_BG_mask_path.joinpath(mask_name)))
            mask[mask>0] = 1
            mask_png = Image.fromarray(mask.astype(np.uint8))
            mask_png.save(noiseless_FG_BG_mask_path.joinpath(f"{mask_name.split('.')[0]}.{file_extension}"))

    noise_params = params["noise_params"]
    noise_params_probs = [key for key in noise_params if str(key).startswith("prob")]

    for nl_idx, nl in enumerate(params["noise_levels"]): 
        start = time.time()
        
        mod_prob = nl/(100*4) # divide by 100 to make probability, divide by for tow split prob equally across possible modifications 
        print(f"generating noise level {nl}: each cell is modified with probability: {nl/100}")

        for npp_idx, npp in enumerate(noise_params_probs): 
            noise_params[npp] = min(mod_prob, 1)

        start = time.time()
        for split in ["train", "val"]: 
            if split=="val" and len(val_idx)==0:  
                continue

            source_folder = data_location.joinpath(f"{split}/masks/FG_BG_noise_0")
            single_mask_folder = data_location.joinpath(f"{split}/masks/single_masks")
            if not os.path.isdir(single_mask_folder): 
                 print("creating single masks")
                 make_noisy_train_val_single(data_location, val_idx=val_idx)
            target_folder = data_location.joinpath(f"{split}/masks/FG_BG_noise_{nl}")
            

            setup_folder(target_folder, overwrite=overwrite, makedir=True)

            all_mask_names = os.listdir(source_folder)
            all_mask_idx = [int(''.join(filter(str.isdigit, filename))) for filename in all_mask_names]

            for idx, mask_idx in enumerate(all_mask_idx): 
                print(f"{(idx+1):03d}/{len(all_mask_idx)}", end="\r")
                noisy_mask_array = build_noisy_FG_BG_mask(single_mask_folder.joinpath(f"{mask_idx}"), **noise_params)

                noisy_mask = Image.fromarray(noisy_mask_array.astype(np.uint8))
                noisy_mask.save(target_folder.joinpath(f"{mask_idx}.{file_extension}"))

        stop = time.time()
        print(f"this took {np.round(stop-start, 4)} seconds")




def build_noisy_semantic_mask(data_path, 
                              mask_idx, 
                              clean_semantic_mask_folder = "semantic_noise_0",
                              clean_instance_mask_folder = "instance_indexing",
                              global_flip_prob=0.1, # probability to flip any class label with any other class label. Overwritten by class_flip_dict
                              class_flip_dict = {}, # assigns each cell typ a probability to be flipped and which cell types it can be flipped with
                              file_extension="png"):
    
    '''
    Example class_flip_dict {1: (0.1, [2]), 3: (0.3, [2,4])}: 
        - class 1 is flipped to class 2 with probability 0.1
        - class 3 is flipped to class 2 or class 4 with probability 0.3
    '''

    true_semantic_mask = np.array(Image.open(Path(data_path).joinpath(f"{clean_semantic_mask_folder}/{mask_idx}.{file_extension}")))
    true_instance_mask = np.array(Image.open(Path(data_path).joinpath(f"{clean_instance_mask_folder}/{mask_idx}.{'tif'}")))
    
    noisy_semantic_mask = true_semantic_mask.copy()
    
    object_classes = list(np.unique(true_semantic_mask))# all classes present in semantic mask 
    object_classes.remove(0) # background cannot be flipped 

    if (class_flip_dict is not None): 
        if len(class_flip_dict)>0: 
            flippable_classes = list(class_flip_dict.keys()) # classes that can be flipped 
    else: 
        
        if global_flip_prob > 1: 
            global_flip_prob = global_flip_prob/100
        #print(f"Flipping class labels randomly with probability {global_flip_prob}")
        flippable_classes = object_classes # if no class_flip_dict is defined all classes can be flipped with all others
        class_flip_dict = {cl:(global_flip_prob, [f for f in flippable_classes if not f==cl]) for cl in flippable_classes}
        #print(class_flip_dict)


    for object_id in np.unique(true_instance_mask): # iterate through all visible cells 

        if object_id==0: # ignore background 
            continue

        single_object_mask = (true_instance_mask==object_id).astype(int) # pick out single cell             
        single_object_class_mask = np.multiply(true_semantic_mask, single_object_mask) # identify cells class 
        single_object_class = np.max(single_object_class_mask) 

        if single_object_class in flippable_classes: # if cell is flippable
            if  len(class_flip_dict[single_object_class][1])>0:
                if np.random.rand() < class_flip_dict[single_object_class][0]: # it is flipped with the probability defined in class_flip_dict
                    noisy_semantic_mask[single_object_mask==1] = np.random.choice(class_flip_dict[single_object_class][1]) # and assigned one of the potential alternative classes 


    return noisy_semantic_mask




def sef_setup_empty_mask(single_mask_array): 
    assert(len(np.unique(single_mask_array))<=2) # check if single mask has max. two values 
    base_mask  = np.zeros_like(single_mask_array)
    return base_mask


def get_instance(stack, id):
    instance_stack = (stack==id)
    return instance_stack.sum(dim=2).clip(0,1)


def make_individual_masks(mask_folder_3d, file_extension="png"): 

    mask_folder_3d = Path(mask_folder_3d)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    all_sample_idx = [int(idx) for idx in os.listdir(mask_folder_3d)]

    for j, sample_idx in enumerate(all_sample_idx): 
        print(f"{(j+1):03d}/{len(all_sample_idx)}", end="\r")
        stack = np.load(mask_folder_3d.joinpath(f"{sample_idx}/{sample_idx}.npy"))
        n_vals = len(np.unique(stack))

        stack = np.array(stack, dtype=np.int16)
        stack = torch.tensor(stack, dtype=torch.int16, device=device).permute(1,2,0)
        instances = [np.array(get_instance(stack, i).detach().cpu().numpy(), dtype=np.uint8) for i in range(1, n_vals)]

        target_path = mask_folder_3d.parent.joinpath(f"single_masks/{sample_idx}")
        os.makedirs(target_path, exist_ok=True)


        for instance_id, instance_mask in enumerate(instances): 
            mask_png = Image.fromarray(instance_mask)
            mask_png.save(target_path.joinpath(f"{instance_id}.{file_extension}"))




def shift_mask(mask, shift_x=10, shift_y=10): 
    """
    Shifts the input mask by specified values in x and y direction. 
    """
    shifted_image = ndimage.shift(mask, [shift_x, shift_y])

    return shifted_image.astype(int)

def scale_mask(image, increase=True, struct=np.ones((3,3)), iters=5): 
    """
    Scales the input mask by iteratively applying 'struct' to each pixel. 
    For further info see doczmentation of ndimage.binary_dilation; ndimage.binary_erosion
    Note max value has to be 1
    """
    
    #assert(np.max(image)<=1)
    image = image    
    if increase: 
        scaled_image = ndimage.binary_dilation(image, structure=struct, iterations=iters)
    else: 
        scaled_image = ndimage.binary_erosion(image, structure=struct, iterations=iters)
    return scaled_image.astype(int)


def elastic_transform(image, alpha, sigma, random_state=None):
    """Elastic deformation of images as described in [Simard2003]_.
    .. [Simard2003] Simard, Steinkraus and Platt, "Best Practices for
       Convolutional Neural Networks applied to Visual Document Analysis", in
       Proc. of the International Conference on Document Analysis and
       Recognition, 2003.

       Code copied from: https://gist.github.com/chsasank/4d8f68caf01f041a6453e67fb30f8f5a
    """
    assert len(image.shape)==2
    #assert(np.max(image)<=1)

    if random_state is None:
        random_state = np.random.RandomState(None)

    shape = image.shape

    dx = ndimage.filters.gaussian_filter((random_state.rand(*shape) * 2 - 1), sigma, mode="constant", cval=0) * alpha
    dy = ndimage.filters.gaussian_filter((random_state.rand(*shape) * 2 - 1), sigma, mode="constant", cval=0) * alpha

    x, y = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]), indexing='ij')
    indices = np.reshape(x+dx, (-1, 1)), np.reshape(y+dy, (-1, 1))

    transformed_image = ndimage.interpolation.map_coordinates(image, indices, order=1).reshape(shape)
    transformed_image[transformed_image>0.5]  = 1
    transformed_image[transformed_image<=0.5] = 0
    
    return transformed_image.astype(int)



def build_noisy_FG_BG_mask(single_mask_path, max_pixel_val=255.,
                    prob_shift=0.1, max_shift_strength=10, 
                    prob_scale=0.1, max_scale_strength=7, 
                    prob_elastic=0.1, min_alpha_elastic=5,  max_alpha_elastic=10, min_sigma_elastic=0.5, max_sigma_elastic=2,   
                    prob_leave_out=0.05, verbose=False): 
    """
    Combines multiple masks of individual cell objects to a semantic mask (differentiating only foreground and backgorund). 
    Depending on the probabilities prob_shift and prob_scale any individual cell-mask can be shifted and/or scaled before being 
    added to the full mask. 
    """
    single_mask_files = os.listdir(single_mask_path)
    scaling_factor = max_pixel_val

    for i in range(len(single_mask_files)): 
        single_mask_file = Path(single_mask_path).joinpath(single_mask_files[i]) 
        single_mask = np.array(Image.open(single_mask_file)).astype("float")
        if i == 0:
            summed_FG_BG_mask = sef_setup_empty_mask(single_mask)

        if np.max(single_mask) > 1: 
            single_mask = (single_mask/scaling_factor).astype(int)

        if np.max(single_mask)==0: 
            continue

        x_shift, x_scale, x_elastic, x_leave_out = np.random.rand(4)
        if verbose:
            print(f"shift:{x_shift < prob_shift}, scale: {x_scale < prob_scale}, elastic: {x_elastic < prob_elastic}, leave out:{ x_leave_out < prob_leave_out}")
            
        if x_shift < prob_shift: 
            shift_x, shift_y = np.random.randint(-max_shift_strength, max_shift_strength, size=2)
            if verbose: 
                print("shift values:", shift_x, shift_y)
            single_mask = shift_mask(single_mask, shift_x=shift_x, shift_y=shift_y)

        if x_scale < prob_scale: 
            scale_direction = np.random.rand() < 0.5 # increase or decrease 
            iters = np.random.randint(1, max_scale_strength) # increase/decrease factor
            single_mask = scale_mask(single_mask, increase=scale_direction, iters=iters)
            
        if x_elastic < prob_elastic: 
            alpha = min_alpha_elastic + (max_alpha_elastic-min_alpha_elastic)*np.random.rand() 
            sigma = min_sigma_elastic + (max_sigma_elastic-min_sigma_elastic)*np.random.rand()#min(alpha*0.1, max_sigma_elastic)
            if verbose:
                print(f"alpha = {np.round(alpha, 2)}, sigma = {np.round(sigma, 2)}")
            single_mask = elastic_transform(single_mask, alpha=alpha, sigma=sigma)
            
        if x_leave_out < prob_leave_out: 
            single_mask = np.zeros_like(single_mask)

        if verbose: 
            vals, counts = np.unique(single_mask, return_counts=True)
            print(vals, counts)
            if (len(vals) == 2) and (counts[1]>counts[0]): 
                break 
        summed_FG_BG_mask += single_mask

    return summed_FG_BG_mask.clip(0,1).astype(np.uint8)



# def get_list_of_neighboring_cells(single_cell_folder): 

#     single_mask_names = os.listdir(single_cell_folder)
#     #print(single_mask_names)

#     Neighbors = []
#     for mask_id, mask_name in enumerate(single_mask_names): 
#         single_mask = (np.array(Image.open(single_cell_folder.joinpath(mask_name)))).astype(int)

#         if mask_id == 0: 
#             FG_BG = np.zeros_like(single_mask)
#             Instance_mask = np.zeros_like(single_mask) # this should be the true instance mask 

#         FG_BG += single_mask
        
#         if 2 in np.unique(FG_BG): 
#             neighbors = [mask_id+1] # this should be the true mask id read from filename 
#             overlap = FG_BG==2

#             overlapping_instances = list(np.unique(Instance_mask[overlap])) 
#             #overlapping_instances = [o  for o in overlapping_instances if o<=mask_id]
#             neighbors.extend(list(overlapping_instances))
#             Neighbors.append(tuple(neighbors))

#         FG_BG = FG_BG.clip(0, 1)
#         Instance_mask += single_mask*(mask_id+1)
#         Instance_mask = Instance_mask.clip(0, mask_id+1)


#     return Neighbors, Instance_mask



def get_list_of_neighboring_cells(label_instance_mask_orig): 
    
    label_instance_mask_orig = np.array(label_instance_mask_orig, dtype = int)
    label_instance_mask = label_instance_mask_orig.copy()

    label_instance_mask[label_instance_mask==0] = 10000

    label_instance_mask_diff_x = np.logical_and(np.diff(label_instance_mask, axis=0, append=0)!=0, abs(np.diff(label_instance_mask, axis=0, append=0))<9000)
    label_instance_mask_diff_y = np.logical_and(np.diff(label_instance_mask, axis=1, append=0)!=0, abs(np.diff(label_instance_mask, axis=1, append=0))<9000)

    label_instance_mask_diff = label_instance_mask_diff_x+label_instance_mask_diff_y

    x, y = np.where(label_instance_mask_diff!=0)

    NEIGBORS = []
    for i in range(len(x)):     
        neighs = np.unique(label_instance_mask_orig[x[i]-1:x[i]+1, y[i]-1:y[i]+1])
        neighs = list(np.sort(neighs))

        if 0 in neighs: 
            neighs.remove(0)
        if len(neighs) >1: 
            NEIGBORS.append(tuple(neighs))

    NEIGBORS = list(set(NEIGBORS))
    return NEIGBORS, label_instance_mask_orig


#def build_merged_instance_mask(single_mask_path, merge_percent=0.1, three_class = False, border_width=2): 
def build_merged_instance_mask(instance_label_mask_path, merge_percent=0.1, three_class = False, border_width=2, neighbors = []): 
    """
    neighbors: a list of tuples of neighboring cells, 
    """

    instance_mask = np.array(Image.open(instance_label_mask_path), dtype=int)
    #print(f"use merge_percent: {merge_percent}")
    if len(neighbors)==0: 
        neighbors, _ = get_list_of_neighboring_cells(instance_mask) #get_list_of_neighboring_cells(single_mask_path)

    merged_instance_mask = instance_mask.copy()
    neighbors.sort(key=len)
    neighbors = neighbors[:int(merge_percent*len(neighbors))]
    #print(len(neighbors))


    for tup in neighbors: 
        min_val = min(tup)
        for val in tup: 
            merged_instance_mask[merged_instance_mask==val] = min_val

    if three_class:     
        three_class_merged_instance_mask = build_instance_3class_mask(merged_instance_mask, border_width=border_width)
        return three_class_merged_instance_mask
    else: 
        return merged_instance_mask
    

def build_instance_3class_mask(label_instance_mask, border_width=2): 
    '''
    label_instance_mask: np.array
    '''

    #label_instance_mask = np.array(Image.open(label_instance_mask_file), dtype=int)
    three_class_instance_mask = np.zeros_like(label_instance_mask, dtype=int)

    for obj_id in np.unique(label_instance_mask): 
        if obj_id==0: 
            continue

        single_object = (label_instance_mask == obj_id).astype(int)
        if border_width<=1: 
            single_object_large = single_object
            decr_factor=1
        else: 
            incr_factor = int(border_width)//2
            decr_factor = int(border_width) - incr_factor
            single_object_large = scale_mask(single_object, increase=1, iters=incr_factor)

        single_object_small = scale_mask(single_object, increase=0, iters=decr_factor)
        object_border = single_object_large - single_object_small
        three_class_instance_mask += (2*object_border + single_object_small)

    three_class_instance_mask = three_class_instance_mask.clip(0,2)
    return three_class_instance_mask
