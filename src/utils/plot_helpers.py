import numpy as np
from PIL import Image
from pathlib import Path


def reduce_single_masks(source_folder, file_names): 
    """
    Reduces the original RGB-A single masks by using only the alpha channel
 
    Args:
        source_folder (int): The folder where the masks can be found
        file_names (list): List of names of the single masks
 
    """

    for file_name in file_names: 
        mask_png = Image.open(Path(source_folder).joinpath(file_name))
        mask_np = np.array(mask_png)
        if mask_np.shape[2] != 4: 
            raise ValueError("Expected masks to have alpha channel")

        for i in range(4): 
            print(file_name, "  ", i, "  ", np.unique(mask_np[:, :, i]))
        mask_np = mask_np[:, :, -1] # keep only alpha channel
        

        mask_png_small = Image.fromarray(mask_np)
        mask_png_small.save(Path(source_folder).joinpath(file_name))
    

def build_semantic_mask(source_folder, file_names, class_identifiers): 
    """
    Combines individual cell masks to a semantic map of the full scene. 
 
    Args:
        source_folder (int): The folder where the masks can be found
        file_names (list): List of names of the single masks
        class_identifiers (list): LIst of strings identifying the different cell types as shown in file names  

 
    Returns:
        semantic_map: A png file where cells of the same type have the same pixel value
    """

    pass
