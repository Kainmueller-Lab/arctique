import numpy as np
from PIL import Image
from pathlib import Path


def reduce_single_masks(source_folder, file_names): 
    """
    Reduces the original RGB-A single masks by using only the alpha channel. 

    This is a workaround because when masks are rendered the black background sometimes has pixel values >0
    and the cell object sometimes has pixel values <255. Therefore we set the background to transparent 
    and thereby guarantee that the alpha channel will only have two pixel values: 0 for background and 255 for the cell object. 
 
    Args:
        source_folder (int): The folder where the masks can be found
        file_names (list): List of names of the single masks
 
    """

    for file_name in file_names: 
        mask_png = Image.open(Path(source_folder).joinpath(file_name))
        mask_np = np.array(mask_png)
        if mask_np.shape[2] != 4: 
            raise ValueError("Expected masks to have alpha channel")

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



def build_instance_mask(source_folder, file_names): 
    """
    Combines individual cell masks to a semantic map of the full scene. 
 
    Args:
        source_folder (int): The folder where the masks can be found
        file_names (list): List of names of the single masks
        class_identifiers (list): LIst of strings identifying the different cell types as shown in file names  

 
    Returns:
        instance_mask: A png file where each cell object has a different pixel value
    """

    for file_idx, file_name in enumerate(file_names):   

        mask_png = Image.open(Path(source_folder).joinpath(file_name))
        mask_np = np.array(mask_png)

        if file_idx == 0: 
            # instance_mask is the same shape as the masks
            instance_mask = np.zeros_like(mask_np) 
    
        instance_mask += mask_np*(file_idx/255.)

        np.save(instance_mask, file=str(Path(source_folder).joinpath("instance_mask.npy")))
        
        instance_mask_png = Image.fromarray(instance_mask)
        instance_mask.save(Path(source_folder).joinpath("instance_mask.png"))
