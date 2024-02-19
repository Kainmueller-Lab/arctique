import numpy as np
from PIL import Image
from pathlib import Path

import os


def make_color_palette(n_colors, bg_col = (0,0,0)): 
    colors = [bg_col]
    for i in range(n_colors): 
        not_unique = True
        while not_unique: # generate random RGB color and ensure no repetition
            col = tuple(np.random.randint(0, 255, size=3)) 
            if col not in colors: 
                colors.append(col)
                not_unique =False
    palette = []
    for color in colors:
        palette.extend(color)

    return palette

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
        mask_png = Image.open(file_name)
        mask_np = np.array(mask_png)
        if mask_np.shape[2] != 4: 
            raise ValueError("Expected masks to have alpha channel")

        mask_np = mask_np[:, :, -1] # keep only alpha channel
        
        
        mask_png_small = Image.fromarray(mask_np)
        mask_png_small.save(file_name)
    

def build_semantic_mask(source_folder, cell_info_tuples): 
    """
    Combines individual cell masks to a semantic map of the full scene. 
 
    Args:
        source_folder (int): The folder where the masks can be found
        cell_info_tuples (list): List of tuples each of which contains: cell_id, cell_type, cell_mask_filename 

 
    Returns:
        semantic_map: A png file where cells of the same type have the same pixel value
    """

    unique_cell_types = set([cit[1] for cit in cell_info_tuples]) # identify unique cell types 
    cell_type_dict = {uct : (i+1) for i, uct in enumerate(unique_cell_types)} # assign unique id to each cell type


    for cell_counter, cell_info_tuple in enumerate(cell_info_tuples): 
        cell_mask_file = cell_info_tuple[2]
        cell_type = cell_info_tuple[1]
        mask_png = Image.open(cell_mask_file)
        mask_np = np.array(mask_png).astype(np.float64)

        if cell_counter == 0:
            semantic_mask = np.zeros_like(mask_np) # semantic_mask has same shape as the masks
    
        # assign pixel value based on cell class
        semantic_mask += mask_np*(cell_type_dict[cell_type]/255.) # in mask the object has pixel value 255 (backgrund is 0)

        np.save(str(Path(source_folder).joinpath("semantic_mask.npy")), semantic_mask)

        # generate unique colors for each class (background is black by default)
        palette = make_color_palette(len(cell_type_dict.keys()))
        # generate Image object from array, needs to be converted to uint8 to avoid aliasing
        colored_instance_mask = Image.fromarray(semantic_mask.astype(np.uint8))
        # assign color palette to image
        colored_instance_mask.putpalette(palette)
        # save to png
        colored_instance_mask.save(str(Path(source_folder).joinpath("semantic_mask.png")))



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
        mask_png = Image.open(file_name) # open mask png
        mask_np = np.array(mask_png).astype(np.float64) # convert to numpy


        if file_idx == 0:
            instance_mask = np.zeros_like(mask_np) # instance_mask has same shape as the masks
    
        instance_mask += mask_np*((file_idx+1)/255.) # in mask the object has pixel value 255 (backgrund is 0)

    np.save(str(Path(source_folder).joinpath("instance_mask.npy")), instance_mask)

    # generate unique colors for each instance (background is black by default)
    palette = make_color_palette(len(np.unique(instance_mask)))
    # generate Image object from array, needs to be converted to uint8 to avoid aliasing
    colored_instance_mask = Image.fromarray(instance_mask.astype(np.uint8))
    # assign color palette to image
    colored_instance_mask.putpalette(palette)
    # save to png
    colored_instance_mask.save(str(Path(source_folder).joinpath("instance_mask.png")))




def remove_single_masks(file_names): 
    '''
    Removes individual cell masks.
    '''

    for file_name in file_names: 
        os.remove(file_name)

