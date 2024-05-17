import numpy as np
from PIL import Image
from pathlib import Path

import os

def make_color_palette(n_colors, bg_col = (0,0,0), seed=42): 
    np.random.seed(seed)
    colors = [bg_col]
    for i in range(n_colors): 
        not_unique = True
        while not_unique: # generate random RGB color and ensure no repetition
            col = tuple(np.random.randint(0, 255, size=3)) 
            if col not in colors: 
                colors.append(col)
                not_unique =False
    # palette = []
    # for color in colors:
    #     palette.extend(color)

    return colors

def put_palette(image_arr, palette, ids=None):
    '''
    image_arr: np.array, image array, shape (H, W) with pixel values in [0, 255]
    palette: list, list of RGB values
    '''
    if ids is None:
        ids = np.unique(image_arr)
    new_image_arr = np.zeros((image_arr.shape[0], image_arr.shape[1], 3), dtype=np.uint8)
    for i, idx in enumerate(ids):
            new_image_arr[image_arr == idx] = palette[i]
    return new_image_arr

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
        mask_png_small.close()

    

def build_semantic_mask(source_folder, cell_info_dicts, file_name:int=0, palette = None): 
    """
    Combines individual cell masks to a semantic map of the full scene. 
 
    Args:
        source_folder (int): The folder where the masks can be found
        cell_info_tuples (list): List of tuples each of which contains: cell_id, cell_type, cell_mask_filename 

 
    Returns:
        semantic_map: A png file where cells of the same type have the same pixel value
    """

    unique_cell_types = set([cit["Type"] for cit in cell_info_dicts]) # identify unique cell types 
    cell_type_dict = {uct : (i+1) for i, uct in enumerate(unique_cell_types)} # assign unique id to each cell type
        
    if not os.path.exists(source_folder + f'/train/{file_name}/masks/'):
        raise TypeError('The creation of masks might have been unsuccessful')
    
    for cell_counter, cell_info_tuple in enumerate(cell_info_dicts): 
        cell_mask_file = cell_info_tuple["Filename"]
        cell_type = cell_info_tuple["Type"]
        mask_png = Image.open(cell_mask_file)
        mask_np = np.array(mask_png).astype(np.float64)

        if cell_counter == 0:
            semantic_mask = np.zeros_like(mask_np) # semantic_mask has same shape as the masks
    
        # assign pixel value based on cell class
        semantic_mask += mask_np*(cell_type_dict[cell_type]/255.) # in mask the object has pixel value 255 (backgrund is 0)
        
        # generate unique colors for each class (background is black by default)
        if palette is None: 
            palette = make_color_palette(len(cell_type_dict.keys()))
        # generate Image object from array, needs to be converted to uint8 to avoid aliasing
        colored_instance_mask = Image.fromarray(semantic_mask.astype(np.uint8))
        # assign color palette to image
        colored_instance_mask.putpalette(palette)
        # save to png
    
    new_source_folder = source_folder + '/train_combined_masks/semantic/'
    if not os.path.exists(new_source_folder):
        os.makedirs(new_source_folder)
    np.save(str(Path(new_source_folder).joinpath(f"{file_name}.npy")), semantic_mask)
    colored_instance_mask.save(str(Path(new_source_folder).joinpath(f"{file_name}.png")))
    # colored_instance_mask.close()


# this works for 2d but must be adapted for 3d
# def build_instance_mask(source_folder, file_names): 
#     """
#     Combines individual cell masks to a semantic map of the full scene. 
 
#     Args:
#         source_folder (int): The folder where the masks can be found
#         file_names (list): List of names of the single masks
#         class_identifiers (list): LIst of strings identifying the different cell types as shown in file names  

 
#     Returns:
#         instance_mask: A png file where each cell object has a different pixel value
#     """

#     for file_idx, file_name in enumerate(file_names):   
#         mask_png = Image.open(file_name) # open mask png
#         mask_np = np.array(mask_png).astype(np.float64) # convert to numpy


#         if file_idx == 0:
#             instance_mask = np.zeros_like(mask_np) # instance_mask has same shape as the masks
    
#         instance_mask += mask_np*((file_idx+1)/255.) # in mask the object has pixel value 255 (backgrund is 0)

#     np.save(str(Path(source_folder).joinpath("instance_mask.npy")), instance_mask)

#     # generate unique colors for each instance (background is black by default)
#     palette = make_color_palette(len(np.unique(instance_mask)))
#     # generate Image object from array, needs to be converted to uint8 to avoid aliasing
#     colored_instance_mask = Image.fromarray(instance_mask.astype(np.uint8))
#     # assign color palette to image
#     colored_instance_mask.putpalette(palette)
#     # save to png
#     colored_instance_mask.save(str(Path(source_folder).joinpath("instance_mask.png")))
        

def build_instance_mask(source_folder, cell_info_dicts, file_name:int=0, palette=None): 
   # source_folder, cell_info_tuples, file_name="semantic_mask", palette = None): 
    """
    Combines individual cell masks to an instance map of the full scene. 
 
    Args:
        source_folder (int): The folder where the masks can be found
        file_names (list): List of names of the single masks
        class_identifiers (list): LIst of strings identifying the different cell types as shown in file names  

 
    Returns:
        instance_mask: A png file where each cell object has a different pixel value
    """
    cell_ID_list = [c["ID"] for c in cell_info_dicts] # make list of all cell ids 
    cell_ID_dict = {cid : (i+1) for i, cid in enumerate(cell_ID_list)} # assign unique integer to each cell id

    if not os.path.exists(source_folder + f'/train/{file_name}/masks/'):
        raise TypeError('The creation of masks might have been unsuccessful')
    
    #print(cell_info_dicts)
    #for file_idx, file_name in enumerate(file_names):   ´
    for cell_counter, cell_info_tuple in enumerate(cell_info_dicts): 
        cell_id = cell_info_tuple["ID"]
        cell_mask_file = cell_info_tuple["Filename"]
        mask_png = Image.open(cell_mask_file) # open mask png
        mask_np = np.array(mask_png).astype(np.float64) # convert to numpy

        if cell_counter==0: 
            instance_mask = np.zeros_like(mask_np) # instance_mask has same shape as the masks
    
        instance_mask += mask_np*((cell_ID_dict[cell_id])/255.) # in mask the object has pixel value 255 (backgrund is 0)

    new_source_folder = source_folder + "/train_combined_masks/instance/"
    if not os.path.exists(new_source_folder):
        os.makedirs(new_source_folder)
    np.save(str(Path(new_source_folder).joinpath(f"{file_name}.npy")), instance_mask)

    # generate unique colors for each instance (background is black by default)
    # palette = make_color_palette(len(np.unique(instance_mask)))
    # generate Image object from array, needs to be converted to uint8 to avoid aliasing
    colored_instance_mask = Image.fromarray(instance_mask.astype(np.uint8))
    # assign color palette to image 
    colored_instance_mask.putpalette(palette)
    # save to png
    colored_instance_mask.save(str(Path(new_source_folder).joinpath(f"{file_name}.png")))
    # colored_instance_mask.close()


def remove_single_masks(file_names): 
    """
    Removes individual cell masks.
    """

    for file_name in file_names: 
        os.remove(file_name)


def build_gif(png_file_names, gif_file_name): 
    """
    Combines multiple PNG images to a GIF
    """

    frames = []
    for png_file in png_file_names: 
        frames.append(Image.open(png_file))

    frame_one = frames[0]
    frame_one.save(gif_file_name,
                   format="GIF", 
                   append_images=frames[1:], #list of images to append as additional frames.
                   save_all=True, 
                   duration=500, #display duration of each frame, in milliseconds
                   loop=0)#Number of times to repeat the animation'