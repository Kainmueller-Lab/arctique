import os
import numpy as np
from PIL import Image
from tqdm import tqdm
import shutil
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
import json



def get_index(name, folder=True):
    '''
    gets the index of the file or folder, e.g. 'file_1.txt' 1 or '1'
    Args:
        name (str): file or folder name
        folder (bool): True if input is a folder, False otherwise
    Returns:
        int: index of the file or folder
    '''
    if folder:
        return int(name)
    else:
        return int(name.split('.')[-2].split('_')[-1])


def check_folder(file):
    '''
    checks is the file is a folder or not
    Args:
        file (str): file name
    Returns:
        bool: True if input is a folder, False otherwise
    '''
    return len(file.split('.')) == 1


def check_complete_pairs(file_path, dirs, folder_files=11):
    '''
    for all givien path checks for each sample if all files are present
    Args:
        file_path (str): path to the files
        dirs (list): list of directories
        folder_files (int): number of files in each folder
    Returns:
        list: list of indices that are complete
        list: list of indices that are incomplete
        dict: dictionary of complete paths per sample
    '''
    n_dirs = len(dirs)
    index_counts = {}
    for dir in dirs:
        files = os.listdir(file_path+dir)
        for f in files:
            is_folder = check_folder(f)
            index = get_index(f, folder=is_folder)
            if is_folder:
                n_subfiles = len(os.listdir(file_path+dir+f))
                if n_subfiles == folder_files:
                    index_counts[index] = index_counts.get(index, 0) + 1
            else:
                index_counts[index] = index_counts.get(index, 0) + 1
    complete_indices = [k for k, v in index_counts.items() if v == n_dirs]
    incomplete_indices = [k for k, v in index_counts.items() if v < n_dirs]
    return complete_indices, incomplete_indices, index_counts


def remove_files(dir, indices, check=True):
    '''
    removes files from a directory based on the indices
    Args:
        dir (str): directory path
        indices (list): list of indices to be keep
        check (bool): if True, asks for confirmation before removing files
    '''
    files = os.listdir(dir)
    folder = check_folder(files[0])
    files_to_removed = [f for f in files if get_index(f, folder=folder) not in indices]

    if check:
        remove = input(f'Remove {len(files_to_removed)} Files? (y/n) ')
    if check and remove == 'y':
        for f in files_to_removed:
            path = os.path.join(dir, f)
            if folder:
                for f in os.listdir(path):
                    os.remove(os.path.join(path, f))
                os.rmdir(path)
            else:
                os.remove(path)
    
    print(f'Removed {len(files_to_removed)} Files in {dir}')


def check_missing_renders(
        index_range,
        file_path,
        dirs=[
            'images/', 
            'masks/cytoplasm_indexing/', 
            'masks/instance_3d_indexing/', 
            'masks/instance_indexing/',
            'masks/semantic_indexing/',
            'metadata/',
            'parameters/'],
        save_path='missing_renders.txt',
        folder_files=11):
    '''
    checks for missing renders in the given range
    Args:
        index_range (list): range of indices to check, e.g. [1, 100]
        file_path (str): path to the files
        dirs (list): list of directories
        folder_files (int): number of files in each folder
    Returns:
        list: list of indices that are missing
    '''
    missing_indices = []
    complete_indices, _, _ = check_complete_pairs(file_path, dirs, folder_files)
    for i in range(index_range[0], index_range[1]+1):
        if i not in complete_indices:
            missing_indices.append(i)
    if save_path:
        with open(file_path+save_path, 'w') as f:
            for i in missing_indices:
                f.write(f'{i}\n')
    return missing_indices


class CleanArtifacts(object):
    '''
    Class to clean up artifacts in the data
    '''
    def __init__(
            self,
            file_path,
            dirs=[
                'images/', 
                'masks/cytoplasm_indexing/', 
                'masks/instance_3d_indexing/', 
                'masks/instance_indexing/',
                'masks/semantic_indexing/',
                'metadata/',
                'parameters/'],
            folder_files=11):
        '''
        Args:
            file_path (str): path to the files
            dirs (list): list of directories
            folder_files (int): number of files in each folder
        '''
        self.file_path = file_path
        self.dirs = dirs
        self.path_images = file_path+dirs[0]
        self.path_masks = file_path+dirs[4]
        self.path_masks_cyto = file_path+dirs[1]
        self.folder_files = folder_files

    def execute(self, dir_images=None, show_dir='artifacts/', check=True):
        '''
        pipeline to detect, show and remove artifacts
        Args:
            dir_images (str): path to the images
            new_dir (str): path to save the images with artifacts
            check (bool): if True, asks for confirmation before removing files
        '''
        if dir_images is None:
            dir_images = self.path_images
        remove_indices = self.detect_artifacts(dir_images)
        self.show_artifacts(remove_indices, show_dir)
        #remove_files(self.file_path+'images/', remove_indices, check=check)

    def remove(self, artifacts_dir='artifacts/', dir_images=None, check=True):
        '''
        removes the images with artifacts based on sorted out images in artifacts_dir
        Args:
            artifacts_dir (str): path to the images with artifacts
            dir_images (str): path to the images
            check (bool): if True, asks for confirmation before removing files
        '''
        if dir_images is None:
            dir_images = self.path_images
        files = os.listdir(os.path.join(self.file_path, artifacts_dir))
        img_files = os.listdir(dir_images)
        remove_indices = [get_index(f, folder=False) for f in files]
        keep_indices = [get_index(f, folder=False) for f in img_files]
        keep_indices = list(set(keep_indices) - set(remove_indices))
        remove_files(dir_images, keep_indices, check=check)
    
    #################################################
    ### pipeline functions for artifact detection ###
    #################################################
    
    def detect_artifacts(self, dir_images=None):
        '''
        detects artifacts in the images
        Args:
            dir_images (str): path to the images
        Returns:
            list: list of indices of images with artifacts
        '''
        if dir_images is None:
            dir_images = self.path_images
        
        files = os.listdir(dir_images)
        remove_indices = []
        for f in (pbar := tqdm(files)):
            path = os.path.join(dir_images, f)
            img = Image.open(path)
            img = np.array(img)
            f_mask = str(get_index(f, folder=False)) + '.tif'
            mask = Image.open(os.path.join(self.path_masks, f_mask))
            mask = np.array(mask)
            f_mask_cyto = str(get_index(f, folder=False)) + '.tif'
            mask_cyto = Image.open(os.path.join(self.path_masks_cyto, f_mask_cyto))
            mask_cyto = np.array(mask_cyto)

            if (self.detect_black_artifact(img)
                or self.detect_missing_stroma_cells(mask) 
                or self.detect_missing_crypts(img, mask) 
                or self.detect_cell_overlap(mask)
                or self.detect_missing_stroma(img, mask, mask_cyto)
                or self.detect_irregular_crypts(img, mask)):

                idx = get_index(f, folder=False)
                if idx not in remove_indices:
                    remove_indices.append(idx)
                pbar.set_description(f'Found Artifact in {f}')
        
        return remove_indices

    def show_artifacts(self, remove_indices, new_dir='artifacts/', dir_images=None):
        '''
        shows the artifacts in the images
        Args:
            remove_indices (list): list of indices of images with artifacts
            new_dir (str): path to save the images with artifacts
        '''
        if dir_images is None:
            dir_images = self.path_images
        files = os.listdir(dir_images)
        path = os.path.join(self.file_path, new_dir)
        os.makedirs(path, exist_ok=True)
        for f in (pbar := tqdm(files)):
            index = get_index(f, folder=False)
            if index in remove_indices:
                # copy the image to the new directory
                shutil.copyfile(os.path.join(dir_images, f), os.path.join(path, f))
                pbar.set_description(f'Copied {f} to {new_dir}')

    #################################
    ### detect specific artifacts ###
    #################################

    def detect_black_artifact(self, img, color_threshold=15, area_threshold=0.001):
        '''
        detects black (low color) artifacts in a given image
        Args:
            img (np.array): image to be checked, HxWx...
            color_threshold (int): threshold for the color value
            area_threshold (float): threshold for the area of the artifact
        Returns:
            bool: True if black artifact is present, False otherwise
        '''
        img = np.array(img)
        min_val = np.min(img)
        area = np.sum(img == min_val)
        size = img.shape[0] * img.shape[1]
        if min_val < color_threshold and area/size > area_threshold:
            return True
        else:
            return False

    def detect_missing_stroma_cells(self, mask, ids_stroma=[2, 3, 4, 5]):
        '''
        detects missing stroma cells in the masks
        checks if at least one stroma cell type is present
        Args:
            mask (np.array): mask to be checked, HxW
            ids_stroma (list): list of stroma cell types
        Returns:
            bool: True if stroma cells are missing, False otherwise
        '''
        mask = np.array(mask)
        classes = np.unique(mask)
        if len(set(ids_stroma).intersection(set(classes))) == 0:
            return True
        else:
            return False

    def detect_missing_crypts(self, img, mask, light_factor=0.8, ids_stroma=[2, 3, 4, 5], ids_crypts=[1]):
        '''
        detects missing crypts in the images via the light difference
        between epithelial cells and stroma cells
        Args:
            img (np.array): image to be checked, HxWx3
            mask (np.array): mask to be checked, HxW
            light_factor (float): factor for the light difference
            ids_stroma (list): list of stroma cell types
            ids_crypts (list): list of crypt cell types
        Returns:
            bool: True if crypts are missing, False otherwise
        '''
        mask = np.array(mask)
        classes = np.unique(mask)

        if len(set(ids_crypts).intersection(set(classes))) != 0:
            if len(set(ids_stroma).intersection(set(classes))) != 0:
                mask_stroma = np.zeros_like(mask)
                for id in ids_stroma:
                    mask_stroma += mask == id
                mask_crypts = np.zeros_like(mask)
                for id in ids_crypts:
                    mask_crypts += mask == id
                mean_stroma = np.mean(img[mask_stroma==1])
                mean_crypts = np.mean(img[mask_crypts==1])
                if mean_stroma/mean_crypts < light_factor:
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

    def detect_cell_overlap(
            self, mask, n_touching=2, ids_stroma=[2, 3, 4, 5], ids_crypts=[1]):
        '''
        detects if touching epithelial and stroma cells are present
        Args:
            mask (np.array): mask to be checked, HxW
            n_touching (int): number of touching pixels
            ids_stroma (list): list of stroma cell types
            ids_crypts (list): list of crypt cell types
        Returns:
            bool: True if touching cells are present, False otherwise
        '''
        mask = np.array(mask)
        mask_stroma = np.zeros_like(mask)
        for id in ids_stroma:
            mask_stroma += mask == id
        mask_crypts = np.zeros_like(mask)
        for id in ids_crypts:
            mask_crypts += mask == id
        mask_stroma = 1-ndi.binary_erosion(1-mask_stroma, border_value=1)
        overlap = np.sum(mask_stroma * mask_crypts)
        if overlap > n_touching:
            return True
        else:
            return False  

    def detect_missing_stroma(self, img, mask, mask_cyto, ids_stroma=[2, 3, 4, 5]):
        '''
        by measuring the maximal value around the edges of the stroma cells
        checks if there is stroma around the edges
        Args:
            img (np.array): image to be checked, HxWx3
            mask (np.array): mask to be checked, HxW
            mask_cyto (np.array): cytoplasm mask to be checked, HxW
            ids_stroma (list): list of stroma cell types
        Returns:
            bool: True if there is stroma around the edges, False otherwise
        '''
        mask = np.copy(np.array(mask))
        mask[mask_cyto != 0] = mask_cyto[mask_cyto != 0]
        mask_stroma = np.zeros_like(mask)
        for id in ids_stroma:
            mask_stroma += mask == id

        # check if there is stroma around the edges
        mask_stroma = ((1-ndi.binary_erosion(1-mask_stroma, border_value=1, iterations=2)) - mask_stroma)==1
        max_vals = np.max(img, axis=(0,1))
        cut_img = img[mask_stroma]
        max_value = np.zeros(cut_img.shape[0])
        for i in range(3):
            max_value = max_value + (cut_img[:, i] > (max_vals[i]-15))
        max_value = max_value == 3
        fraction = np.mean(max_value)
        if fraction > 0.9:
            return True
        else:
            return False
        
    def detect_irregular_crypts(self, img, mask, ids_crypts=[1]):
        '''
        measures minimal value around epithelial cells compared to median value
        Args:
            img (np.array): image to be checked, HxWx3
            mask (np.array): mask to be checked, HxW
            ids_crypts (list): list of crypt cell types
        Returns:
            bool: True if there are high irregularities, False otherwise
        '''
        mask = np.array(mask)
        mask_crypts = np.zeros_like(mask)
        classes = np.unique(mask)
        if len(set(ids_crypts).intersection(set(classes))) == 0:
            return False
        for id in ids_crypts:
            mask_crypts += mask == id
        mask_crypts = ((1-ndi.binary_erosion(1-mask_crypts, border_value=1, iterations=1)) - mask_crypts)==1
        cut_img = img[mask_crypts]
        std = np.mean(np.min(cut_img, axis=0)[:3]/np.median(cut_img, axis=0)[:3])
        if std < 0.1:
            return True
        else:
            return False


def correct_semantic_assignment(
        file_path,
        dirs_semantic='masks/semantic_indexing/',
        dir_cyto='masks/cytoplasm_indexing/',
        dir_metadata='metadata/',
        assignment={'EPI': 1, 'PLA': 2, 'LYM': 3, 'EOS': 4, 'FIB': 5}):
    '''
    corrects the indexing of the semantic masks based on the metadata
    and assigns consistent values to the classes as given in the assignment
    Args:
        file_path (str): path to the files
        dirs_semantic (str): path to the semantic masks
        dir_cyto (str): path to the cytoplasm masks
        dir_metadata (str): path to the metadata
        assignment (dict): dictionary to assign the classes to the new values
    '''
    path_semantic = file_path+dirs_semantic
    path_cyto = file_path+dir_cyto
    path_metadata = file_path+dir_metadata
    suffix='_corrected'
    new_path_semantic = path_semantic[:-1]+suffix
    new_path_cyto = path_cyto[:-1]+suffix
    os.makedirs(new_path_semantic, exist_ok=True)
    os.makedirs(new_path_cyto, exist_ok=True)
    files = os.listdir(path_semantic)

    for f in tqdm(files):
        idx = get_index(f, folder=False)

        # construct look-up table how indices are assigned
        false_assignment = {}
        with open(os.path.join(path_metadata, f'metadata_{idx}.json'), 'r') as file:
            metadata = json.load(file)
            for m in metadata:
                id_type_idx = int(metadata[m]['ID_Type'])
                id_name = metadata[m]['Type']
                if id_name not in false_assignment.values():
                    false_assignment[id_type_idx] = id_name
            false_assignment = dict(sorted(false_assignment.items()))
            mapping = {}
            for key, val in false_assignment.items():
                if val in assignment.keys():
                    mapping[key] = assignment[val]
        
        # load the masks and map the values
        mask = Image.open(os.path.join(path_semantic, f))
        mask = np.array(mask)
        mask_corrected = np.zeros_like(mask)
        for c in np.unique(mask):
            if c in mapping.keys():
                mask_corrected[mask == c] = mapping[c]

        # same for cytoplasm masks
        mask_cyto = Image.open(os.path.join(path_cyto, f))
        mask_cyto = np.array(mask_cyto)
        mask_cyto_corrected = np.zeros_like(mask_cyto)
        for c in np.unique(mask_cyto):
            if c in mapping.keys():
                mask_cyto_corrected[mask_cyto == c] = mapping[c]
        
        # save the corrected masks
        mask_corrected = Image.fromarray(mask_corrected)
        mask_corrected.save(os.path.join(new_path_semantic, f))
        mask_cyto_corrected = Image.fromarray(mask_cyto_corrected)
        mask_cyto_corrected.save(os.path.join(new_path_cyto, f))
        