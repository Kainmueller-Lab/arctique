import os


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
        indices (list): list of indices to be removed
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