import os
import quality_control as qc


def split_data(indices, working_dir, rel_paths, dataset_name='render_v', n_test=100, check=True):
    '''
    splits the data into train and test set based on n_test images of indices
    Args:
        indices (list): list of indices to split
        working_dir (str): path to the working directory
        rel_paths (list): list of relative paths to the images
        dataset_name (str): name of the dataset
        n_test (int): number of test images
        check (bool): if True, asks for confirmation before creating the split
    '''
    indices.sort()
    test_indices = indices[:n_test]
    train_indices = indices[n_test:]

    if check:
        create_split = input(f'Create split {dataset_name} with {len(train_indices)} and {len(test_indices)} ? (y/n)')

    if create_split == 'y':
        for path in rel_paths:
            # create train and test folders
            old_path = os.path.join(working_dir, path)
            new_path = os.path.join(working_dir, dataset_name)
            new_path_train = os.path.join(new_path, 'train', path)
            new_path_test = os.path.join(new_path, 'test', path)
            if not os.path.exists(new_path_train):
                os.makedirs(new_path_train, exist_ok=True)
            if not os.path.exists(new_path_test):
                os.makedirs(new_path_test, exist_ok=True)

            # copy split into train and test paths
            files = os.listdir(old_path)
            for f in files:
                is_folder = qc.check_folder(f)
                idx = qc.get_index(f, folder=is_folder)
                if idx in test_indices:
                    source = os.path.join(old_path, f)
                    target = os.path.join(new_path_test, f)
                    os.rename(source, target)
                else:
                    source = os.path.join(old_path, f)
                    target = os.path.join(new_path_train, f)
                    os.rename(source, target)