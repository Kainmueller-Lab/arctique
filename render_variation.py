import bpy
import sys
import os
import argparse
import glob
from tqdm import tqdm
import json

# IMPORT SOURCES
dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir)

import render as r


def parse_dataset_args():
    parser = argparse.ArgumentParser()

    # rendering parameters
    parser.add_argument("--gpu-device", type=int, default=0, help="List of GPU devices to use for rendering")
    parser.add_argument("--gpu", type=bool, default=True, help="Use GPU for rendering")
    parser.add_argument("--output-dir", type=str, default="rendered", help="Set output folder")
    parser.add_argument("--variation-dir", type=str, default="variation", help="Set output folder")                                                                                                                                                # add argument with list of all gpu devices
    parser.add_argument("--folder-parameters", type=str, default="parameters_og", help="Folder with the original rendereing parameters")
    parser.add_argument("--folder-changed-parameters", type=str, default="parameters_changed", help="Folder with the parameters to change")
    parser.add_argument("--file-changed-parameters", type=str, default="rips_0.json", help="Folder with the parameters to change")
    parser.add_argument("--render-masks", type=bool, default=False, help="Use GPU for rendering")
    
    args = parser.parse_args()
    for arg in vars(args):
        print(f"- {arg}: {getattr(args, arg)}")
    return args


def create_render_dir(args):
    '''
    Create the output directory for the rendered images
    '''
    name_variation = args.file_changed_parameters.split('.')[0]
    render_path = os.getcwd() + '/' + args.output_dir + '/' + name_variation
    dir = render_path + '/train_combined_masks/semantic'
    dir_parameters = render_path + '/parameters'
    if not os.path.exists(dir):
        os.makedirs(dir)
    if not os.path.exists(dir_parameters):
        os.makedirs(dir_parameters)
    return render_path, dir_parameters



def main():
    args = parse_dataset_args() 
    render_path, new_dir_parameters = create_render_dir(args)

    # setup parameters and load changed parameters
    dir_parameters = os.getcwd() + '/' + args.output_dir +'/' + args.folder_parameters
    path_changed_parameters = os.getcwd() + '/' + args.output_dir +'/' + args.folder_changed_parameters + '/' + args.file_changed_parameters
    parameter_files = os.listdir(dir_parameters)
    changed_parameters = json.load(open(path_changed_parameters, 'r'))

    # render individual samples
    for f_name in tqdm(parameter_files):
        
        # load original parameters
        paramters = json.load(open(dir_parameters+f'/{f_name}', 'r'))
        i = int(f_name.split('_')[1].split('.')[0])

        # change parameters with 1) changed parameters and 2) arguments and add seed
        for key, value in changed_parameters.items():
            paramters[key] = value
        for key, value in args.__dict__.items():
            paramters[key] = value
        print(f"Rendering image {i} with parameters: {paramters}")

        # save changed parameters
        with open(new_dir_parameters+f'/{f_name}', 'w') as outfile:
            json.dump(paramters, outfile)
        
        # recreate scene and render
        my_scene = r.create_scene(**paramters)
        my_manipulation = r.MainpulateScene(my_scene)
        my_manipulation.inflate_epi_cells(factor=1+paramters['epi_rescaling'])
        r.render_scene(my_scene, render_path, i, gpu=args.gpu, device=args.gpu_device, render_masks=args.render_masks)
        bpy.ops.wm.read_factory_settings(use_empty=True)

if __name__ == "__main__":
    main()