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

import src.arrangement.arrangement as arr 
import src.objects.cells as cells
import src.objects.tissue as tissue
import src.shading.materials as materials
import src.scene as scene
import src.utils as utils
import src.utils.geometry as geom
import src.utils.surface_filling as sf
import src.objects.tissue_architecture as arch
import src.utils.helper_methods as hm
import src.utils.geometry as geo

# this next part forces a reload in case you edit the source after you first start the blender session
#import imp
import importlib as imp # imp module is deprecated since python 3.12
imp.reload(arr)
imp.reload(cells)
imp.reload(tissue)
imp.reload(materials)
imp.reload(scene)
imp.reload(utils)
imp.reload(geom)
imp.reload(sf)
imp.reload(arch)
imp.reload(hm)
imp.reload(geo)



def parse_dataset_args():
    parser = argparse.ArgumentParser()
    
    # RENDERING PARAMETERS                                                                                                                                                        # add argument with list of all gpu devices
    parser.add_argument("--gpu-device", type=int, default=0, help="List of GPU devices to use for rendering")
    parser.add_argument("--gpu", type=bool, default=True, help="Use GPU for rendering")
    parser.add_argument("--output-dir", type=str, default="rendered", help="Set output folder")
    parser.add_argument("--start-idx", type=int, default=400, help="Dataset size")
    parser.add_argument("--n-samples", type=int, default=200, help="Dataset size")

    # DATASET PARAMETERS
    # tissue
    parser.add_argument("--tissue-thickness", type=float, default=0.05, help="Tissue thickness")
    parser.add_argument("--tissue-size", type=float, default=1.28, help="Tissue size")
    parser.add_argument("--tissue-location", type=tuple, default=(0, 0, 0.5), help="Tissue location")
    parser.add_argument("--tissue-padding", type=float, default=0.2, help="Tissue padding")
    parser.add_argument("--tissue-rips", type=float, default=-0.5, help="Degree of rip like structures in tissue")
    parser.add_argument("--tissue-rips-std", type=float, default=0.2, help="Degree of rip like structures in tissue")
    parser.add_argument("--stroma-intensity", type=float, default=1, help="Degree of rip like structures in tissue")
    parser.add_argument("--noise-seed-shift", type=float, default=0, help="Degree of rip like structures in tissue")
    
    # nuclei
    parser.add_argument("--epi-number", type=int, default= 150, help="number of surface cells") # 150
    parser.add_argument("--filler-scale", type=float, default= 0.8, help="Scale of the size of smaller filler nuclei w.r.t to the original nuclei size")
    parser.add_argument("--stroma-density", type=int, default= 0.7, help="density in stroma") # 1200
    parser.add_argument("--ratios", type=list, default=[0, 0.2, 0.4, 0.2, 0.2], help="ratios of different cell types")
    parser.add_argument("--surf_scale", type=tuple, default=(0.8, 0.5, 1), help="Surface scale")
    parser.add_argument("--delete-fraction", type=list, default=[0, 0, 0, 0, 0], help="ratios of different cell types")
    parser.add_argument("--nuclei-intensity", type=float, default=1, help="overall intensity of nuclei")
    parser.add_argument("--mix-factor", type=float, default=0, help="overall intensity of nuclei")

    # TODO add manipulation paramter of same scene

    #other default value for --output_dir: "/Volumes/ag_kainmueller/vguarin/synthetic_HE" via internal VPN
    
    args = parser.parse_args()
    
    for arg in vars(args):
        print(f"- {arg}: {getattr(args, arg)}")
    return args


def create_scene(
        tissue_thickness = 0.05, tissue_size = 1.28, tissue_location = (0, 0, 0.5),
        tissue_rips = -0.5, tissue_rips_std = 0.1,
        tissue_padding = 0.5, epi_count = 80, stroma_density = 0.5, 
        ratios = [0, 0.3, 0.4, 0.2, 0.1],
        seed=0, **kwargs):
    '''
    creates a tissue crop with cells and nuclei
    Args:
        tissue_thickness: float, thickness of the tissue [unit: 10 micrometers]
        tissue_size: float, size of the tissue [unit: 10 micrometers]
        tissue_location: tuple, location of the tissue [unit: 10 micrometers]
        tissue_padding: float, padding of the tissue [unit: 10 micrometers]
        surf_number: int, number of surface cells
        filler_scale: float, scale of size of smaller filler nuclei w.r.t to original nuclei size
        number: int, number of volume cells
        ratios: list, ratios of different cell types
        vol_scale: tuple, volume scale
        surf_scale: tuple, surface scale
    Returns:
        my_scene: BioMedicalScene object
    '''
    scene.BioMedicalScene.clear()
        
    # 1) initialize microscope objects and add to scene
    params_cell_shading = {
        'PLA': {
            'Nucleus': {'name': 'Nucleus_PLA', 'color': (0.315, 0.003, 0.631, 1), 'staining_intensity': 170},
            'Cytoplasm': {'name': 'Cytoplasm_PLA', 'color': (0.456, 0.011, 0.356, 1), 'staining_intensity': 140}},
        'LYM': {
            'Nucleus': {'name': 'Nucleus_LYM', 'color': (0.315, 0.003, 0.631, 1), 'staining_intensity': 400}},
        'EOS': {
            'Nucleus': {'name': 'Nucleus_EOS', 'color': (0.315, 0.003, 0.631, 1), 'staining_intensity': 170},
            'Cytoplasm': {'name': 'Cytoplasm_EOS', 'color': (0.605, 0.017, 0.043, 1), 'staining_intensity': 140}},
        'FIB': {
            'Nucleus': {'name': 'Nucleus_FIB', 'color': (0.315, 0.003, 0.631, 1), 'staining_intensity': 200}},
        'EPI': {
            'Nucleus': {'name': 'Nucleus_EPI', 'color': (0.315, 0.003, 0.631, 1), 'staining_intensity': 100}}}
    my_materials = materials.Material(
        seed=seed, cell_type_params=params_cell_shading, tissue_rips=tissue_rips, tissue_rips_std=tissue_rips_std)
    my_tissue = tissue.Tissue(
        my_materials.muscosa, thickness=tissue_thickness,
        size=tissue_size, location=tissue_location)
    my_light_source = scene.LightSource(material=my_materials.light_source)
    my_camera = scene.Camera()
    my_scene = scene.BioMedicalScene(my_light_source, my_camera)
    my_scene.add_cell_params(params_cell_shading)
    my_scene.add_tissue(tissue=my_tissue.tissue)

    # 2) create macrostructures in tissue block, rotate and scale them and cut them
    tissue_arch = arch.TissueArch(seed=seed)
    tissue_arch.random_crop(my_tissue.tissue)
    macro_structure = tissue_arch.get_architecture()
    crypt, crypt_vol_1, crypt_vol_2, mucosa = macro_structure
    my_scene.bound_architecture(
        volumes=[crypt_vol_1, crypt_vol_2, mucosa], surfaces=[crypt],
        padding=tissue_padding)
    vol_goblet = hm.copy_object(crypt_vol_2, 'vol_goblet')
    extended_stroma = hm.copy_object(mucosa, 'extended_stroma')
    hm.add_boolean_modifier(extended_stroma, crypt_vol_1, operation='UNION', name='add epi to stroma', apply=True)
    hm.add_boolean_modifier(vol_goblet, extended_stroma, name='Remove inner volume', apply=True)

    # 3) populate scene with nuclei/cells
    # add epi volume filling
    crypt_goblet = arr.VoronoiFill(vol_goblet, extended_stroma, cells.CellType.GOB)
    crypt_fill = arr.VoronoiFill(crypt_vol_1, mucosa, cells.CellType.EPI)
    my_scene.add_arrangement(crypt_fill) # NOTE: 200 nuclei take about 40 s
    my_scene.add_arrangement(crypt_goblet)

    # Add volume filling
    # add tissue padding befor filling
    MIX_TYPES = [
        cells.CellType.MIX,
        cells.CellType.PLA, 
        cells.CellType.LYM, 
        cells.CellType.EOS, 
        cells.CellType.FIB]
    stroma_fill = hm.copy_object(mucosa, name='stroma_fill')
    hm.convert2mesh(stroma_fill)
    volume_fill = arr.VolumeFill(
        mucosa, stroma_density, MIX_TYPES, ratios, strict_boundary=True, seed=seed)
    my_scene.add_arrangement(volume_fill, bounding_mesh=stroma_fill) # NOTE: 240 nuclei take about 20 s
    #my_scene.cut_cytoplasm_nuclei()

    # 4) cut objects and add staining
    my_scene.add_cell_params(params_cell_shading)
    my_scene.delete_cells()
    my_scene.remove_goblet_volume(crypt_vol_2)
    my_scene.remove_cells_volume(mucosa)
    my_scene.cut_cells()
    my_scene.cut_tissue()
    my_scene.add_tissue_staining(materials=[my_materials.muscosa, my_materials.crypt_staining])
    my_scene.add_nuclei_mask(material=my_materials.nuclei_mask)
    my_scene.add_staining_to_cell(materials=my_materials.cell_staining)

    # 5) hide non cell objects
    goblet_cells = []
    for cell in my_scene.cell_objects:
        cell_type = cell.name.split('_')[-2]
        if cell_type == 'GOB':
            goblet_cells.append(cell)
    for obj in [crypt, crypt_vol_1, stroma_fill, vol_goblet, extended_stroma]+goblet_cells:
        obj.hide_viewport = True
        obj.hide_render = True

    return my_scene


class mainpulate_scene(object):
    def __init__(self, my_scene):
        self.my_scene = my_scene
    
    def change_tissue_thickness(self, thickness):
        pass

    def delete_objects(self, objects):
        pass

    def change_staining(self, cell_object, staining_intensity, staining_color):
        pass

    def change_tissue_staining(self, materials):
        pass

    def inflate_cells(self):
        pass
    # TODO 
    # fix indexing problem



def recreate_scene(**kwargs):
    '''
    recreates a scene from parameters
    Args:
        kwargs: dict, parameters of the scene
    Returns:
        my_scene: BioMedicalScene object
    '''
    my_scene = create_scene(**kwargs)
    return my_scene


def render_scene(my_scene, render_path, sample_name, gpu=True, device=0, output_shape=(512, 512), max_samples=256):
    '''
    renders a scene
    Args:
        my_scene: BioMedicalScene object
        render_path: str, path to save renders
        sample_name: int, sample name
        gpu: bool, use gpu for rendering
        devices: list, list of gpu devices to use for rendering
        output_shape: tuple, dimensions of output
        max_samples: int, number of samples for rendering
    '''
    
    # set render engine
    bpy.context.scene.render.engine = 'CYCLES'
    if gpu:
        bpy.context.scene.cycles.device = 'GPU'
    bpy.context.preferences.addons['cycles'].preferences.compute_device_type = "CUDA"
    bpy.context.preferences.addons["cycles"].preferences.get_devices()
    print(bpy.context.preferences.addons["cycles"].preferences.compute_device_type)
    for j, d in enumerate(bpy.context.preferences.addons["cycles"].preferences.devices):
        if j == device:
            d["use"] = True
        else:
            d["use"] = False
        print(d["name"], d["use"])
    
    my_scene.sample_name = sample_name
    
    my_scene.render(filepath = render_path,  # where to save renders
                    scene = True, # if true scene is rendered
                    single_masks = True, # if true single cell masks are rendered
                    semantic_mask = True, # if true semantic mask is generated
                    instance_mask = True, # if true instance mask is generated
                    depth_mask = False, # if true depth mask is generated
                    obj3d = False, # if true scene is saved as 3d object
                    output_shape = output_shape, # dimensions of output
                    max_samples = max_samples) # number of samples for rendering. Fewer samples will render more quickly. Default is 1024


def main():
    args = parse_dataset_args() 
    
    # create output directory
    if args.output_dir == 'rendered':
        render_path = os.getcwd() + '/rendered'
    else:
        render_path = args.output_dir
    print(render_path)
    dir = render_path + '/train_combined_masks/semantic'
    dir_parameters = render_path + '/parameters'
    if not os.path.exists(dir):
        os.makedirs(dir)
    if not os.path.exists(dir_parameters):
        os.makedirs(dir_parameters)
    

    # render individual samples
    for i in tqdm(range(args.start_idx, args.start_idx + args.n_samples)):
        # paramters = {
        #     'tissue_thickness': args.tissue_thickness, 'tissue_size': args.tissue_size,
        #     'tissue_location': args.tissue_location, 'tissue_padding': args.tissue_padding,
        #     'tissue_rips': args.tissue_rips,
        #     'epi_count': args.epi_number, 'stroma_density': args.stroma_density, 'ratios': args.ratios,
        #     'seed': i}
        paramters = {'seed': i}
        for key, value in args.__dict__.items():
            if key not in paramters.keys():
                paramters[key] = value
        with open(dir_parameters+f'/parameters_{i}.json', 'w') as outfile:
            json.dump(paramters, outfile)
        my_scene = create_scene(**paramters)
        render_scene(my_scene, render_path, i+1, gpu=args.gpu, device=args.gpu_device)
        bpy.ops.wm.read_factory_settings(use_empty=True)

if __name__ == "__main__":
    main()