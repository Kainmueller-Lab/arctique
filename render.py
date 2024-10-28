import bpy
import sys
import os
import argparse
import glob
from tqdm import tqdm
import json
import time
import numpy as np

# IMPORT SOURCES
dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir)

import src.arrangement.arrangement as arr 
import src.arrangement.surface_filling as sf
import src.objects.cells as cells
import src.objects.tissue as tissue
import src.shading.materials as materials
import src.scene as scene
import src.utils as utils
import src.utils.geometry as geom
import src.objects.tissue_architecture as arch
import src.utils.helper_methods as hm
import src.utils.geometry as geo
import src.objects.cells as cells
import src.arrangement.deformation as defo

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
    parser.add_argument("--start-idx", type=int, default=0, help="Dataset size")
    parser.add_argument("--index-list", type=str, default='', help="ony considered if len>1, List of indices to render")
    parser.add_argument("--n-samples", type=int, default=100, help="Dataset size")
    parser.add_argument("--base-16bit", type=int, default=55, help="Base for 16 bit images (55 or 257)")

    # DATASET PARAMETERS
    # tissue parameters (when adaptive they orientate at a default tissue thickness of 0.05 and a default tissue size of 1.28)
    parser.add_argument("--tissue-thickness", type=float, default=0.05, help="Tissue thickness")
    parser.add_argument("--tissue-thickness_lb", type=float, default=0.025, help="Tissue thickness")
    parser.add_argument("--tissue-size", type=float, default=1.28, help="Tissue size in 100 microns")  # 1.28
    parser.add_argument("--scale_scene", type=tuple, default=1.15, help="Scale scene (to adapt to different cell scalings)")
    parser.add_argument("--tissue-color", type=tuple, default=(0.64, 0.347, 0.642, 1), help="Tissue location")
    parser.add_argument("--nucleus-color", type=tuple, default=(0.315, 0.003, 0.48, 1), help="Tissue location")
    parser.add_argument("--color-variation", type=tuple, default=(0.0, 0.0, 0.0), help="std of HSV color variation")
    parser.add_argument("--red-base", type=tuple, default=(0.605, 0.017, 0.216, 1), help="Tissue location")
    parser.add_argument("--red-shift", type=tuple, default=(1, 1), help="(min, max) Tissue location")
    parser.add_argument("--tissue-location", type=tuple, default=(0, 0, 0.5), help="Tissue location")
    parser.add_argument("--tissue-padding", type=float, default=0.15, help="Tissue padding")  # 0.2
    parser.add_argument("--tissue-rips", type=float, default=0.5, help="Degree of rip like structures in tissue")
    parser.add_argument("--tissue-rips-std", type=float, default=0.25, help="Degree of rip like structures in tissue")
    parser.add_argument("--tissue-rips-curl", type=tuple, default=(0, 1), help="(min ,max) Degree of rip curl")
    parser.add_argument("--stroma-intensity", type=float, default=0.7, help="Degree of rip like structures in tissue")
    parser.add_argument("--noise-seed-shift", type=float, default=0, help="Degree of rip like structures in tissue")
    parser.add_argument("--light-source-brightness", type=float, default=32, help="Degree of rip like structures in tissue")
    parser.add_argument("--adaptiv-brightness", type=bool, default=True, help="Use GPU for rendering")
    parser.add_argument("--focal-offset", type=float, default=0, help="Degree of rip like structures in tissue")
    parser.add_argument("--over-staining", type=tuple, default=(0.2, 1), help="Degree of overstaining")
    parser.add_argument("--goblet-intensity", type=tuple, default=(0, 1), help="Degree goblet cell staining")
    parser.add_argument("--nuclei-deformation", type=float, default=1, help="Degree of nuclei deformation")

    # nuclei
    parser.add_argument("--epi-number", type=int, default=300, help="number of surface cells") # 150
    parser.add_argument("--filler-scale", type=float, default=0.8, help="Scale of the size of smaller filler nuclei w.r.t to the original nuclei size")
    parser.add_argument("--stroma-density", type=int, default= 1, help="density in stroma") # 0.5, 1200
    parser.add_argument("--ratios", type=list, default=[0, 0.2, 0.5, 0.2, 0.1], help="ratios of different cell types (MIX, PLA, LYM, EOS, FIB); LYM should be at least 0.8 for best results")
    parser.add_argument("--surf_scale", type=tuple, default=(0.8, 0.5, 1), help="Surface scale")
    parser.add_argument("--delete-fraction", type=list, default=[0, 0, 0, 0, 0], help="ratios of different cell types")
    parser.add_argument("--nuclei-intensity", type=float, default=0.7, help="overall intensity of nuclei") # TODO
    parser.add_argument("--mix-factor", type=float, default=0, help="overall intensity of nuclei") # TODO
    parser.add_argument("--epi-rescaling", type=int, default=0, help="overall intensity of nuclei") # TODO
    parser.add_argument("--mix-cyto", type=float, default=0, help="overall intensity of nuclei") # TODO
    parser.add_argument("--red-points-strength", type=float, default=0, help="Degree of rip like structures in tissue")
    
    args = parser.parse_args()
    
    for arg in vars(args):
        print(f"- {arg}: {getattr(args, arg)}")
    return args


def interpolate(alpha, t1, t2 =(0.409, 0.215, 0.430, 1)):
    return tuple([t1[i]*alpha + t2[i]*(1-alpha) for i in range(len(t1))])


def uniform_sample(min, max, seed=0):
    np.random.seed(seed)
    return np.random.uniform(min, max)


def create_scene(
        tissue_thickness = 0.05, tissue_size = 1.28, tissue_location = (0, 0, 0.5),
        tissue_thickness_lb = 0.05, scale_scene = 1.15, 
        light_source_brightness = 60, adaptiv_brightness = True, tissue_color = (0.409, 0.215, 0.430, 1),
        nucleus_color = (0.315, 0.003, 0.531, 1), red_points_strength = 0, red_base = (0.605, 0.017, 0.216, 1),
        tissue_rips = -0.5, tissue_rips_std = 0.1, tissue_rips_curl = (0, 1),
        nuclei_intensity = 1, mix_cyto = 0, over_staining = (0, 1), goblet_intensity = (0.5, 2),
        tissue_padding = 0.5, epi_count = 80, stroma_density = 0.5, mix_factor = 0, stroma_intensity = 1,
        ratios = [0, 0.1, 0.8, 0.06, 0.04], focal_offset = 0, 
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

    params = {}

    # 0) parameters for variations
    base_intensity = 120*(1-nuclei_intensity)
    tissue_size = tissue_size*scale_scene
    cells.initialize_mixing_attribute(mix_factor)
    if tissue_thickness != tissue_thickness_lb:
        tissue_thickness = uniform_sample(tissue_thickness_lb, tissue_thickness, seed=seed)
        print(f"Adaptiv thickness: {tissue_thickness}")
    if adaptiv_brightness:
        light_source_brightness = (light_source_brightness)**(0.01/0.05*tissue_thickness/0.01)
        print(f"Adaptiv brightness: {light_source_brightness}")
    tissue_rips_curl = uniform_sample(tissue_rips_curl[0], tissue_rips_curl[1], seed=seed)
    over_staining = uniform_sample(over_staining[0], over_staining[1], seed=seed)
    goblet_intensity = (0.7 + uniform_sample(0, goblet_intensity[1]-goblet_intensity[0], seed=seed)*1.3)*over_staining + 2*(1-over_staining)

    # 1) initialize microscope objects and add to scene
    start = time.time()
    params_cell_shading = {
        'PLA': {
            'Nucleus': {'name': 'Nucleus_PLA', 'color': nucleus_color, 'staining_intensity': 250*nuclei_intensity+base_intensity},
            'Cytoplasm': {'name': 'Cytoplasm_PLA', 'color': tissue_color, 'staining_intensity': 230}},
        'LYM': {
            'Nucleus': {'name': 'Nucleus_LYM', 'color': interpolate(nuclei_intensity, nucleus_color), 'staining_intensity': 400*nuclei_intensity+base_intensity},},
        'EOS': {
            'Nucleus': {'name': 'Nucleus_EOS', 'color': nucleus_color, 'staining_intensity': 350*nuclei_intensity+base_intensity},
            'Cytoplasm': {'name': 'Cytoplasm_EOS', 'color': interpolate(1-mix_cyto, red_base, tissue_color), 'staining_intensity': 200}},
        'FIB': {
            'Nucleus': {'name': 'Nucleus_FIB', 'color': interpolate(nuclei_intensity, nucleus_color), 'staining_intensity': 290*nuclei_intensity+base_intensity},},
        'EPI': {
            'Nucleus': {'name': 'Nucleus_EPI', 'color': interpolate(nuclei_intensity, nucleus_color), 'staining_intensity': 290*nuclei_intensity+base_intensity}}}
    my_materials = materials.Material(
        over_staining=over_staining,
        seed=seed, cell_type_params=params_cell_shading, tissue_rips=tissue_rips, tissue_rips_curl=tissue_rips_curl, red_base=red_base,
        tissue_rips_std=tissue_rips_std, stroma_intensity=stroma_intensity, goblet_intensity=goblet_intensity, stroma_color=tissue_color,
        brightness=light_source_brightness, red_points_strength=red_points_strength)#over_staining)
    print(tissue_location)
    my_tissue = tissue.Tissue(
        my_materials.muscosa, thickness=tissue_thickness,
        size=tissue_size+tissue_padding/2, location=tissue_location)
    my_light_source = scene.LightSource(material=my_materials.light_source)
    my_camera = scene.Camera(
        focus_pos=0.6221+(tissue_thickness*(1+focal_offset)-0.05),
        size=tissue_size)
    my_scene = scene.BioMedicalScene(my_light_source, my_camera)
    my_scene.add_cell_params(params_cell_shading)
    my_scene.add_tissue(tissue=my_tissue.tissue)
    end = time.time()
    print(f"Initialization took {end - start} s")

    # 2) create macrostructures in tissue block, rotate and scale them and cut them
    start = time.time()
    tissue_arch = arch.TissueArch(seed=seed)
    elapsed = time.time() - start
    print(f"Architecture init took {elapsed} s")
    tissue_arch.random_crop(my_tissue.tissue)
    elapsed_old = elapsed
    elapsed = time.time() - start 
    print(f"Architecture crop took {elapsed-elapsed_old} s")
    macro_structure = tissue_arch.get_architecture()   # NOTE crypt is bad news, dont touch it
    for obj in macro_structure:
        obj.location[2] += 0.5  # move up to the tissue surface
    for obj in macro_structure[:-1]:
        smooth = obj.modifiers.new(name='smooth', type='SMOOTH')
        smooth.iterations = 20
        bpy.context.view_layer.objects.active = obj
        bpy.ops.object.modifier_apply(modifier=smooth.name)
    elapsed_old = elapsed
    elapsed = time.time() - start 
    print(f"Architecture get took {elapsed-elapsed_old} s")
    crypt, crypt_vol_1, crypt_vol_2, vol_goblet, mucosa = macro_structure
    mucosa_fill = hm.copy_object(mucosa, 'muscosa_fill')
    my_scene.bound_architecture(
       volumes=[mucosa_fill, crypt_vol_1, vol_goblet, crypt_vol_2], surfaces=[crypt],
       padding=tissue_padding)
    ext_stroma = hm.copy_object(mucosa_fill, 'ext_stroma')
    elapsed_old = elapsed
    elapsed = time.time() - start
    print(f"Architecture bound took {elapsed-elapsed_old} s")
    # for obj in [crypt_vol_2]:
    #     if obj.type == 'MESH':
    #         obj_name = obj.name
    #         defo.elastic_deform(obj_name, deformation_strength=0.01, noise_scale=8, seed=seed)
    hm.add_boolean_modifier(mucosa_fill, crypt_vol_2, name='add epi to stroma', apply=True)
    elapsed_old = elapsed
    elapsed = time.time() - start 
    print(f"Architecture boolean took {elapsed-elapsed_old} s")
    end = time.time()
    print(f"Architecture creation took {end - start} s")

    # 3) populate scene with nuclei/cells
    # add epi volume filling
    start = time.time()
    crypt_goblet = arr.VoronoiFill(vol_goblet, ext_stroma, cells.CellType.GOB)
    crypt_fill = arr.VoronoiFill(crypt_vol_1, ext_stroma, cells.CellType.EPI)
    my_scene.add_arrangement(crypt_fill, my_scene.tissue_empty) # NOTE: 200 nuclei take about 40 s
    my_scene.add_arrangement(crypt_goblet, my_scene.tissue_empty)
    end = time.time()
    print(f"Voronoi filling took {end - start} s")

    # Add volume filling
    # add tissue padding befor filling
    start = time.time()
    MIX_TYPES = [
        cells.CellType.MIX,
        cells.CellType.PLA, 
        cells.CellType.LYM, 
        cells.CellType.EOS, 
        cells.CellType.FIB]
    volume_fill = arr.VolumeFill(mucosa_fill, stroma_density, MIX_TYPES, ratios, seed=seed)
    my_scene.add_arrangement(volume_fill, bounding_mesh=mucosa_fill) # NOTE: 240 nuclei take about 20 s
    end = time.time()
    print(f"Volume filling took {end - start} s")

    # Deform tissue
    # apply simple subdivision to smooth the surface
    # hm.subdivide_list([mucosa_fill, crypt_vol_2], 1, type='SIMPLE')
    # hm.subdivide_object(mucosa_fill, 1, type='SIMPLE')
    # hm.subdivide_object(crypt_vol_2, 1, type='SIMPLE')
    # for obj in [mucosa_fill, crypt_vol_2]:
    #     if obj.type == 'MESH':
    #         obj_name = obj.name
    #         defo.elastic_deform(obj_name, deformation_strength=0.01, noise_scale=20, seed=seed)

    # Deform goblet cells
    for obj in bpy.data.objects:
        if obj.type == 'MESH' and obj.name.startswith("Goblet_Type_GOB"):
            obj_name = obj.name
            defo.elastic_deform(obj_name, deformation_strength=0.01, noise_scale=15, seed=seed)
            defo.elastic_deform(obj_name, deformation_strength=0.02, noise_scale=5, seed=seed)
    
    # Deform epithelial cells
    for obj in bpy.data.objects:
        if obj.type == 'MESH' and obj.name.startswith("Nucleus_Type_EPI"):
            obj_name = obj.name
            defo.elastic_deform(obj_name, deformation_strength=0.01, noise_scale=20, seed=seed)
            defo.elastic_deform(obj_name, deformation_strength=0.0025, noise_scale=100, seed=seed)

    # 4) cut objects and add staining

    # my_scene.cut_tissue()
    # elapsed_old = elapsed
    # elapsed = time.time() - start  # TODO switch back

    start = time.time()
    my_scene.add_cell_params(params_cell_shading)
    elapsed = time.time() - start
    print(f"Adding cell params took {elapsed} s")
    my_scene.delete_cells()
    elapsed_old = elapsed
    elapsed = time.time() - start
    print(f"Deleting cells took {elapsed-elapsed_old} s")
    my_scene.cut_cytoplasm_nuclei()
    elapsed_old = elapsed
    elapsed = time.time() - start
    print(f"Cutting cytoplasm and nuclei took {elapsed-elapsed_old} s")
    #my_scene.remove_goblet_volume(crypt_vol_2)
    #my_scene.cut_tissue() # TODO switch back
    my_scene.remove_cells_volume(crypt_vol_2, tolerance=0.01, types=('GOB'))
    elapsed_old = elapsed
    elapsed = time.time() - start
    print(f"Removing goblet volume took {elapsed-elapsed_old} s")
    my_scene.remove_cells_volume(mucosa_fill)
    my_scene.remove_cells_volume(crypt_vol_2, types=('EPI'))
    elapsed_old = elapsed
    elapsed = time.time() - start
    print(f"Removing cells volume took {elapsed-elapsed_old} s")
    my_scene.cut_cells()
    elapsed_old = elapsed
    elapsed = time.time() - start
    print(f"Cutting cells took {elapsed-elapsed_old} s")
    my_scene.cut_tissue() # TODO switch back
    elapsed_old = elapsed
    elapsed = time.time() - start
    print(f"Cutting tissue took {elapsed-elapsed_old} s")
    my_scene.add_tissue_staining(materials=[my_materials.muscosa, my_materials.crypt_staining])
    elapsed_old = elapsed
    elapsed = time.time() - start
    print(f"Adding tissue staining took {elapsed-elapsed_old} s")
    my_scene.add_nuclei_mask(material=my_materials.nuclei_mask)
    elapsed_old = elapsed
    elapsed = time.time() - start
    print(f"Adding nuclei mask took {elapsed-elapsed_old} s")
    my_scene.add_staining_to_cell(materials=my_materials.cell_staining)
    my_scene.add_staining_to_cell(materials=[my_materials.goblet_staining])
    elapsed_old = elapsed
    elapsed = time.time() - start
    print(f"Adding cell staining took {elapsed-elapsed_old} s")
    mucosa_fill.location.z = mucosa_fill.location.z - 0.0005
    end = time.time()
    print(f"Cutting and staining took {end - start} s")

    # 6) add elastic deformations
    start = time.time()
    for obj in bpy.data.objects:
        # Check if the object is a cube (type is 'MESH' and name starts with 'Cube')
        print(obj.name)
        if obj.type == 'MESH' and obj.name.startswith("Plane")==False and obj.name.startswith("tissue")==False:
            obj_name = obj.name
            print(obj.name)
            defo.elastic_deform(obj_name, seed=seed, deformation_strength=0.0125)
    end = time.time()
    print(f"Deforming objects took {end - start} s")

    # 5) hide non cell objects
    start = time.time()
    goblet_cells = []
    for cell in my_scene.cell_objects:
        cell_type = cell.name.split('_')[-2]
        if cell_type == 'GOB':
            goblet_cells.append(cell)
    for obj in [crypt, crypt_vol_1, mucosa, ext_stroma, vol_goblet]:#+goblet_cells:
        obj.hide_viewport = True
        obj.hide_render = True
    end = time.time()
    print(f"Hiding non cell objects took {end - start} s")

    # write all randomized parameters to a file
    params = {
        'seed': seed,
        'tissue_thickness': tissue_thickness,
        'tissue_size': tissue_size,
        'tissue_location': tissue_location,
        'tissue_padding': tissue_padding,
        'tissue_rips': tissue_rips,
        'tissue_rips_std': tissue_rips_std,
        'tissue_rips_curl': tissue_rips_curl,
        'stroma_intensity': stroma_intensity,
        'nuclei_intensity': nuclei_intensity,
        'mix_factor': mix_factor,
        'stroma_density': stroma_density,
        'ratios': ratios,
        'focal_offset': focal_offset,
        'light_source_brightness': light_source_brightness,
        'adaptiv_brightness': adaptiv_brightness,
        'over_staining': over_staining,
        'red_points_strength': red_points_strength}

    return my_scene, params


class MainpulateScene(object):
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

    def inflate_epi_cells(self, factor):
        for cell in self.my_scene.cell_objects:
            cell_type = cell.name.split('_')[-2]
            if cell_type == 'EPI':
                cell.scale = tuple([cell.scale[i]*factor for i in range(3)])



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


def render_scene(
        my_scene, render_path, sample_name, gpu=True, device=0, 
        output_shape=(512, 512), max_samples=1024, render_masks=True, base_16bit=55,
        additional_info=None):
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
    bpy.context.scene.cycles.use_denoising = False
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
                    single_masks = render_masks, # if true single cell masks are rendered
                    semantic_mask = render_masks, # if true semantic mask is generated
                    instance_mask = render_masks, # if true instance mask is generated
                    cyto_mask = render_masks, # if true cytoplasm mask is generated
                    depth_mask = False, # if true depth mask is generated
                    obj3d = False, # if true scene is saved as 3d object
                    output_shape = output_shape, # dimensions of output
                    max_samples = max_samples,
                    base_16bit = base_16bit,
                    additional_info=additional_info) # number of samples for rendering. Fewer samples will render more quickly. Default is 1024


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

    # if index list is given, render only these indices
    if len(args.index_list) > 1:
        with open(render_path + '/' + args.index_list, 'r') as f:
            indices = f.read().splitlines()
        indices = [int(i) for i in indices]
        indices = [i-1 for i in indices]  # NOTE tempory fix, due to 1-indexing
        indices = indices[args.start_idx: args.start_idx + args.n_samples]
    else:
        indices = list(range(args.start_idx, args.start_idx + args.n_samples))
    print(f'Indices: {indices}')

    # render individual samples
    for i in tqdm(indices):
        paramters = {'seed': i}
        for key, value in args.__dict__.items():
            if key not in paramters.keys():
                paramters[key] = value
        with open(dir_parameters+f'/parameters_{i+1}.json', 'w') as outfile:
            json.dump(paramters, outfile)
        my_scene, rand_params = create_scene(**paramters)
        render_scene(
            my_scene, render_path, i+1, gpu=args.gpu, device=args.gpu_device,
            base_16bit=args.base_16bit, additional_info=rand_params)
        bpy.ops.wm.read_factory_settings(use_empty=True)

if __name__ == "__main__":
    main()