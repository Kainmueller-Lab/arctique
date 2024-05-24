import bpy
import sys
import os
import argparse
import glob
from tqdm import tqdm

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



def parse_dataset_args():
    parser = argparse.ArgumentParser()
    
    # RENDERING PARAMETERS                                                                                                                                                        # add argument with list of all gpu devices
    parser.add_argument("--gpu_devices", type=list, default=[1], help="List of GPU devices to use for rendering")
    parser.add_argument("--gpu", type=bool, default=True, help="Use GPU for rendering")
    parser.add_argument("--output_dir", type=str, default="rendered", help="Set output folder")
    parser.add_argument("--start_idx", type=int, default=0, help="Dataset size")
    parser.add_argument("--n_samples", type=int, default=500, help="Dataset size")

    # DATASET PARAMETERS
    # tissue
    parser.add_argument("--tissue_thickness", type=float, default=0.05, help="Tissue thickness")
    parser.add_argument("--tissue_size", type=float, default=2, help="Tissue size")
    parser.add_argument("--tissue_location", type=tuple, default=(0, 0, 0.5), help="Tissue location")
    parser.add_argument("--tissue_padding", type=float, default=1, help="Tissue padding")
    
    # nuclei
    parser.add_argument("--surf_number", type=int, default=80, help="number of surface cells")
    parser.add_argument("--filler_scale", type=float, default=0.8, help="Scale of the size of smaller filler nuclei w.r.t to the original nuclei size")
    parser.add_argument("--number", type=int, default=800, help="number of volume cells")
    parser.add_argument("--ratios", type=list, default=[0.1, 0.3, 0.4, 0.1, 0.1], help="ratios of different cell types")
    parser.add_argument("--surf_scale", type=tuple, default=(0.8, 0.5, 1), help="Surface scale")

    #other default value for --output_dir: "/Volumes/ag_kainmueller/vguarin/synthetic_HE" via internal VPN
    
    args = parser.parse_args()
    
    for arg in vars(args):
        print(f"- {arg}: {getattr(args, arg)}")
    return args


def create_scene(
        tissue_thickness = 0.05, tissue_size = 2, tissue_location = (0, 0, 0.5),
        tissue_padding = 0.5, surf_number = 80, number = 80, 
        ratios = [0.1, 0.3, 0.4, 0.1, 0.1],
        seed=0):
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
    print(number)
    scene.BioMedicalScene.clear()
        
    # 1) initialize microscope objects and add to scene
    my_materials = materials.Material(seed=seed)
    my_tissue = tissue.Tissue(
        my_materials.muscosa, thickness=tissue_thickness,
        size=tissue_size, location=tissue_location)
    my_light_source = scene.LightSource(material=my_materials.light_source)
    my_camera = scene.Camera()
    my_scene = scene.BioMedicalScene(my_light_source, my_camera)
    my_scene.add_tissue(tissue=my_tissue.tissue)

    # 2) create macrostructures in tissue block, rotate and scale them and cut them
    tissue_arch = arch.TissueArch(seed=seed)
    tissue_arch.random_crop(my_tissue.tissue)
    macro_structure = tissue_arch.get_architecture()
    crypt, crypt_vol_1, crypt_vol_2, mucosa = macro_structure
    my_scene.bound_architecture(
        volumes=[crypt_vol_1, crypt_vol_2, mucosa], surfaces=[crypt],
        padding=tissue_padding)

    # 3) populate scene with nuclei/cells
    # add bounding volumes
    # NOTE: MIX_VOL bounds the volume for the mixed cell types
    # NOTE: EPI_VOL bounds the volume for the epithelial cell types.
    # MIX_VOL, EPI_VOL = utils.geometry.add_dummy_volumes(my_tissue, tissue_padding)

    # # add epi volume filling
    # EPI_COUNT = 200
    # EPI_ATTRIBUTE = cells.CellAttributeEpi(size=0.1, scale=(1, 0.5, 0.5))
    # crypt_fill = arr.VoronoiFill(EPI_VOL, EPI_COUNT, EPI_ATTRIBUTE)
    # my_scene.add_arrangement(crypt_fill) # NOTE: 200 nuclei take about 30 s

    # Add volume filling
    # add tissue padding befor filling
    MIX_TYPES = [
        cells.CellType.MIX,
        cells.CellType.PLA, 
        cells.CellType.LYM, 
        cells.CellType.EOS, 
        cells.CellType.FIB]
    volume_fill = arr.VolumeFill(
        mucosa, number, MIX_TYPES, ratios, strict_boundary=False, seed=seed)
    my_scene.add_arrangement(volume_fill)
    my_scene.cut_cells(boolean_object=mucosa)
    #volume_fill = arr.VolumeFill(MIX_VOL, MIX_COUNT, MIX_TYPES, RATIOS, strict_boundary=True)

    # 4) cut objects and add staining
    my_scene.cut_cells()
    my_scene.cut_tissue()
    my_scene.add_tissue_staining(materials=[my_materials.muscosa, my_materials.crypt_staining])
    my_scene.add_staining(material=my_materials.nuclei_mask)
    my_scene.add_staining(material=my_materials.nuclei_staining)

    # 5) hide non cell objects
    for obj in [crypt, crypt_vol_1]:
        obj.hide_viewport = True
        obj.hide_render = True

    return my_scene


def recreate_scene(parameters):
    '''
    recreates a scene from parameters
    Args:
        parameters: dict, parameters of the scene
    Returns:
        my_scene: BioMedicalScene object
    '''
    my_scene = create_scene(
        tissue_thickness = parameters['tissue_thickness'], tissue_size = parameters['tissue_size'], 
        tissue_location = parameters['tissue_location'], tissue_padding = parameters['tissue_padding'],
        surf_number = parameters['surf_number'], filler_scale = parameters['filler_scale'], number = parameters['number'], 
        ratios = parameters['ratios'], vol_scale = parameters['vol_scale'], surf_scale = parameters['surf_scale'])
    return my_scene


def render_scene(my_scene, render_path, sample_name, gpu=True, devices=[0], output_shape=(512, 512), max_samples=1024):
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
        if j in devices:
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
    if not os.path.exists(dir):
            os.makedirs(dir)

    # render individual samples
    for i in tqdm(range(args.start_idx, args.start_idx + args.n_samples)):
        my_scene = create_scene(
            tissue_thickness = args.tissue_thickness, tissue_size = args.tissue_size, 
            tissue_location = args.tissue_location, tissue_padding = args.tissue_padding,
            surf_number = args.surf_number, number = args.number, 
            ratios = args.ratios,
            seed=i)
        render_scene(my_scene, render_path, i+1, gpu=args.gpu, devices=args.gpu_devices)
        bpy.ops.wm.read_factory_settings(use_empty=True)


if __name__ == "__main__":
    main()