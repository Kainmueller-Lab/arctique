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
import src.shading.shading as shading
import src.scene as scene
import src.utils as utils
import src.utils.geometry as geom
import src.utils.surface_filling as sf

# this next part forces a reload in case you edit the source after you first start the blender session
#import imp
import importlib as imp # imp module is deprecated since python 3.12
imp.reload(arr)
imp.reload(cells)
imp.reload(tissue)
imp.reload(shading)
imp.reload(scene)
imp.reload(utils)
imp.reload(geom)
imp.reload(sf)



def parse_dataset_args():
    parser = argparse.ArgumentParser()
    
    # RENDERING PARAMETERS
    parser.add_argument("--gpu_device")
    parser.add_argument("--output_dir", type=str, default="J:/jannik/GitHub/rendered_HE/rendered", help="Set output folder")
    parser.add_argument("--n_samples", type=int, default=10, help="Dataset size")
   
    # DATASET PARAMETERS
    # tissue
    parser.add_argument("--tissue_thickness", type=float, default=0.2, help="Tissue thickness")
    parser.add_argument("--tissue_size", type=float, default=2, help="Tissue size")
    parser.add_argument("--tissue_location", type=tuple, default=(0, 0, 0.5), help="Tissue location")
    parser.add_argument("--tissue_padding", type=float, default=0.5, help="Tissue padding")
    
    # nuclei
    parser.add_argument("--surf_number", type=int, default=80, help="number of surface cells")
    parser.add_argument("--filler_scale", type=float, default=0.8, help="Scale of the size of smaller filler nuclei w.r.t to the original nuclei size")
    parser.add_argument("--number", type=int, default=80, help="number of volume cells")
    parser.add_argument("--ratios", type=list, default=[0.6, 0.2, 0.2], help="ratios of different cell types")
    parser.add_argument("--vol_scale", type=tuple, default=(1, 0.7, 1), help="Volume scale")
    parser.add_argument("--surf_scale", type=tuple, default=(0.8, 0.5, 1), help="Surface scale")

    #other default value for --output_dir: "/Volumes/ag_kainmueller/vguarin/synthetic_HE" via internal VPN
    
    args = parser.parse_args()
    
    for arg in vars(args):
        print(f"- {arg}: {getattr(args, arg)}")
    return args


def create_scene(
        tissue_thickness = 0.2, tissue_size = 2, tissue_location = (0, 0, 0.5),
        tissue_padding = 0.5, surf_number = 80, filler_scale = 0.8, number = 80, 
        ratios = [0.6, 0.2, 0.2], vol_scale = (1, 0.7, 1), surf_scale = (0.8, 0.5, 1)):
    '''
    creates a tissue crop with cells and nuclei
    Args:
        tissue_thickness: float, thickness of the tissue [unit: 10 micrometers]
        tissue_size: float, size of the tissue [unit: 10 micrometers]
        tissue_location: tuple, location of the tissue [unit: 10 micrometers]
        tissue_padding: float, padding of the tissue [unit: 10 micrometers]
        surf_number: int, number of surface cells
        filler_scale: float, scale of the size of smaller filler nuclei w.r.t to the original nuclei size
        number: int, number of volume cells
        ratios: list, ratios of different cell types
        vol_scale: tuple, volume scale
        surf_scale: tuple, surface scale
    Returns:
        my_scene: BioMedicalScene object
    '''
    scene.BioMedicalScene.clear()
    
    # 1) initialize microscope objects and add to scene
    my_materials = shading.Material()
    my_tissue = tissue.Tissue(my_materials.tissue_staining, thickness=tissue_thickness, size=tissue_size, location=tissue_location) # thickness and location of tissue should encapsulate min and max z-coordinates of cells 
    my_light_source = scene.LightSource(material=my_materials.light_source)
    my_camera = scene.Camera()
    my_scene = scene.BioMedicalScene(my_light_source, my_camera)

    # 2) create macrostructures in tissue block, rotate and scale them and cut them
    VOL_OBJ, SURF_OBJ = my_scene.add_dummy_objects(my_tissue, tissue_padding, vol_scale, surf_scale)
    my_scene.add_tissue(tissue=my_tissue.tissue)

    # 3) populate scene with nuclei/cells
    # Add surface filling -> crypts
    SURF_ATTRIBUTE = cells.CellAttributeEpi()
    surface_fill = arr.SurfaceFill(SURF_OBJ, surf_number, SURF_ATTRIBUTE, filler_scale)  # NOTE: For some very weird reason you need to create the surface filling before the volume filling. TODO: Fix that
    my_scene.add_arrangement(surface_fill)

    # Add volume filling
    ATTRIBUTES = [cells.CellAttributeA(), cells.CellAttributeB(), cells.CellAttributeC()]
    volume_fill = arr.VolumeFill(VOL_OBJ, number, ATTRIBUTES, ratios, strict_boundary=False)
    my_scene.add_arrangement(volume_fill)

    # 4) cut objects and add staining
    my_scene.cut_cells()
    my_scene.add_staining(material=my_materials.nuclei_staining)

    # 5) hide non cell objects
    VOL_OBJ.hide_viewport = True
    VOL_OBJ.hide_render = True
    SURF_OBJ.hide_viewport = True
    SURF_OBJ.hide_render = True

    return my_scene


def render_scene(my_scene, render_path, sample_name, output_shape=(500, 500), max_samples=1024):
    
    # set render engine
    bpy.context.scene.render.engine = 'CYCLES'
    bpy.context.scene.cycles.device = 'GPU'
    bpy.context.preferences.addons['cycles'].preferences.compute_device_type = "CUDA"
    bpy.context.preferences.addons["cycles"].preferences.get_devices()
    print(bpy.context.preferences.addons["cycles"].preferences.compute_device_type)
    for j, d in enumerate(bpy.context.preferences.addons["cycles"].preferences.devices):
        if j in (0, 1): # Using only the first two devices
            d["use"] = True
        else:
            d["use"] = False
        print(d["name"], d["use"])
    
    my_scene.sample_name = sample_name
    
    my_scene.render(filepath = render_path,  # where to save renders
                    scene = True, # if true scene is rendered
                    single_masks = False, # if true single cell masks are rendered
                    semantic_mask = True, # if true semantic mask is generated
                    instance_mask = True, # if true instance mask is generated
                    depth_mask = False, # if true depth mask is generated
                    obj3d = False, # if true scene is saved as 3d object
                    output_shape = output_shape, # dimensions of output
                    max_samples = max_samples) # number of samples for rendering. Fewer samples will render more quickly. Default is 1024


def main():
    args = parse_dataset_args() 
    
    # create output directory
    render_path = args.output_dir #render scene
    print(render_path) 
    dir = render_path + '/train_combined_masks/semantic'
    if not os.path.exists(dir):
            os.makedirs(dir)
    max_n_samples = len(glob.glob(f"{dir}/*.png")) 

    # render individual samples
    for i in tqdm(range(max_n_samples, max_n_samples + args.n_samples)):
        my_scene = create_scene(
            tissue_thickness = args.tissue_thickness, tissue_size = args.tissue_size, 
            tissue_location = args.tissue_location, tissue_padding = args.tissue_padding,
            surf_number = args.surf_number, filler_scale = args.filler_scale, number = args.number, 
            ratios = args.ratios, vol_scale = args.vol_scale, surf_scale = args.surf_scale)
        render_scene(my_scene, render_path, i+1)
        bpy.ops.wm.read_factory_settings(use_empty=True)


if __name__ == "__main__":
    main()