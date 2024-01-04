import bpy
import numpy as np
import random
import sys
import os
from math import radians, sin, cos, pi
from mathutils import Matrix, Vector
from pathlib import Path
from PIL import Image

# IMPORT SOURCES
dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir )

import src.arrangement.arrangement as arr 
import src.objects.cells as cells
import src.objects.tissue as tissue
import src.shading.shading as shading
import src.scene as scene


# this next part forces a reload in case you edit the source after you first start the blender session
#import imp
import importlib as imp # imp module is deprecated since python 3.12
imp.reload(arr)
imp.reload(cells)
imp.reload(tissue)
imp.reload(shading)
imp.reload(scene)


###################  PARAMETER  #####################
# args_camera = {'pos'} # no change just test


###################  MAIN  METHOD  #####################
# create the necessary objects
scene.BioMedicalScene.clear()
    
# add microscope objects
my_materials = shading.Material()
my_tissue = tissue.Tissue(my_materials.tissue_staining, thickness=0.2, size=2, location=(0, 0, 0.5)) # thickness and location of tissue should encapsulate min and max z-coordinates of cells 
my_light_source = scene.LightSource(material=my_materials.light_source)
my_camera = scene.Camera()

# create scene
my_scene = scene.BioMedicalScene(my_light_source, my_camera)

# define cell arrangements
cell_distribution_A = arr.CellDistribution(
    cell_attributes = cells.CellAttributeA(),
    num_cells = 10,
    #num_cells = 100,
    min_coords = Vector([-1, -1, 0.4]),
    max_coords = Vector([1, 1, 0.6])
)
cell_distribution_B = arr.CellDistribution(
    cell_attributes = cells.CellAttributeB(),
    num_cells = 3,
    #num_cells = 30,
    min_coords = Vector([-1, -1, 0.4]),
    max_coords = Vector([1, 1, 0.6])
)

# add cell arrangements to scene
my_scene.add_arangement(cell_distribution_A)
my_scene.add_arangement(cell_distribution_B)
my_scene.add_tissue(tissue=my_tissue.tissue)


# loop through slices: 
slice_thickness = 0.05
my_scene.tissue_empty.scale.z = slice_thickness

RENDER_PATH = 'C:/Users/cwinklm/Documents/Alpacathon/rendered_HE/renders3d/'


for idx, loc in enumerate(np.arange(0.4, 0.6, slice_thickness)): 
    my_scene.tissue_empty.location.z = loc
    my_scene.cut_cells()
    my_scene.add_staining(material=my_materials.nuclei_staining)

    my_scene.render3d(filepath = RENDER_PATH,  # where to save renders
                scene = False, 
                semantic_mask = True, 
                instance_mask = False, 
                obj3d = False,     
                depth_mask = True, 
                output_shape = (500, 500), 
                semantic_mask_name = f"Semantic_mask_{idx}") 
    my_scene.uncut_cells()#



frames = []
for idx in range(len(np.arange(0.4, 0.6, slice_thickness))): 
    frames.append(Image.open(Path(RENDER_PATH).joinpath(f"Semantic_mask_{idx}.png")))

frame_one = frames[0]
frame_one.save(Path(RENDER_PATH).joinpath("download.gif"),
               format="GIF", 
               append_images=frames[1:], #list of images to append as additional frames.
               save_all=True, 
               duration=500, #display duration of each frame, in milliseconds
               loop=0)#Number of times to repeat the animation'





# def make_gif(frame_folder):
#     frames = [Image.open(image) for image in glob.glob(f"{frame_folder}/*")]
#     frame_one = frames[0]
#     frame_one.save(Path(RENDER_PATH).joinpath("download.gif"), 
#                    format="GIF", 
#                    append_images=frames,
#                    save_all=True, 
#                    duration=1000, 
#                    loop=1)











##### remove and apply modifiers: 
##https://blender.stackexchange.com/questions/44514/how-to-remove-apply-all-modifiers-of-one-object-in-python

# check depth of field option in render settings 
    

