import bpy
import random
import sys
import os
from mathutils import Vector

# IMPORT SOURCES
dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir )

import src.arrangement.arrangement as arr 
import src.objects.cells as cells
import src.objects.tissue as tissue
import src.shading.shading as shading
import src.scene as scene
import src.utils as utils
from src.utils.helper_methods import generate_lattice_parameters

# this next part forces a reload in case you edit the source after you first start the blender session
#import imp
import importlib as imp # imp module is deprecated since python 3.12
imp.reload(arr)
imp.reload(cells)
imp.reload(tissue)
imp.reload(shading)
imp.reload(scene)
imp.reload(utils)

###################  PARAMETER  #####################
# args_camera = {'pos'} # no change just test


NUMBER = 200
TYPES = ["A", "B", "C"]
RADII = [0.1, 0.05, 0.03]
RATIOS = [0.1, 0.3, 0.6]

# Example torus mesh
bpy.ops.mesh.primitive_torus_add()
MESH = bpy.context.active_object




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


# TODO
# - Add different sizes
# - Add different scales
# - Add differen orientations
# - Add blowup algorithm (Monte Carlo)

# add cell arrangement
# TODO: deal with attributes: size, scale, typename, deformation etc.
volume_fill = arr.VolumeFill(MESH, NUMBER, ATTRIBUTES, RATIOS)
my_scene.add_arrangement(volume_fill)


# Add tissue
my_scene.add_tissue(tissue=my_tissue.tissue)
my_scene.cut_cells()
my_scene.add_staining(material=my_materials.nuclei_staining)

# render scene
RENDER_PATH = 'C:/Users/cwinklm/Documents/Alpacathon/rendered_HE/renders/'
#RENDER_PATH = 'renders/'

my_scene.render(filepath = RENDER_PATH,  # where to save renders
               scene = True, # if true scene is rendered
               masks = True, # if true singel cell masks are rendered
               semantic_mask = True, # if true semantic mask is generated
               instance_mask = True, # if true instance mask is generated
               depth_mask = True, # if true depth mask is generated
               obj3d = True, # if true scene is saved as 3d object
               remove_single_masks = False, # if False single cell masks are deleted 
               output_shape = (500, 500), # dimensions of output
               max_samples = 10) # number of samples for rendering. Fewer samples will render more quickly. Default is 1024








