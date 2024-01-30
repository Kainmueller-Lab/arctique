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
my_scene.cut_cells()
my_scene.add_staining(material=my_materials.nuclei_staining)

RENDER_PATH = 'C:/Users/cwinklm/Documents/Alpacathon/rendered_HE/renders2/'

my_scene.render(filepath = RENDER_PATH,  # where to save renders
               scene = True, # if true scene is rendered
               single_masks = True, # if true singel cell masks are rendered
               semantic_mask = True, # if true semantic mask is generated
               instance_mask = True, # if true instance mask is generated
               depth_mask = True, # if true depth mask is generated
               obj3d = True, # if true scene is saved as 3d object
               output_shape = (500, 500), # dimensions of output
               max_samples = 10) # number of samples for rendering. Fewer samples will render more quickly. Default is 1024








