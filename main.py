import bpy
import numpy as np
import random
import sys
import os
from math import radians, sin, cos, pi
from mathutils import Matrix, Vector
from pathlib import Path

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
import imp
imp.reload(arr)
imp.reload(cells)
imp.reload(tissue)
imp.reload(shading)
imp.reload(scene)

###################  PARAMETER  #####################
# args_camera = {'pos'}


###################  MAIN  METHOD  #####################
# create the necessary objects
scene.BioMedicalScene.clear()
    
# add microscope objects
my_materials = shading.Material()
my_tissue = tissue.Tissue(my_materials.tissue_staining)
my_light_source = scene.LightSource(material=my_materials.light_source)
my_camera = scene.Camera()

# create scene
my_scene = scene.BioMedicalScene(my_light_source, my_camera)

# define cell arrangements
cell_distribution = arr.CellDistribution(
    cell_attributes = cells.CellAttributeA(),
    num_cells = 200,
    min_coords = Vector([-1, -1, 0.4]),
    max_coords = Vector([1, 1, 0.6])
)
locations = [Vector([-0.5, -0.5, 0.5]), Vector([0.5, -0.5, 0.5]), Vector([-0.5, 0.5, 0.5]), Vector([0.5, 0.5, 0.5])]
cell_list = arr.CellList(
    cell_attributes = cells.CellAttributeB(),
    locations = locations
)

# add cell arrangements to scene
my_scene.add_arangement(cell_distribution)
my_scene.add_arangement(cell_list)


# render scene
my_scene.render()



# Setup a folder called 3d_outputs and export scene as obj 
current_folder = os.path.dirname(os.path.realpath(__file__))
FOLDER = Path(current_folder).joinpath("Images")
FOLDER = str(FOLDER)
try:
    if not os.path.exists(FOLDER):
        os.makedirs(FOLDER)
except OSError as error:
    print("Directory '%s' can not be created")


bpy.ops.export_scene.obj(filepath=FOLDER+"//my_scene.obj")

