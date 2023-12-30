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

# OPTION 1: Distributions & Deformations ############################
# define cell arrangements
# TODO: Add use_unstable_deform flag to toggle deformation mode
# cell_distribution_A = arr.CellDistribution(
#     cell_attributes = cells.CellAttributeA(),
#     num_cells = 50,
#     min_coords = Vector([-1, -1, 0.4]),
#     max_coords = Vector([1, 1, 0.6])
# )
# cell_distribution_B = arr.CellDistribution(
#     cell_attributes = cells.CellAttributeB(),
#     num_cells = 10,
#     min_coords = Vector([-1, -1, 0.4]),
#     max_coords = Vector([1, 1, 0.6])
# )
# # add cell arrangements to scene
# my_scene.add_arangement(cell_distribution_A)
# my_scene.add_arangement(cell_distribution_B)

# OPTION 2: Voronoi diagram ############################
# NOTE: Per cell attribute type set the number of cells you want.
# Those will be distributed rndomly in a box given by min/max_coord.
# TODO: Add cell type to nucleus object name
# TODO: Add type specific scale to the nuclei

# Create 3D point lists per cell attribute
# NOTE: you can pass deterministic lists or sample points randomly
cell_count_A = 60
cell_count_B = 40
min_coords = Vector([-1, -1, 0.4])
max_coords = Vector([1, 1, 0.6])
points_A = [list(map(random.uniform, min_coords, max_coords)) for _ in range(cell_count_A)]
points_B = [list(map(random.uniform, min_coords, max_coords)) for _ in range(cell_count_B)]
distribution_dict = {}
distribution_dict[cells.CellAttributeA()] = points_A
distribution_dict[cells.CellAttributeB()] = points_B 
#padding = (max_coords - min_coords) * 0.1
padding = Vector([0,0,0])
# define cell arrangements
voronoi_arr = arr.VoronoiDiagram(distribution_dict, min_coords-padding, max_coords+padding)
# add cell arrangements to scene
my_scene.add_arangement(voronoi_arr)

my_scene.add_tissue(tissue=my_tissue.tissue)
my_scene.cut_cells()
my_scene.add_staining(material=my_materials.nuclei_staining)



# render scene
my_scene.render(filepath='renders/')



# Setup a folder called 3d_outputs and export scene as obj 
# current_folder = os.path.dirname(os.path.realpath(__file__))
FOLDER = Path(dir).joinpath("Images")#Path(current_folder).joinpath("Images")
FOLDER = str(FOLDER)
try:
    if not os.path.exists(FOLDER):
        os.makedirs(FOLDER)
except OSError as error:
    print("Directory '%s' can not be created")


bpy.ops.export_scene.obj(filepath=FOLDER+"//my_scene.obj")
