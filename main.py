import bpy
import numpy as np
import random
import sys
import os
from math import radians, sin, cos, pi
from mathutils import Matrix, Vector, geometry
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
# my_scene.add_arrangement(cell_distribution_A)
# my_scene.add_arrangement(cell_distribution_B)



# OPTION 2: Voronoi cell distributions ############################
# Given a list of 3D points (could be a randomly created list or a deterministic list of positions) per cell attribute
# a Voronoi diagram is created with these points as seeds.
# In each Voronoi region a nucleus object is created with size and scale corresponding to the cell attribute.
# The Voronoi diagram ensures that the objects to not intersect.
# TODO: Add bending if necessary

# Create 3D point lists per cell attribute
# cell_count_A = 60
# cell_count_B = 40
# min_coords = Vector([-1, -1, 0.4])
# max_coords = Vector([1, 1, 0.6])
# points_A = [list(map(random.uniform, min_coords, max_coords)) for _ in range(cell_count_A)]
# points_B = [list(map(random.uniform, min_coords, max_coords)) for _ in range(cell_count_B)]
# # Store the point lists per attribute in a dict
# distribution_dict = {}
# distribution_dict[cells.CellAttributeA()] = points_A
# distribution_dict[cells.CellAttributeB()] = points_B 
# # define cell arrangements and add to scene
# voronoi_arr = arr.VoronoiDiagram(distribution_dict)
# my_scene.add_arrangement(voronoi_arr)


# OPTION 3: Epithelial crypts ############################
# - This script models tissue cuts through epithelial crypts.
# - The nuclei are created from a Voronoi diagram whose seeds are placed onto the mesh of a scaled icosphere.
# - The seed points are not chosen randomly but set as the vertices and face centers of the mesh.
# - The resulting Voronoi regions are intersected with an inner and outer icosphere to create small nuclei objects along the surface.
# - In order to reduce the computing time not the whole surface and nuclei are rendered but only a slice of changeable thickness.
# - The resulting objects form a cut though such a crypt with nuclei at its boundary.
# - In addition two hull objects are created that can be used to define the inside and outside of the crypt.
# - It is possible to manually set the scale, location and rotation in xy-plane of the resulting crypt cut.
# NOTE: Generating 4 crypts currently take about 1-2 minutes

### PARAMETERS
# TODO: Create parameters (max 4?) based on lattice
ico_scales = [(0.5, 0.22), (0.4, 0.17), (0.45, 0.2), (0.45, 0.15)]
angles = [40, 60, 70, 70]
centers = [(-0.5,-0.5,0.5), (0.5,0.5,0.5), (-0.5,0.5,0.5), (0.5,-0.5,0.5)]

outer_hulls = []
count = 0
for ico_scale, angle, center in zip(ico_scales, angles, centers):
    param_dict = {}
    param_dict["ico_xy_scale"] = ico_scale # Scale of the icosphere w.r.t. to the x-y-axes.
    param_dict["z_rot_angle"] = angle # Rotation along z-axis of crypt in degrees.
    param_dict["center_loc"] = center # Center of crypt cut in world coordinates.
    # Define cell arrangements and add to scene
    epi_arr = arr.EpithelialArrangement(param_dict)
    # TODO: Hide hulls and extract them as objects to feed into distribution arrangement
    my_scene.add_arrangement(epi_arr)
    outer_hulls.append(epi_arr.outer_hull)
    count += 1
    if count>2:
        break

# Add nuclei distribution
# Create 3D point lists per cell attribute
cell_count_A = 200
min_coords = Vector([-1, -1, 0.45])
max_coords = Vector([1, 1, 0.55])
points_A = [list(map(random.uniform, min_coords, max_coords)) for _ in range(cell_count_A)]
distribution_dict = {}
distribution_dict[cells.CellAttributeA()] = points_A
voronoi_arr = arr.VoronoiDiagram(distribution_dict)
# Define which regon should not be populated by the distribution
voronoi_arr.empty_regions = outer_hulls
my_scene.add_arrangement(voronoi_arr)

# TODO: Name the nuclei in a reasonable way
# TODO: Refactor the add method of EpithelialArr.
# TODO: Refactor the optional examples in main.py

# my_scene.add_tissue(tissue=my_tissue.tissue)
# my_scene.cut_cells()
# my_scene.add_staining(material=my_materials.nuclei_staining)

my_scene.hide_auxiliary_objects()

# # render scene
# my_scene.render(filepath='renders/')



# # Setup a folder called 3d_outputs and export scene as obj 
# # current_folder = os.path.dirname(os.path.realpath(__file__))
# FOLDER = Path(dir).joinpath("Images")#Path(current_folder).joinpath("Images")
# FOLDER = str(FOLDER)
# try:
#     if not os.path.exists(FOLDER):
#         os.makedirs(FOLDER)
# except OSError as error:
#     print("Directory '%s' can not be created")


# bpy.ops.export_scene.obj(filepath=FOLDER+"//my_scene.obj")
