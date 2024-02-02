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
# TODO: Shorten this:
# - The nuclei are created based on Voronoi regions intersected with a box,
# i.e. those base meshes have all the same height and have planar boundary faces on top
# and bottom -> this is too synthetic/regular
# - The size of nuclei can be controlled using REGION_SCALE and NUCLEI_SCALE
# - However it is not possible to set an absoulte scale to the NUCLEI_SCALE
# This might lead to issues when we need specific differing scales between nuclei of different types
# - The Voronoi diagram still highly depends on the number and location of "boundary" seeds
# We need to fix a sensible number and structure for that
# - One option to get rid of too regular boundary nuclei is to set the
# bounding distributio box larger than necessary. This way, when
# cropping to camera we only have visible "real inner" nuclei.
# - Find out which of the auxiliary points (cube, octogon, lattice) is best suited
# So far it seems they affect the digram structure but not substantially the regions after intersection.
# - Add padding to the distribution box, which gets cut away after intersecting. This leads to "real inner" nuclei.

### PARAMETERS
# TODO: Create parameters (max 4?) based on lattice
# NOTE: Generating 4 crypts currently take about 1-2 minutes
ico_scales = [(0.5, 0.25), (0.5, 0.2), (0.6, 0.3), (0.45, 0.15)]
angles = [40, 60, 50, 70]
centers = [(-0.5,-0.5,0.5), (0.5,0.5,0.5), (-0.5,0.5,0.5), (0.5,-0.5,0.5)]

# outer_hulls = []
for ico_scale, angle, center in zip(ico_scales, angles, centers):
    param_dict = {}
    param_dict["ico_xy_scale"] = ico_scale # Scale of the icosphere w.r.t. to the x-y-axes.
    param_dict["z_rot_angle"] = angle # Rotation along z-axis of crypt in degrees.
    param_dict["center_loc"] = center # Center of crypt cut in world coordinates.
    # Define cell arrangements and add to scene
    epi_arr = arr.EpithelialArrangement(param_dict)
    # TODO: Hide hulls and extract them as objects to feed into distribution arrangement
    #outer_hull = ...
    #outer_hulls.append(outer_hull)
    my_scene.add_arrangement(epi_arr)

# TODO: Create distribution of other cell type nuclei that are only generated outside the outer hulls.
# TODO: Name the nuclei in a reasonable way
# TODO: Refactor the add method of EpithelialArr.

# my_scene.add_tissue(tissue=my_tissue.tissue)
# my_scene.cut_cells()
# my_scene.add_staining(material=my_materials.nuclei_staining)



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
