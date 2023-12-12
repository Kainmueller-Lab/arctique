import bpy
import numpy as np
import random
import sys
import os
import bpy
from math import radians, sin, cos, pi
from mathutils import Matrix, Vector
from pathlib import Path

# IMPORT SOURCES
dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir )

import src.helper_methods as hm
import src.cell_classes as cell_classes

# this next part forces a reload in case you edit the source after you first start the blender session
import imp
imp.reload(hm)
imp.reload(cell_classes)


###################  PARAMETERS  #####################
ANGLE = 20

# CELL PARAMETERS
NUM_FIELD_CELLS = 25
FIELD_CELL_SCALE = (1,1,1)
FIELD_CELL_SIZE = 0.04

NUM_CURVE_CELLS = 15
CURVE_CELL_SCALE = (2.3,1,1)
CURVE_CELL_SIZE = 0.07
DEFORMATION_STRENGTH = 0.007

# CURVE PARAMS
CURVE_INTERVAL = (0,1)
CURVE_NOISE = Vector([0,0,0.04]) # Max absolute displacement of random curve cells along the curve
def CENTRAL_CURVE(t): # Time-parametrized curve in x-y-plane, t lying in CURVE_INTERVAL
    t *= 0.5*pi
    X = 0.5*cos(t)
    Y = 0.9*sin(t)
    Z = 0
    return Vector([X,Y,Z])

def SCALE_WEIGHT(t): # Time parametrized scale factor of cells along the curve
    t *= 0.5*pi
    return 0.7*(1-t) + 1.2*t

# GENERAL PARAMS
MIN_COORDS = Vector([0,0,0])
MAX_COORDS = Vector([1,1,0.04])
SEED = 123

###################  MAIN  METHOD  #####################
current_folder = os.path.dirname(os.path.realpath(__file__))
FOLDER = Path(current_folder).joinpath( "Images")
FOLDER = str(FOLDER)

try:
    if not os.path.exists(FOLDER):
        os.makedirs(FOLDER)
except OSError as error:
    print("Directory '%s' can not be created")


hm.delete_objects()

# set up scene
scene = bpy.context.scene
scene.render.engine = 'CYCLES'
scene.view_layers["ViewLayer"].use_pass_object_index =True
scene.use_nodes = True
scene.render.filepath =  FOLDER 
scene.render.image_settings.file_format = 'PNG'

# set up tree
tree = scene.node_tree
links = tree.links
render_nodes = tree.nodes['Render Layers']
output_file = tree.nodes.new("CompositorNodeOutputFile")
output_file.base_path = FOLDER 
output_file.format.file_format = 'PNG'

# set up camera
bpy.ops.object.camera_add() 
new_camera = bpy.context.object
bpy.context.scene.camera = new_camera

# Create single deformed cell in origin
random_cell = cell_classes.RandomDeformedCell(size=0.05, scale=(2,1,1), deformation_strength=0.02)

for i in range(NUM_FIELD_CELLS):
    cell_instance = random_cell.cell_object.copy()
    bpy.context.collection.objects.link(cell_instance)

    cell_instance.location = Vector([
        random.uniform(MIN_COORDS.x, MAX_COORDS.x),
        random.uniform(MIN_COORDS.y, MAX_COORDS.y),
        random.uniform(MIN_COORDS.z, MAX_COORDS.z)
    ])

    cell_instance.rotation_euler = (random.uniform(0, 2 * pi), random.uniform(0, 2 * pi), random.uniform(0, 2 * pi))
    cell_instance.scale = FIELD_CELL_SCALE 
    #pass_index = bpy.props.IntProperty(name="Pass Index", subtype='UNSIGNED')
    bpy.context.object.pass_index = i+1 

    maskid = tree.nodes.new('CompositorNodeIDMask')
    maskid.index = i+1
    
    scene.node_tree.links.new(render_nodes.outputs['IndexOB'], maskid.inputs[0])
    
    output_file.file_slots.new(f"mask_output_{i+1}")
    scene.node_tree.links.new(maskid.outputs['Alpha'], output_file.inputs[i+1])
    
bpy.ops.view3d.camera_to_view_selected()
bpy.context.scene.render.border_min_x = 0
bpy.context.scene.render.border_max_x = 1
bpy.context.scene.render.border_min_y = 0
bpy.context.scene.render.border_max_y = 1


bpy.ops.render.render(write_still=True)
