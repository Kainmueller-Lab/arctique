import bpy
import numpy as np
import random
import sys
import os
from math import radians, sin, cos, pi
from mathutils import Matrix, Vector

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

# Set random seed
# random.seed(SEED)

        
###################  MAIN  METHOD  #####################

hm.delete_objects()

# Create single deformed cell in origin
random_cell = cell_classes.RandomDeformedCell(size=0.05, scale=(2,1,1), deformation_strength=0.02)

# Create uniformly distributed cell field
rc_field = cell_classes.RandomCellField(NUM_FIELD_CELLS, MIN_COORDS, MAX_COORDS)
rc_field.generate_cells(FIELD_CELL_SIZE, FIELD_CELL_SCALE, DEFORMATION_STRENGTH)

c_curve = cell_classes.CellCurve(NUM_CURVE_CELLS, CENTRAL_CURVE, CURVE_NOISE, CURVE_INTERVAL, SCALE_WEIGHT)
c_curve.generate_cells(CURVE_CELL_SIZE, CURVE_CELL_SCALE, DEFORMATION_STRENGTH)


# TODO render same config in multiple variants 
# 1) full shading
# 2) classical instance masks
# 3) depth masks (or 3D labels)
# TODO add cell shading properties
