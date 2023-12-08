import bpy
import numpy as np
import random

from math import radians, sin, cos, pi
from mathutils import Matrix, Vector

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
#random.seed(SEED)

################# HELPER METHODS ##############

def move_selection(offset_vector):
    selection = bpy.context.selected_objects
    for obj in selection:
        obj.location += offset_vector
        
        
def delete_pos_z_objects():
    # Select all objects with z location greater than 0
    for obj in bpy.context.scene.objects:
        if obj.location.z >= 0:
            obj.select_set(True)

    # Delete the selected objects
    bpy.ops.object.delete()
    

def lerp(a, b, t):
    return a*(1-t) + b*t
    

# Function to set the local x-axis orientation of an object along a direction vector
def set_orientation(obj, direction_vector):
    # Calculate the rotation matrix to align the object  with the direction vector
    rotation_matrix = Matrix.Translation(obj.location) @ direction_vector.to_track_quat('X').to_matrix().to_4x4()
    # Apply the rotation to the object
    obj.matrix_world = rotation_matrix

def compute_normal(curve_function, t): # Only in x-y-plane
    epsilon = 1e-4  # Small value to avoid division by zero
    # Compute the tangent vector
    tangent_vector = (curve_function(t + epsilon) - curve_function(t - epsilon))
    # Normalize the tangent vector
    tangent_vector /= np.linalg.norm(tangent_vector)
    assert tangent_vector.z == 0

    # Compute the normal vector (perpendicular to the tangent)
    normal_vector = Vector([tangent_vector.y, -tangent_vector.x, 0])
    return normal_vector

    
################ CLASSES #############
        

class RandomDeformedCell:
    def __init__(self, size=1.0, scale=(1,1,1), deformation_strength=0.1, orientation=None):
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.orientation = orientation

        # Create a cube
        bpy.ops.mesh.primitive_cube_add(size=self.size, location=(0, 0, 0))
        self.cell_object = bpy.context.active_object
        
        # Apply two levels of subdivision surface
        modifier = self.cell_object.modifiers.new("Subsurface Modifier", "SUBSURF")
        modifier.levels = 2
        
        # Orient
        if orientation:
            set_orientation(self.cell_object, orientation)

        # Scale 
        self.cell_object.scale = self.scale
        
        # Deform the vertices randomly
        self.deform_vertices()
        

    def deform_vertices(self):
        # Iterate through each vertex and deform its position
        for vertex in self.cell_object.data.vertices:
            original_position = vertex.co.copy()
            deformation_vector = Vector([
                random.uniform(-self.deformation_strength, self.deformation_strength),
                random.uniform(-self.deformation_strength, self.deformation_strength),
                random.uniform(-self.deformation_strength, self.deformation_strength)
            ])*Vector(self.scale)
            vertex.co = original_position + deformation_vector
            
            
class RandomCellField:
    def __init__(self, num_cells, min_coords, max_coords):
        self.num_cells = num_cells
        self.min_coords = min_coords
        self.max_coords = max_coords
        
    def generate_cells(self, cell_size, cell_scale, deformation_strength):
        for _ in range(self.num_cells):
            cell = RandomDeformedCell(cell_size, cell_scale, deformation_strength)
            cell.cell_object.location = Vector([
                random.uniform(self.min_coords.x, self.max_coords.x),
                random.uniform(self.min_coords.y, self.max_coords.y),
                random.uniform(self.min_coords.z, self.max_coords.z)
            ])
            
            
class CellCurve:
    def __init__(self, num_cells, curve_function, curve_noise, curve_interval=(0,1), scale_weight=None):
        self.num_cells = num_cells
        self.curve_function = curve_function
        self.curve_noise = curve_noise
        self.curve_interval = curve_interval
        self.scale_weight = scale_weight
        self.increment = (curve_interval[1]-curve_interval[0]) / num_cells
        
    def generate_cells(self, cell_size, cell_scale, deformation_strength):
        t = self.curve_interval[0]
        for idx in range(self.num_cells):
            # Change orientation of cells such that x-axis points to curve normal
            normal = compute_normal(self.curve_function, t)
            # Scale cells along the curve according to scale_weight
            current_scale = tuple([self.scale_weight(t) * s for s in cell_scale]) if self.scale_weight!=None else cell_scale
            cell = RandomDeformedCell(cell_size, current_scale, deformation_strength, normal)
            cell.cell_object.location = Vector([
                self.curve_function(t).x + random.uniform(-self.curve_noise.x, self.curve_noise.x),
                self.curve_function(t).y + random.uniform(-self.curve_noise.y, self.curve_noise.y),
                # TODO: Change this once we do not delete cells with zloc > 0
                self.curve_function(t).z + random.uniform(0, self.curve_noise.z) 
            ])
            # TODO: Handle increment weight better such that it can be tweaked by parameters
            increment_weight = self.scale_weight(t)*self.scale_weight(t) if self.scale_weight else 1
            current_increment = self.increment * increment_weight
            t += current_increment

        
        
###################  MAIN  METHOD  #####################

delete_pos_z_objects()

# Create single deformed cell in origin
random_cell = RandomDeformedCell(size=0.05, scale=(2,1,1), deformation_strength=0.02)

# Create uniformly distributed cell field
rc_field = RandomCellField(NUM_FIELD_CELLS, MIN_COORDS, MAX_COORDS)
#rc_field.generate_cells(FIELD_CELL_SIZE, FIELD_CELL_SCALE, DEFORMATION_STRENGTH)

c_curve = CellCurve(NUM_CURVE_CELLS, CENTRAL_CURVE, CURVE_NOISE, CURVE_INTERVAL, SCALE_WEIGHT)
c_curve.generate_cells(CURVE_CELL_SIZE, CURVE_CELL_SCALE, DEFORMATION_STRENGTH)




################### OLD CODE ########################

#bpy.ops.mesh.primitive_cube_add()
#so = bpy.context.active_object

## move object
#so.location[0]=2

## Rotate
#so.rotation_euler[0] += radians(ANGLE)

## create modifier
#mod_subsurf = so.modifiers.new("MyModifier", "SUBSURF")
#mod_subsurf.levels = 1

#mesh = so.data

## Define the range for random values (adjust as needed)
#min_value = -0.5
#max_value = 0.5

## Loop through all vertices and add a random 3D vector
#for vert in mesh.vertices:
#    random_vector = Vector([random.uniform(min_value, max_value) for _ in range(3)])
#    z = vert.co[2]
#    diff = random_vector
#    vert.co += diff

## Update the mesh with the new vertex positions
#mesh.update()




        
#move_selection(-Vector([1,1,1]))
