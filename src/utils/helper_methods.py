import bpy
import numpy as np
from mathutils import Matrix, Vector


def move_selection(offset_vector):
    selection = bpy.context.selected_objects
    for obj in selection:
        obj.location += offset_vector
        
        
def delete_objects():
    # Select all objects 
    for obj in bpy.context.scene.objects:
        obj.hide_viewport = False
    bpy.ops.object.select_all(action='SELECT')

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

def random_unit_vector():
    # Generate random values for x, y, and z
    x, y, z = np.random.uniform(-1, 1, 3)

    # Create a vector
    vector = np.array([x, y, z])

    # Normalize the vector to make it a unit vector
    unit_vector = vector / np.linalg.norm(vector)

    return unit_vector
