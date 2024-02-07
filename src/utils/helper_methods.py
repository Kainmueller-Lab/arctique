import bpy
import math
import numpy as np
import random
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
    # Remove orphaned meshes
    for mesh in bpy.data.meshes:
        bpy.data.meshes.remove(mesh, do_unlink=True)


def add_point_cloud(locations, radius):
    # Create a small sphere object for each base point
    for idx, location in enumerate(locations):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=location)
        sphere = bpy.context.active_object
        sphere.name = f"Point_{idx}"


def get_cube_points(min_coords, max_coords, padding = 0):
    points = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                point = [max_coords[0]+padding if i else min_coords[0]-padding,
                        max_coords[1]+padding if j else min_coords[1]-padding,
                        max_coords[2]+padding if k else min_coords[2]-padding]
                points.append(point)
    return points

def get_octogon_points(min_coords, max_coords, padding = 0):
    points = []
    cube_center = [0.5 * (i + j) for i, j in zip(min_coords, max_coords)]
    for i in range(3):
        point1, point2 = cube_center.copy(), cube_center.copy()
        point1[i] = max_coords[i] + padding
        points.append(point1)
        point2[i] = min_coords[i] - padding
        points.append(point2)
    return points

def get_lattice_points(min_coords, max_coords):
    points = []
    cube_center = Vector([0.5 * (i + j) for i, j in zip(min_coords, max_coords)])
    diameter = Vector([j - i for i, j in zip(min_coords, max_coords)])
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for k in [-1, 0, 1]:
                if not (i==0 and j==0 and k==0):
                    point = cube_center + Vector([a*b for (a,b) in zip(diameter, [i,j,k])])
                    points.append(point)
    return points

def add_box(min_coords, max_coords):
    bpy.ops.mesh.primitive_cube_add(size=1, enter_editmode=False, align='WORLD', location=((max_coords[0]+min_coords[0])/2, (max_coords[1]+min_coords[1])/2, (max_coords[2]+min_coords[2])/2))
    bpy.context.active_object.dimensions = [(max_coords[i] - min_coords[i]) for i in range(3)]
    bpy.context.active_object.location = [(max_coords[i] + min_coords[i]) / 2 for i in range(3)]
    return bpy.context.active_object

def intersect_with_object(target_objects, box_object):
    # Iterate through each target object
    for target_object in target_objects:
        boolean = target_object.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
        boolean.operation = 'INTERSECT'
        boolean.object = box_object

        # Apply the boolean modifier and remove the cube object
        bpy.ops.object.modifier_apply({"object": target_object}, modifier="Boolean Modifier")
    return target_objects

def subtract_object(target_objects, subtract_object):
    # Iterate through each target object
    for target_object in target_objects:
        boolean = target_object.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
        boolean.operation = 'DIFFERENCE'
        boolean.use_self = True
        boolean.object = subtract_object

        # Apply the boolean modifier and remove the cube object
        bpy.ops.object.modifier_apply({"object": target_object}, modifier="Boolean Modifier")
    return target_objects

def shrinkwrap(cell_objects, nuclei_scale=1):
    nucleus_objects = []
    for cell_object in cell_objects:
        bpy.ops.mesh.primitive_cube_add(enter_editmode=False, align='WORLD', location=cell_object.location)
        nucleus_object = bpy.context.active_object
        shrinkwrap = nucleus_object.modifiers.new(name="Shrinkwrap Modifier", type='SHRINKWRAP')
        shrinkwrap.target = cell_object
        bpy.ops.object.modifier_apply(modifier="Shrinkwrap Modifier")
        subsurf = nucleus_object.modifiers.new("Subsurface Modifier", type='SUBSURF')
        subsurf.levels = 2
        bpy.ops.object.modifier_apply(modifier="Subsurface Modifier")
        nucleus_object.scale = (nuclei_scale, nuclei_scale, nuclei_scale)
        nucleus_objects.append(nucleus_object)
    return nucleus_objects
        
def remove_objects(object_list):
    for obj in object_list:
        bpy.data.objects.remove(obj, do_unlink=True)

def get_objects_with(string):
    # Create a list to store matching objects
    object_list = []
    # Iterate through all objects
    for obj in bpy.data.objects:
        # Check if the object name starts with the defined pattern
        if obj.name.startswith(string):
            object_list.append(obj)
    return object_list



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

def rotate_objects(objects, alpha):
    # Convert angle to radians
    alpha_rad = math.radians(alpha)
    # Select the objects and activate the context
    bpy.context.view_layer.objects.active = objects[0]
    bpy.ops.object.select_all(action='DESELECT')
    for obj in objects:
        obj.select_set(True)
    # Apply the rotation using bpy.ops.transform.rotate
    bpy.ops.transform.rotate(value=alpha_rad, orient_axis='Z', orient_type='GLOBAL',
                             orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                             orient_matrix_type='GLOBAL', constraint_axis=(False, False, True),
                             mirror=False, use_proportional_edit=False, proportional_edit_falloff='SMOOTH',
                             proportional_size=1, use_proportional_connected=False, use_proportional_projected=False,
                             snap=False, snap_elements={'INCREMENT'}, use_snap_project=False, snap_target='CLOSEST',
                             use_snap_self=True, use_snap_edit=True, use_snap_nonedit=True, use_snap_selectable=False)
    bpy.ops.object.select_all(action='DESELECT')
    return  

def translate_objects(objects, location):
    # Select the objects and activate the context
    bpy.context.view_layer.objects.active = objects[0]
    bpy.ops.object.select_all(action='DESELECT')
    for obj in objects:
        obj.select_set(True)
    # Apply translation using bpy.ops.transform.translate
    bpy.ops.transform.translate(value=location, orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                                orient_matrix_type='GLOBAL', constraint_axis=(False, False, True),
                                mirror=False, use_proportional_edit=False, proportional_edit_falloff='SMOOTH',
                                proportional_size=1, use_proportional_connected=False, use_proportional_projected=False,
                                snap=False, snap_elements={'INCREMENT'}, use_snap_project=False, snap_target='CLOSEST',
                                use_snap_self=True, use_snap_edit=True, use_snap_nonedit=True, use_snap_selectable=False)
    # Deselect all objects
    bpy.ops.object.select_all(action='DESELECT')
    return  

def is_point_inside_mesh(mesh, point):
    # Transform the point to the mesh's local coordinates
    local_point = mesh.matrix_world.inverted() @ point
    # Calculate the normalized direction vector from the point to the reference position
    ray_direction = Vector([0.0, 0.0, 1.0])
    # Set up the ray casting parameters
    do_intersect, _, _, _ = mesh.ray_cast(local_point, local_point + 100 * ray_direction)
    return do_intersect


def do_intersect(obj1, obj2):
    # Set the active object for the boolean operation
    bpy.context.view_layer.objects.active = obj1
    # Apply the Boolean modifier with INTERSECT operation
    bpy.ops.object.modifier_add(type='BOOLEAN')
    bpy.context.object.modifiers["Boolean"].operation = 'INTERSECT'
    bpy.context.object.modifiers["Boolean"].use_self = False
    bpy.context.object.modifiers["Boolean"].object = obj2
    # Execute the Boolean operation
    bpy.ops.object.modifier_apply({"object": obj1}, modifier="Boolean")
    # Check if the result has geometry (intersection occurred)
    if len(obj1.data.vertices) > 0: 
        print(f"Intersecting verts: {len(obj1.data.vertices)}")
    intersection_exists = bool(obj1.data.vertices)
    # Remove the boolean modifier
    bpy.ops.object.modifier_remove({"object": obj1}, modifier="Boolean")
    return intersection_exists

def get_objects_with(string):
    # Create a list to store matching objects
    object_list = []
    # Iterate through all objects
    for obj in bpy.data.objects:
        # Check if the object name starts with the defined pattern
        if obj.name.startswith(string):
            object_list.append(obj)
    return object_list

def generate_lattice_parameters(theta):
   '''
   Theta defines a lattice of 4 isosceles parallelograms of area 1.
   Each such parallelogram defines a rotated ellipsoid fitting into it.
   '''
   theta = np.radians(theta)
   length = 1 / np.sin(theta)
   mu = 0.6 # NOTE: This is handpicked and steers the distance between crypts in the lattice. - ck
   x_axis = mu * length * np.cos(theta/2)
   y_axis = mu * length * np.sin(theta/2)
   v = Vector([length,0,0])
   w = Vector([length*np.cos(theta), length*np.sin(theta), 0])
   max_delta_scale = 0.05
   max_delta_angle = 10
   ico_scales = [(x_axis + max_delta_scale*random.uniform(-1,1), y_axis + max_delta_scale*random.uniform(-1,1)) for _ in range(4)]
   angles = [np.degrees(theta)/2 + max_delta_angle*random.uniform(-1,1) for _ in range(4)]
   base_pt = Vector([-0.5, -0.5, 0.5])
   lattice = [Vector([0,0,0]), v, w, v+w]
   centers = [base_pt + vec for vec in lattice]
   return ico_scales, angles, centers
