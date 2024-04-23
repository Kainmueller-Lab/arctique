import bmesh
import bpy
import math
import numpy as np
import random

from itertools import combinations
from mathutils import Vector, Matrix
from src.utils.helper_methods import *

def lerp(a, b, t):
    return a*(1-t) + b*t
 
def set_orientation(obj, direction_vector):
    """
    Rotates an object such that its local x-axis is oriented along a direction vector.

    Parameters:
    - obj: the object to orient
    - direction_vector: the vector defining the orientation direction

    Return:
    - None
    """
    # Calculate the rotation matrix to align the object  with the direction vector
    rotation_matrix = Matrix.Translation(obj.location) @ direction_vector.to_track_quat('X').to_matrix().to_4x4()
    # Apply the rotation to the object
    obj.matrix_world = rotation_matrix

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
    """
    Check if a point is inside a mesh.

    Args:
        mesh: The mesh to check against.
        point: The point to check.

    Returns:
        bool: True if the point is inside the mesh, False otherwise.
    """
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

def compute_mean_scale(object):
    '''
    Given a mesh object this method returns the mean scale of the object.
    The mean scale is defined as the geometric mean of the x, y and z-diameter of the mesh.      
    '''
    mesh = object.data
    min_coords = (min(vert.co[0] for vert in mesh.vertices),
                    min(vert.co[1] for vert in mesh.vertices),
                    min(vert.co[2] for vert in mesh.vertices))
    max_coords = (max(vert.co[0] for vert in mesh.vertices),
                    max(vert.co[1] for vert in mesh.vertices),
                    max(vert.co[2] for vert in mesh.vertices))
    diameter = [max - min for max, min in zip(max_coords, min_coords)]
    return (diameter[0]*diameter[1]*diameter[2]) ** (1/3)

def remove_objects_inside_mesh(objects, empty_regions):
    nuclei_to_remove = []
    nuclei_to_keep = []
    for nucleus_object in objects:
        for region in empty_regions:
            if is_point_inside_mesh(region, nucleus_object.location):
                #print(f"Removing nucleus {nucleus_object.name}")
                nuclei_to_remove.append(nucleus_object)
            else:
                nuclei_to_keep.append(nucleus_object)
    remove_objects(nuclei_to_remove)
    return nuclei_to_keep

def add_point_cloud(locations, radius):
    """
    Generates a point cloud of small spheres at the given locations.

    Parameters:
        locations (List[Tuple[float, float, float]]): A list of 3D coordinates where the spheres will be placed.
        radius (float): The radius of the spheres.

    Returns:
        None
    """
    for idx, location in enumerate(locations):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=location)
        sphere = bpy.context.active_object
        sphere.name = f"Point_{idx}"


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

def move_selection(offset_vector):
    selection = bpy.context.selected_objects
    for obj in selection:
        obj.location += offset_vector
          
def pos_value(x):
    return x if x>0 else 0

def perturb_vertices_old(mesh, deform_strength):
    VERTS_TO_MOVE = 8 # TODO: Take care of magical number
    PROPORTIONAL_SIZE = 0.5 # TODO: Necessary? If yes implement
    for _ in range(VERTS_TO_MOVE):
    # TODO: Test different deformations
        #source_v = mesh.vertices[random.randint(0, len(mesh.vertices) - 1)]
        #direction = source_v.normal
        direction = Vector(random_unit_vector())
        for w in mesh.vertices:
            mu = deform_strength * pos_value(np.dot(direction, w.normal))
            w.co += direction * mu

def perturb_vertices(mesh, scaled_deform_strength):
    VERTS_TO_MOVE = 2 # TODO: Take care of magical number
    PROPORTIONAL_SIZE = 0.7 # TODO: Necessary? If yes implement
    # TODO: Test different deformations
    for _ in range(VERTS_TO_MOVE):
        #source_v = mesh.vertices[random.randint(0, len(mesh.vertices) - 1)]
        #direction = source_v.normal
        direction = Vector(random_unit_vector())
        for w in mesh.vertices:
            mu = scaled_deform_strength * pos_value(np.dot(direction, w.normal)-PROPORTIONAL_SIZE)
            w.co -= direction * mu
    
def rescale_obj(obj, size, scale):
    max_radius = max(np.linalg.norm(v.co) for v in obj.data.vertices)
    for v in obj.data.vertices:
        v.co *= size/max_radius
        v.co = Vector([v/sc for v, sc in zip(v.co, scale)])
            
def deform_mesh(obj, attribute):
    size = attribute.radius
    scale = attribute.scale
    deform_strength = attribute.deformation_strength
    scaled_deform_strength = Vector(s*size*deform_strength for s in scale) 
    mesh = obj.data
    perturb_vertices(mesh, scaled_deform_strength)
    # TODO: Old deform
    # # TODO: Rescale such that diameter of mesh lies in bounding ball
    # #rescale_obj(obj, size, scale)
    # abs_deform_strength = size*deform_strength  
    # mesh = obj.data
    # perturb_vertices(mesh, abs_deform_strength)
    # # TODO: Rescale such that diameter of mesh lies in bounding ball
    # rescale_obj(obj, size, attribute.scale)


def remove_top_and_bottom_faces(obj):
    bm = bmesh.new()
    bm.from_mesh(obj.data)

    epsilon = 1e-3
    # Deselect all faces
    for face in bm.faces:
        face.select_set(False)
    
    # Select faces with vertical normals
    for face in bm.faces:
        normal = face.normal
        # Check if the z-component of the normal is close to -1 or 1
        if abs(normal.z) > 1-epsilon:
            face.select_set(True)
    
    # Delete selected faces
    bmesh.ops.delete(bm, geom=[f for f in bm.faces if f.select], context='FACES')
    
    # Update the mesh data
    bm.to_mesh(obj.data)
    bm.free()