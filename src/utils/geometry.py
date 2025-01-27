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


def diameter(coords):
    return max([np.linalg.norm(v-w) for v,w in combinations(coords,2)])

def centroid(coords):
    return sum(coords, Vector((0.0, 0.0, 0.0))) / len(coords)

# Function to set the local x-axis orientation of an object along a direction vector
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

def random_unit_vector(seed=None):
    if seed is not None:
        np.random.seed(seed)
    # Generate random values for x, y, and z
    x, y, z = np.random.uniform(-1, 1, 3)
    # Create a vector
    vector = np.array([x, y, z])
    # Normalize the vector to make it a unit vector
    unit_vector = vector / np.linalg.norm(vector)
    return unit_vector

def rotate_objects(objects, alpha):
    """
    Rotate a list of objects by a given angle around the Z-axis.

    Parameters:
    objects (list): a list of objects to rotate
    alpha (float): the angle of rotation in degrees

    Returns:
    None
    """
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
    """
    Translate a list of objects by a given location.

    Parameters:
    objects (list): a list of objects to translate
    location (tuple of 3 floats): the location to translate the objects to

    Returns:
    None
    """
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
    """
    Determines if two Blender objects intersect by applying a Boolean INTERSECT operation.
    
    Parameters:
    obj1 (bpy.types.Object): The first Blender object to check for intersection.
    obj2 (bpy.types.Object): The second Blender object to check for intersection.
    
    Returns:
    bool: True if the objects intersect, False otherwise.
    """

    bpy.context.view_layer.objects.active = obj1
    # Apply the Boolean modifier with INTERSECT operation
    mod = bpy.ops.object.modifier_add(type='BOOLEAN')
    mod.show_viewport = False
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
    '''
    Removes all objects in the given list that are inside any of the given meshes.
    
    Parameters:
    objects (list of bpy.types.Object): The list of objects to check for intersection with the given meshes.
    empty_regions (list of bpy.types.Object): The list of meshes to check against.
    
    Returns:
    list of bpy.types.Object: The objects that were not inside any of the given meshes.
    '''
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
    """
    Generates a box with the given minimum and maximum coordinates.

    Parameters:
        min_coords (list or tuple of 3 floats): The minimum x, y, z coordinates of the box.
        max_coords (list or tuple of 3 floats): The maximum x, y, z coordinates of the box.

    Returns:
        bpy.types.Object: The generated box object.
    """
    bpy.ops.mesh.primitive_cube_add(size=1, enter_editmode=False, align='WORLD', location=((max_coords[0]+min_coords[0])/2, (max_coords[1]+min_coords[1])/2, (max_coords[2]+min_coords[2])/2))
    bpy.context.active_object.dimensions = [(max_coords[i] - min_coords[i]) for i in range(3)]
    bpy.context.active_object.location = [(max_coords[i] + min_coords[i]) / 2 for i in range(3)]
    return bpy.context.active_object

def intersect_with_object(target_objects, box_object):
    """
    Intersects each target object with the given box object.

    Parameters:
        target_objects (list of bpy.types.Object): The objects that will be intersected with the box.
        box_object (bpy.types.Object): The box object that will be used for the intersection.

    Returns:
        list of bpy.types.Object: The modified target objects.
    """
    for target_object in target_objects:
        boolean = target_object.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
        boolean.show_viewport = False
        boolean.operation = 'INTERSECT'
        boolean.object = box_object
        boolean.solver = 'FAST'
        bpy.ops.object.modifier_apply({"object": target_object}, modifier="Boolean Modifier")
    return target_objects

def subtract_object(target_objects, subtract_object):
    """
    Subtracts the given object from each target object.

    Parameters:
        target_objects (list of bpy.types.Object): The objects that will have the subtract object subtracted from them.
        subtract_object (bpy.types.Object): The object that will be used for the subtraction.

    Returns:
        list of bpy.types.Object: The modified target objects.
    """
    for target_object in target_objects:
        boolean = target_object.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
        boolean.show_viewport = False
        boolean.operation = 'DIFFERENCE'
        boolean.use_self = True
        boolean.object = subtract_object

        # Apply the boolean modifier and remove the cube object
        bpy.ops.object.modifier_apply({"object": target_object}, modifier="Boolean Modifier")
    return target_objects

def smoothen_object(obj, factor, iterations, apply=True, viewport=True):
    """
    Smooths the given object by applying a Smooth modifier with the given factor and number of iterations.

    Parameters:
        obj (bpy.types.Object): The object to be smoothed.
        factor (float): The smoothness factor.
        iterations (int): The number of iterations the smooth modifier will run.
        apply (bool): Whether to apply the modifier immediately. Defaults to True.
        viewport (bool): Whether to show the modifier in the viewport. Defaults to True.

    Returns:
        None
    """
    smooth_mod = obj.modifiers.new(name="Smooth Modifier", type='SMOOTH')
    if not viewport:
        smooth_mod.show_viewport = False
    smooth_mod.factor = factor
    smooth_mod.iterations = iterations
    if apply:
        bpy.ops.object.modifier_apply(modifier="Smooth Modifier")

def shrinkwrap(target, wrapper, apply=True, viewport=True):
    """
    Applies a Shrinkwrap modifier to the wrapper object to make it conform to the target object.

    Parameters:
        target (bpy.types.Object): The object that the wrapper object will conform to.
        wrapper (bpy.types.Object): The object that will be modified to conform to the target.
        apply (bool): Whether to apply the modifier immediately. Defaults to True.
        viewport (bool): Whether to show the modifier in the viewport. Defaults to True.

    Returns:
        None
    """
    shrinkwrap_mod = wrapper.modifiers.new(name="Shrinkwrap", type='SHRINKWRAP')
    if not viewport:
        shrinkwrap_mod.show_viewport = False
    shrinkwrap_mod.target = target
    shrinkwrap_mod.wrap_method = 'NEAREST_SURFACEPOINT'  # You can choose other methods like 'PROJECT', 'NEAREST_VERTEX', etc.
    bpy.context.view_layer.objects.active = wrapper
    if apply:
        bpy.ops.object.modifier_apply(modifier=shrinkwrap_mod.name)

def subdivide(obj, levels, apply=True, viewport=True):
    """
    Subdivides the given object by adding a Subdivision modifier to it.

    Parameters:
        obj (bpy.types.Object): The object to be subdivided.
        levels (int): The number of times the object will be subdivided.
        apply (bool): Whether to immediately apply the modifier. Defaults to True.
        viewport (bool): Whether to show the modifier in the viewport. Defaults to True.

    Returns:
        None
    """

    obj.modifiers.new(name="Subdivision", type='SUBSURF')
    if not viewport:
        obj.modifiers["Subdivision"].show_viewport = False
    obj.modifiers["Subdivision"].levels = levels
    if apply:
        bpy.ops.object.modifier_apply({"object": obj}, modifier="Subdivision")
    

def move_selection(offset_vector):
    """
    Moves all selected objects by a given offset.

    Parameters:
        offset_vector (Vector): The vector by which to move each selected object.

    Returns:
        None
    """

    selection = bpy.context.selected_objects
    for obj in selection:
        obj.location += offset_vector

def pos_value(x):
    """
    Returns the positive value of the input or zero if the input is negative.

    Parameters:
        x (float or int): The input value to evaluate.

    Returns:
        float or int: The input value if it is positive, otherwise zero.
    """

    return x if x>0 else 0

def perturb_vertices(mesh, scaled_deform_strength):
    """
    Perturbs the vertices of a mesh by moving them in a random direction according to
    their normal vector. This function is used to create a more realistic, irregular
    shape for the tissue surface.

    Parameters:
        mesh (bpy.types.Mesh): The mesh object to be perturbed.
        scaled_deform_strength (float): A scaling factor to control the strength of the
            perturbation. The perturbation amount is proportional to this value.

    Returns:
        None
    """
    VERTS_TO_MOVE = 10
    for _ in range(VERTS_TO_MOVE):
        direction = Vector(random_unit_vector())
        for w in mesh.vertices:
            mu = scaled_deform_strength * pos_value(np.dot(direction, w.normal)-0.8)
            w.co -= direction * mu
    
def rescale_obj(obj, size, scale):
    """
    Rescales the vertices of a mesh object to fit within a specified size and scale.

    This function first normalizes the mesh based on its maximum radius and then
    rescales it according to the provided size and scale factors.

    Parameters:
        obj (bpy.types.Object): The mesh object to be rescaled.
        size (float): The target size for normalization.
        scale (tuple of float): A tuple of scale factors for each axis (x, y, z).

    Returns:
        None
    """

    max_radius = max(np.linalg.norm(v.co) for v in obj.data.vertices)
    for v in obj.data.vertices:
        v.co *= size/max_radius
        v.co = Vector([v/sc for v, sc in zip(v.co, scale)])
            
def deform_mesh(obj, attribute):
    """
    Deforms a mesh object according to the provided attribute parameters.

    The deformation is done by perturbing the vertices of the mesh object in a
    random direction proportional to their normal vector. The amount of
    perturbation is proportional to the deform strength and size of the mesh.

    Parameters:
        obj (bpy.types.Object): The mesh object to be deformed.
        attribute (CellAttribute): An object containing the size, scale, and
            deformation strength of the object.

    Returns:
        None
    """
    size = attribute.size
    scale = attribute.scale
    deform_strength = attribute.deformation_strength
    scaled_deform_strength = Vector(s*size*deform_strength for s in scale) 
    mesh = obj.data
    perturb_vertices(mesh, scaled_deform_strength)

def bend_mesh(obj, bend):
    """
    Bends a mesh object randomly around its loal z-axis by a given factor.

    The bending is done by adding a Simple Deform modifier to the object and
    applying it. The modifier is set to bend the mesh around its local z-axis
    by a random angle between -2*PI*bend and 2*PI*bend.

    Parameters:
        obj (bpy.types.Object): The mesh object to be bent.
        bend (float): The bending factor. A factor of 0 results in no bending.

    Returns:
        None
    """

    modifier = obj.modifiers.new("Simple Deform Modifier", "SIMPLE_DEFORM")
    modifier.show_viewport = False
    modifier.deform_method = 'BEND'
    bpy.ops.object.empty_add(type='ARROWS', align='WORLD', location=obj.location, scale=(1, 1, 1))
    empty = bpy.context.active_object
    # Set the origin and deform axis
    modifier.origin = empty
    modifier.deform_axis = 'Z'
    modifier.angle = 2*np.pi*bend*random.uniform(-1,1)
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)
    bpy.ops.object.modifier_apply(modifier="Simple Deform Modifier")
    bpy.data.objects.remove(empty, do_unlink=True)

def remove_top_and_bottom_faces(obj):
    """
    Removes the top and bottom faces of a mesh object.

    The function removes all faces of the object that have a normal with a z-component close to -1 or 1 (i.e., the faces that are parallel to the xy-plane).

    Parameters:
        obj (bpy.types.Object): The mesh object whose faces are to be removed.

    Returns:
        bpy.types.Object: The modified object with the top and bottom faces removed.
    """
    mesh = bpy.data.meshes.new(name="ModifiedMesh")
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
    bm.to_mesh(mesh)
    bm.free()
    # Create a new object and link it to the scene
    new_obj = bpy.data.objects.new("ModifiedObject", mesh)
    new_obj.location = obj.location
    new_obj.scale = obj.scale
    bpy.context.collection.objects.link(new_obj)
    return new_obj

def remove_vertical_and_horizontal_faces(obj):
    """
    Removes the vertical and horizontal faces of a mesh object.

    The function removes all faces of the object that have a normal with a x- or y-component close to -1 or 1 (i.e., the faces that are parallel to the xz- or yz-plane).

    Parameters:
        obj (bpy.types.Object): The mesh object whose faces are to be removed.

    Returns:
        bpy.types.Object: The modified object with the vertical and horizontal faces removed.
    """
    mesh = bpy.data.meshes.new(name="ModifiedMesh")
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
        if abs(normal.z) > 1-epsilon or abs(normal.x) > 1-epsilon or abs(normal.y) > 1-epsilon:
            face.select_set(True)
    # Delete selected faces
    bmesh.ops.delete(bm, geom=[f for f in bm.faces if f.select], context='FACES')
    # Update the mesh data
    bm.to_mesh(mesh)
    bm.free()
    # Create a new object and link it to the scene
    new_obj = bpy.data.objects.new("ModifiedObject", mesh)
    new_obj.location = obj.location
    new_obj.scale = obj.scale
    bpy.context.collection.objects.link(new_obj)
    return new_obj

# NOTE: Can be deleted in the future. - ck
def add_dummy_objects(tissue, padding, vol_scale, surf_scale):
    # Create temporarily padded tissue
    tissue.tissue.scale = tuple(1+padding for _ in range(3))

    bpy.ops.mesh.primitive_cylinder_add() # Example bounding torus mesh
    cylinder = bpy.context.active_object
    #vol_obj.location = Vector(tissue.tissue.location) + Vector((0, 0, 0.5))
    cylinder.scale = vol_scale
    # NOTE: Necessary to transform the vertices of the mesh according to scale
    # It should be used when the object is created, but maybe there's a better place in the methds for it. ck
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
    # Intersect with tissue
    bpy.ops.mesh.primitive_cube_add(location=tissue.location)
    box = bpy.context.active_object
    box.scale = (1.1, 1.1, tissue.thickness/tissue.size)
    bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)
    vol_obj = subtract_object([box], cylinder)[0]
    remove_top_and_bottom_faces(vol_obj)
    remove_objects([cylinder])
    vol_obj.name = "Volume"

    bpy.ops.mesh.primitive_cylinder_add()
    surf_obj = bpy.context.active_object
    #surf_obj.location = Vector(tissue.tissue.location) + Vector((0, 0, 0.5))
    surf_obj.scale = surf_scale
    # NOTE: Necessary to transform the vertices of the mesh according to scale
    # It should be used when the object is created, but maybe there's a better place in the methds for it. ck
    bpy.ops.object.transform_apply(location=True, rotation=True, scale=True) 
    # Intersect with tissue
    intersect_with_object([surf_obj], tissue.tissue)
    remove_top_and_bottom_faces(surf_obj)
    surf_obj.name = "Surface"

    tissue.tissue.scale = (1,1,1)
    return vol_obj, surf_obj


def add_dummy_volumes(tissue, padding):
    bpy.ops.mesh.primitive_cube_add(size=tissue.size, location=tissue.location) 
    mix_vol = bpy.context.active_object
    mix_vol.name = "Mix_Volume"
    mix_vol.scale = (1 + padding, 1 + padding, tissue.thickness/tissue.size + padding)
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)

    bpy.ops.mesh.primitive_cylinder_add()
    epi_vol = bpy.context.active_object
    epi_vol.name = "Epi_Volume"
    epi_vol.scale = (0.6, 0.9, 1)
    bpy.ops.mesh.primitive_cylinder_add()
    inner_cylinder = bpy.context.active_object
    inner_cylinder.scale = (0.5, 0.7, 1)
    boolean = epi_vol.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
    boolean.operation = 'DIFFERENCE'
    boolean.object = inner_cylinder
    bpy.ops.object.modifier_apply({"object": epi_vol}, modifier="Boolean Modifier")

    # MIX_VOL
    boolean = epi_vol.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
    boolean.operation = 'INTERSECT'
    boolean.object = mix_vol
    bpy.ops.object.modifier_apply({"object": epi_vol}, modifier="Boolean Modifier")

    # EPI_VOL
    bpy.ops.mesh.primitive_cylinder_add()
    outer_cylinder = bpy.context.active_object
    outer_cylinder.scale = (0.6, 0.9, 1)
    boolean = mix_vol.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
    boolean.show_viewport = False
    boolean.operation = 'DIFFERENCE'
    boolean.object = outer_cylinder
    bpy.ops.object.modifier_apply({"object": mix_vol}, modifier="Boolean Modifier")
    remove_objects([inner_cylinder, outer_cylinder])
    return mix_vol, epi_vol

def remove_loose_vertices(obj):
    '''
    Removes all vertices of an object in the list that are not connected to a face
    '''
    bm = bmesh.new()
    mesh = obj.data
    bm.from_mesh(mesh)
    verts = [v for v in bm.verts if not v.link_faces]
    # equiv of bmesh.ops.delete(bm, geom=verts, context='VERTS')
    for v in verts:
        bm.verts.remove(v)
    bm.to_mesh(mesh)
    bm.clear()

def bounding_boxes_intersect(obj1, obj2):
    """
    Determines if the bounding boxes of two Blender objects intersect.

    Parameters:
    obj1 (bpy.types.Object): The first Blender object.
    obj2 (bpy.types.Object): The second Blender object.

    Returns:
    bool: True if the bounding boxes intersect, False otherwise.
    """
    # Get the bounding box coordinates of obj1
    bbox1 = [obj1.matrix_world @ Vector(corner) for corner in obj1.bound_box]
    min_x1 = min(v.x for v in bbox1)
    max_x1 = max(v.x for v in bbox1)
    min_y1 = min(v.y for v in bbox1)
    max_y1 = max(v.y for v in bbox1)
    min_z1 = min(v.z for v in bbox1)
    max_z1 = max(v.z for v in bbox1)

    # Get the bounding box coordinates of obj2
    bbox2 = [obj2.matrix_world @ Vector(corner) for corner in obj2.bound_box]
    min_x2 = min(v.x for v in bbox2)
    max_x2 = max(v.x for v in bbox2)
    min_y2 = min(v.y for v in bbox2)
    max_y2 = max(v.y for v in bbox2)
    min_z2 = min(v.z for v in bbox2)
    max_z2 = max(v.z for v in bbox2)

    # Check for overlap in each dimension
    intersect_x = (min_x1 <= max_x2 and max_x1 >= min_x2)
    intersect_y = (min_y1 <= max_y2 and max_y1 >= min_y2)
    intersect_z = (min_z1 <= max_z2 and max_z1 >= min_z2)

    return intersect_x and intersect_y and intersect_z