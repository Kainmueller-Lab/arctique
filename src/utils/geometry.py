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

def sample_points_on_mesh(obj, nuclei_limit, cut_box, max_translate=None):
    # Generate seeds for Voronoi diagram
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.mode_set(mode='EDIT')
    mesh = bmesh.from_edit_mesh(obj.data)
    bmesh.ops.triangulate(mesh, faces=mesh.faces[:]) # Triangulate # NOTE: Check if this is truly necessary. - ck
    #print(f"Mesh has {len(mesh.faces)} faces and {len(mesh.verts)} vertices.")
    points = []
    for face in mesh.faces:
        centroid = face.calc_center_median()
        points.append(centroid)
        for v in face.verts:
            v_loc = obj.matrix_world @ v.co
            p = 0.5*(centroid + v_loc)
            points.append(p)
    #        q = 0.25*(centroid - v_loc) + centroid # NOTE: Can be used in case the nuclei are not dense enough. - ck
    #        points.append(q)
    for vert in mesh.verts:
        v_loc = obj.matrix_world @ vert.co
        points.append(v_loc)
    bpy.ops.object.mode_set(mode='OBJECT')
    # Retain only points that lie between min and max CUT_BOX
    min_cut_box, max_cut_box = cut_box
    points = [p for p in points if all(min_c <= val <= max_c for val, min_c, max_c in zip(p, min_cut_box, max_cut_box))]
    if len(points) > nuclei_limit:
        print(f"WARNING: Adding only {nuclei_limit} nuclei instead {len(points)} nuclei due to long compute time.\nYou can reduce the slice thickness to generate less nuclei.")
        points = random.sample(points, k=nuclei_limit)
    if max_translate is not None:
        points = [[p[i] + random.uniform(-1, 1)*max_translate for i in range(3)] for p in points]
    return points


def add_point_cloud(locations, radius):
    # Create a small sphere object for each base point
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

        
def get_shrinked_copy(mesh, shrink_value):
    copied_mesh = mesh.data.copy()
    copied_object = bpy.data.objects.new("Copied_Object", copied_mesh)
    bpy.context.collection.objects.link(copied_object)
    # Apply the shrink/fatten operation directly to the copied mesh data
    for vertex in copied_object.data.vertices:
        vertex.co += vertex.normal * shrink_value
    return copied_object
        
        
def generate_points(num_points, min_distance, mesh, padding=True):
    # Create shrinked copy of mesh
    shrink_value = -min_distance / 2
    bounding_mesh = get_shrinked_copy(mesh, shrink_value)  if padding else mesh
    # Compute max number of iterations
    bounding_box = get_bounding_box(bounding_mesh)
    max_radius = min_distance / 2
    box_volume = get_box_volume(bounding_box)
    max_iterations = upper_limit_points(box_volume, max_radius)
    print(f"Max point count: {max_iterations}")

    # Find num_points points inside the mesh that have at least min_distance
    points = []
    iterations = 0
    while len(points) < num_points and iterations < max_iterations:
        new_point = random_point_in_bbox(bounding_box)	
        valid_point = True
        if is_inside(new_point, bounding_mesh):        
            for point in points:
                distance = np.linalg.norm(new_point - point)  # Calculate distance between new point and existing points
                if distance < min_distance:
                    valid_point = False
                    break
        else:
            valid_point = False
        if valid_point:
            points.append(new_point)
        iterations += 1
    if padding:
        bpy.data.objects.remove(bounding_mesh, do_unlink=True)
    return points

def generate_points_per_type(counts, attributes, mesh, padding):
    radii = [attribute.size for attribute in attributes]
    types = [attribute.cell_type for attribute in attributes]
    assert len(counts) == len(radii), "Counts and radii must have the same length"

    points_per_type = []
    # Sort counts and radii by radii, starting with maximal radius
    zipped_data = list(zip(radii, counts, types))
    sorted_data = sorted(zipped_data, key=lambda x: x[0], reverse=True)
    for radius, count, type in sorted_data:
        # Create shrinked copy of mesh
        shrink_value = -radius
        bounding_mesh = get_shrinked_copy(mesh, shrink_value)  if padding else mesh
        # Compute max number of iterations
        bounding_box = get_bounding_box(bounding_mesh) # NOTE: Apparently the bounding box is always centered in the origin
        box_volume = get_box_volume(bounding_box)
        max_iterations = upper_limit_points(box_volume, radius)

        # TODO: test this
        points = []
        iterations = 0
        while len(points) < count and iterations < max_iterations:
            new_point = random_point_in_bbox(bounding_box)
            if is_inside(new_point, bounding_mesh):   
                if is_far_from_points_per_type(new_point, points_per_type, radius):
                    if is_far_from_points(new_point, points, 2*radius):
                        points.append(new_point)
            iterations += 1
        if padding:
            bpy.data.objects.remove(bounding_mesh, do_unlink=True)
        points_per_type.append((points, radius, type))
    assert len(counts) == len(points_per_type), "List points_per_type is not the same length as list counts"
    # NOTE: All points are generated with an origin centered mesh, so we need to translate them back. - ck
    translated_points_per_type = [([pt + mesh.location for pt in points], radius, type) for points, radius, type in points_per_type]
    return translated_points_per_type

def is_far_from_points(pt, points, min_distance):
    for point in points:
        distance = np.linalg.norm(pt - point)  # Calculate distance between new point and existing points
        if distance < min_distance:
            return False
    return True

def is_far_from_points_per_type(pt, points_per_type, radius):
    for points, larger_radius, _ in points_per_type:
        if not is_far_from_points(pt, points, radius + larger_radius):
            return False
    return True

def minimum_distance(points):
    return min([np.linalg.norm(p1 - p2) for p1, p2 in combinations(points, 2)])

def upper_limit_points(volume, radius):
    # NOTE: Using here that the average density in an irregular packing of balls into a box is bounded by 64%.
    # See: https://en.wikipedia.org/wiki/Sphere_packing
    return np.floor((3*0.64*volume)/(4*np.pi*radius*radius*radius))

def is_inside(point, obj):
    _point = point - obj.location
    _, closest, nor, _ = obj.closest_point_on_mesh(_point)
    direction = closest - _point
    return direction.dot(nor) > 0

def get_bounding_box(mesh):
    box_vertices = [Vector(vertex[:]) for vertex in mesh.bound_box]
    min_x = min(vertex[0] for vertex in box_vertices)
    max_x = max(vertex[0] for vertex in box_vertices)
    min_y = min(vertex[1] for vertex in box_vertices)
    max_y = max(vertex[1] for vertex in box_vertices)
    min_z = min(vertex[2] for vertex in box_vertices)
    max_z = max(vertex[2] for vertex in box_vertices)
    return ((min_x, max_x), (min_y, max_y), (min_z, max_z))

def get_box_volume(bounding_box):
    xs, ys, zs = bounding_box
    return (xs[1]-xs[0]) * (ys[1]-ys[0]) * (zs[1]-zs[0])  

def random_point_in_bbox(bounding_box):
    xs, ys, zs = bounding_box
    return Vector((random.uniform(xs[0], xs[1]), random.uniform(ys[0], ys[1]), random.uniform(zs[0], zs[1])))

          
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
    VERTS_TO_MOVE = 3 # TODO: Take care of magical number
    PROPORTIONAL_SIZE = 0.4 # TODO: Necessary? If yes implement
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
    size = attribute.size
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


# NOTE: Can be deleted in the future. - ck
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
    boolean.operation = 'DIFFERENCE'
    boolean.object = outer_cylinder
    bpy.ops.object.modifier_apply({"object": mix_vol}, modifier="Boolean Modifier")
    remove_objects([inner_cylinder, outer_cylinder])
    return mix_vol, epi_vol