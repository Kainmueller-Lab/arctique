import bpy
import numpy as np
import random

from itertools import combinations
from mathutils import Vector


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
        
def add_point_cloud(locations, radius, name=None):
    # Create a small sphere object for each base point
    for idx, location in enumerate(locations):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=location)
        sphere = bpy.context.active_object
        sphere.name = f"Point_{name}_{idx}" if name!=None else f"Point_{idx}"
        
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

def generate_points_per_type(counts, radii, types, mesh):
    assert len(counts) == len(radii), "Counts and radii must have the same length"
    points_per_type = []
    # Sort counts and radii by radii
    zipped_data = list(zip(radii, counts, types))
    sorted_data = sorted(zipped_data, key=lambda x: x[0])
    
    for radius, count, type in sorted_data:
        points = generate_points(count, 2*radius, mesh)
        points_per_type.append((points, radius, type))
    assert len(counts) == len(points_per_type), "List points_per_type is not the same length as list counts"
    return points_per_type



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
    bounding_box = mesh.bound_box
    world_matrix = mesh.matrix_world
    world_vertices = [world_matrix @ Vector(vertex[:]) for vertex in bounding_box]
    
    min_x = min(vertex[0] for vertex in world_vertices)
    max_x = max(vertex[0] for vertex in world_vertices)
    min_y = min(vertex[1] for vertex in world_vertices)
    max_y = max(vertex[1] for vertex in world_vertices)
    min_z = min(vertex[2] for vertex in world_vertices)
    max_z = max(vertex[2] for vertex in world_vertices)
    return ((min_x, max_x), (min_y, max_y), (min_z, max_z))

def get_box_volume(bounding_box):
    xs, ys, zs = bounding_box
    return (xs[1]-xs[0]) * (ys[1]-ys[0]) * (zs[1]-zs[0])  

def random_point_in_bbox(bounding_box):
    xs, ys, zs = bounding_box
    return Vector((random.uniform(xs[0], xs[1]), random.uniform(ys[0], ys[1]), random.uniform(zs[0], zs[1])))
    

#############  MAIN

delete_objects()

# TODO:
# Given a list of n object attributes like size and scale and a list of ratios
# And also given a closed mesh
# Goal: Populate the inside volume with ellipsoids of given atributes that appear according to their ratio
'''
A Monte Carlo code would not be too difficult to write, based on what you have already:

    Insert the desired number of smaller spheres into the volume, at positions chosen at random, checking for overlap as you do at present.
    Conduct Monte Carlo: loop over all N

spheres sequentially. Attempt to move each sphere by a small amount in a random direction. Test for overlap with the other N−1
spheres. If there is any overlap, reject the move and leave the sphere in its original place.
At intervals, attempt to increase the diameter of all the spheres by a very small amount, checking all 12N(N−1)
overlaps. If there are any overlaps, keep the diameter unchanged and continue.
'''

# TODO
# - Add different sizes
# - Add different scales
# - Add differen orientations
# - Add blowup algorithm (Monte Carlo)


NUM_POINTS = 200
TYPES = ["A", "B"]
RADII = [0.1, 0.05]
COUNT_RATIOS = [0.2, 0.8]

# Get counts
sum = sum(COUNT_RATIOS)
normalized_ratios = [ratio/sum for ratio in COUNT_RATIOS]
counts = [int(ratio*NUM_POINTS) for ratio in normalized_ratios]
min_dist = 2 * max(RADII)

# Add a mesh object
bpy.ops.mesh.primitive_torus_add()
mesh = bpy.context.active_object

# Generate points inside mesh with given minimum distance
points_per_type = generate_points_per_type(counts, RADII, TYPES, mesh)

for points, radius, type in points_per_type:
    add_point_cloud(points, radius, type)
#mesh.hide_viewport = True
