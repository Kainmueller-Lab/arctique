import bpy
import numpy as np
import random

from itertools import combinations
from mathutils import Vector

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

def fill_volume(counts, attributes, mesh, padding):
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
