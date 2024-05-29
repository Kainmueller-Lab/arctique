import bpy
import numpy as np
import random

from itertools import combinations
from mathutils import Vector

def get_shrinked_copy(obj, shrink_value):
    copied_obj = obj.copy()
    copied_obj.data = obj.data.copy()
    bpy.context.collection.objects.link(copied_obj)
    for vertex in copied_obj.data.vertices:
        vertex.co -= vertex.normal * shrink_value
    return copied_obj 

def generate_points(num_points, min_distance, mesh, use_strict_boundary=True):
    # Create shrinked copy of mesh
    shrink_value = min_distance / 2
    bounding_mesh = get_shrinked_copy(mesh, shrink_value)  if use_strict_boundary else mesh
    # Compute max number of iterations
    bbox = compute_bbox(bounding_mesh)
    max_radius = min_distance / 2
    box_volume = compute_bbox_volume(bbox)
    max_iterations = upper_limit_points(box_volume, max_radius)
    # Find num_points points inside the mesh that have at least min_distance
    points = []
    random_points = random_points_in_bbox(bbox, max_iterations)
    for new_point in random_points:	
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
        if len(points) == num_points:
            break
    if use_strict_boundary:
        bpy.data.objects.remove(bounding_mesh, do_unlink=True)
    return points

def fill_volume(counts, attributes, volume, use_strict_boundary, seed=None):
    radii = [attribute.size for attribute in attributes]
    types = [attribute.cell_type for attribute in attributes]
    assert len(counts) == len(radii), "Counts and radii must have the same length"

    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    points_per_type = []
    # Sort counts and radii by radii, starting with maximal radius
    zipped_data = list(zip(radii, counts, types))
    sorted_data = sorted(zipped_data, key=lambda x: x[0], reverse=True)
    for radius, count, type in sorted_data:
        # Compute max number of iterations
        bbox = compute_bbox(volume)
        box_volume = compute_bbox_volume(bbox)
        max_iterations = upper_limit_points(box_volume, radius)

        shrink_value = 1.5 * radius if use_strict_boundary else 0 # NOTE: Reduce 1.5 in case nuclei distance to mesh boundary is too large.
        points = []
        random_points = random_points_in_bbox(bbox, max_iterations)
        for new_point in random_points:
            if len(points) == count:
                break
            if is_inside(new_point, volume, shrink_value):   
                if is_far_from_points_per_type(new_point, points_per_type, radius):
                    if is_far_from_points(new_point, points, 2*radius):
                        points.append(new_point)
        points_per_type.append((points, radius, type))
    assert len(counts) == len(points_per_type), "List points_per_type is not the same length as list counts"
    return points_per_type

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

def nearest_vertex(obj, q):
    vertices = [v for v in obj.data.vertices]
    v_coords = [obj.matrix_world @ v.co for v in obj.data.vertices]
    distances = [np.linalg.norm(p-q) for p in v_coords]
    min_index = np.argmin(distances)
    nearest_vertex = vertices[min_index]
    min_distance = distances[min_index]
    return nearest_vertex, min_distance

def is_inside(point, obj, shrink=0):
    '''
    Checks whether a given point lies inside a shrinked version of the given mesh.
    '''
    nearest_vert, _ = nearest_vertex(obj, point)
    normal = nearest_vert.normal
    closest = nearest_vert.co - shrink*normal
    direction = closest - point
    return direction.dot(normal) > 0.01

def random_points_in_bbox(bbox, count):
    xs, ys, zs = bbox
    random_points = []
    for _ in range(int(count)):
        random_point = Vector((random.uniform(xs[0], xs[1]), random.uniform(ys[0], ys[1]), random.uniform(zs[0], zs[1])))
        random_points.append(random_point)
    return random_points

def compute_bbox(obj):
    box_vertices = [obj.matrix_world @ Vector(vertex[:]) for vertex in obj.bound_box]
    min_x = min(vertex[0] for vertex in box_vertices)
    max_x = max(vertex[0] for vertex in box_vertices)
    min_y = min(vertex[1] for vertex in box_vertices)
    max_y = max(vertex[1] for vertex in box_vertices)
    min_z = min(vertex[2] for vertex in box_vertices)
    max_z = max(vertex[2] for vertex in box_vertices)
    return ((min_x, max_x), (min_y, max_y), (min_z, max_z))

def compute_bbox_volume(bbox):
    xs, ys, zs = bbox
    return (xs[1] - xs[0]) * (ys[1] - ys[0]) * (zs[1] - zs[0])