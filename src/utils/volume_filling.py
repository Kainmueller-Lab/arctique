import bpy
import numpy as np
import random

from itertools import combinations
from math import pi, sin, cos
from mathutils import Vector, Matrix

from src.utils.geometry import random_unit_vector

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

def generate_grid_points(bbox, delta):
    """
    Generate a list of all grid points inside a bounding box with a given grid distance.
    
    Parameters:
    bbox (tuple): Bounding box defined as ((min_x, max_x), (min_y, max_y), (min_z, max_z))
    delta (float): Grid distance

    Returns:
    list: List of tuples representing the grid points
    """
    (min_x, max_x), (min_y, max_y), (min_z, max_z) = bbox

    x_points = [x for x in np.arange(min_x, max_x + delta, delta)]
    y_points = [y for y in np.arange(min_y, max_y + delta, delta)]
    z_points = [z for z in np.arange(min_z, max_z + delta, delta)]

    grid_points = [Vector((x, y, z)) for x in x_points for y in y_points for z in z_points]
    
    return grid_points

def bbox_center(bbox):
    (min_x, max_x), (min_y, max_y), (min_z, max_z) = bbox
    return (min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2

def fill_volume(counts, density, attributes, volume, use_strict_boundary, seed=None):
    assert len(counts) == len(attributes), "counts and attributes not same length :(("
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    MAX_COUNT = 900
    M = 16
    N = 12
    # Get a set of directions that cover the sphere
    directions = [Vector((cos(2*pi*i/M) * sin(pi*j/N), sin(2*pi*i/M) * sin(pi*j/N), cos(pi*j/N))) for i in range(M) for j in range(1,N)] + [Vector((0,0,1)), Vector((0,0,-1))]
    density_distance = 0.0
    print(density_distance)
    density = 1
    use_strict_boundary = False
    # TODO:
    # - Make the params input params
    # - Implement strict boundary

    attribute_list = [attribute for attribute, count in zip(attributes, counts) for _ in range(count)]
    random.shuffle(attribute_list)
    
    bounds = volume.bound_box
    world_corners = [volume.matrix_world @ Vector(corner) for corner in bounds]
    min_corner = Vector((min(x for x, y, z in world_corners), min(y for x, y, z in world_corners), min(z for x, y, z in world_corners)))
    max_corner = Vector((max(x for x, y, z in world_corners), max(y for x, y, z in world_corners), max(z for x, y, z in world_corners)))

    def is_inside_bbox(seed_point, min_corner, max_corner):
        return (min_corner[0] <= seed_point[0] <= max_corner[0] and min_corner[1] <= seed_point[1] <= max_corner[1] and min_corner[2] <= seed_point[2] <= max_corner[2])

    def do_intersect(seeds, candidate_seed, min_distance=0):
        current_point, current_attribute = candidate_seed
        for seed_point, seed_attribute in seeds:
            if (seed_point - current_point).length <= current_attribute.size + seed_attribute.size + min_distance:
                return True
        return False
    
        # Find random root inside mesh
    while True:
        pos = Vector((random.uniform(min_corner[0], max_corner[0]), random.uniform(min_corner[1], max_corner[1]), random.uniform(min_corner[2], max_corner[2])))
        if is_inside(pos, volume, 0):
            root = (pos, attributes[0])
            break

    placed_seeds = [root]
    placed_seeds_by_level = {0 : [root]}
    
    while True:
        level = len(placed_seeds_by_level)
        placed_seeds_by_level[level] = []
        if len(placed_seeds_by_level[level-1]) == 0:
            break
        for seed in placed_seeds_by_level[level-1]:
            # rotate all directions by random 3d angle
            rotated_directions = [Matrix.Rotation(random.uniform(0, 2 * pi), 4, Vector((random.random(), random.random(), random.random())).normalized()) @ v for v in directions] 
            for direction in rotated_directions:
                new_attribute = attribute_list[len(placed_seeds)]
                new_position = seed[0] + direction*(seed[1].size + new_attribute.size + density_distance)
                if not do_intersect(placed_seeds, (new_position, new_attribute)) and is_inside(new_position, volume, new_attribute.size) and is_inside_bbox(new_position, min_corner, max_corner):
                    placed_seeds.append((new_position, new_attribute))
                    placed_seeds_by_level[level].append((new_position, new_attribute))
                if len(placed_seeds) >= MAX_COUNT:
                    break
            if len(placed_seeds) >= MAX_COUNT:
                break
        if len(placed_seeds) >= MAX_COUNT:
            break
    print(f"{len(placed_seeds_by_level)} levels used.")
    random.shuffle(placed_seeds)
    placed_seeds = placed_seeds[:int(len(placed_seeds)*density)]

    points_per_attribute = []
    for attribute in attributes: 
        points = [seed[0] for seed in placed_seeds if seed[1] == attribute]
        points_per_attribute.append((points, attribute))
    return points_per_attribute


def fill_volume_old(counts, density, attributes, volume, use_strict_boundary, seed=None):
    assert len(counts) == len(attributes), "counts and attributes not same length :(("
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    MAX_COUNT = 20000 # Default: 5000 / Maximum number of cells placed in the volume
    DELTA = 0.005 # Default: 0.01 / Lattice spacing, smaller values yield more random placements but increases processing time
    density_distance = 0.0
    density = 1
    use_strict_boundary = False
    # TODO:
    # - Make the params input params
    # - Implement strict boundary
    # - Turn seeds into points_per_type

    attribute_list = [attribute for attribute, count in zip(attributes, counts) for _ in range(count)]
    random.shuffle(attribute_list)
    
    bounds = volume.bound_box
    world_corners = [volume.matrix_world @ Vector(corner) for corner in bounds]
    # def add_point_cloud(locations, radius, name=None):
    #     # Create a small sphere object for each base point
    #     for idx, location in enumerate(locations):
    #         bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=location)
    #         sphere = bpy.context.active_object
    #         sphere.name = f"Point_{name}_{idx}" if name!=None else f"Point_{idx}"
    # add_point_cloud(world_corners, 0.1)
    min_corner = Vector((min(x for x, y, z in world_corners), min(y for x, y, z in world_corners), min(z for x, y, z in world_corners)))
    max_corner = Vector((max(x for x, y, z in world_corners), max(y for x, y, z in world_corners), max(z for x, y, z in world_corners)))

    lattice_counts = (int((max_corner[0] - min_corner[0]) / DELTA), int((max_corner[1] - min_corner[1]) / DELTA), int((max_corner[2] - min_corner[2]) / DELTA))
    lattice_points = [min_corner + Vector((x*DELTA, y*DELTA, z*DELTA)) for x in range(lattice_counts[0]) for y in range(lattice_counts[1]) for z in range(lattice_counts[2])]
    midpoint = min_corner + Vector((lattice_counts[0]/2, lattice_counts[1]/2, lattice_counts[2]/2)) * DELTA
    random.shuffle(lattice_points)
    lattice_points = lattice_points[:MAX_COUNT]
    lattice_points = sorted(lattice_points, key=lambda point: (point - midpoint).length)
    #assert (len(lattice_points) == len(attribute_list)), "attribute_list and lattice_points must have same length :("

    def do_intersect(seeds, candidate_seed, min_distance=0):
        current_point, current_attribute = candidate_seed
        for seed_point, seed_attribute in seeds:
            if (seed_point - current_point).length <= current_attribute.size + seed_attribute.size + min_distance:
                return True
        return False

    seeds = []
    for current_point in lattice_points:
        candidate_seed = (current_point, attribute_list[len(seeds)])
        if do_intersect(seeds, candidate_seed, min_distance=density_distance):
            continue
        seeds.append(candidate_seed)
    print("Total Seeds:")
    print(len(seeds))
    if use_strict_boundary:
        seeds = [seed for seed in seeds if is_inside(seed[0], volume, seed[1].size)]
    else:
        seeds = [seed for seed in seeds if is_inside(seed[0], volume)]
    random.shuffle(seeds)
    seeds = seeds[:int(len(seeds)*density)]
    print("Filtered Seeds:")
    print(len(seeds))

    points_per_attribute = []
    for attribute in attributes: 
        points = [seed[0] for seed in seeds if seed[1] == attribute]
        points_per_attribute.append((points, attribute))
    return points_per_attribute



def fill_volume_very_old(counts, density, attributes, volume, use_strict_boundary, seed=None):
    radii = [attribute.size for attribute in attributes]
    types = [attribute.cell_type for attribute in attributes]
    assert len(counts) == len(radii), "Counts and radii must have the same length"

    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    max_count = np.sum(counts)
    types_to_place = [type_ for type_, count in zip(types, counts) for _ in range(count)]
    assert len(types_to_place) == max_count, "Cat error ^.^"
    random.shuffle(types_to_place)
    points_per_type = {type_: [] for type_ in types}

    bbox = compute_bbox(volume)
    max_radius = max(radii) * (2-density) * 0.5
    grid_delta = 0.07 # Allows for more points to be placed near boundary but smaller value take longer

    volume_grid = generate_grid_points(bbox, grid_delta)
    volume_grid = [pt + Vector(random_unit_vector()) * 0.5 * grid_delta for pt in volume_grid]
    root = bbox_center(bbox)
    sorted_volume_grid = sorted(volume_grid, key=lambda point: np.linalg.norm(np.array(point) - np.array(root)))

    shrink_value = 1.2 * max_radius if use_strict_boundary else 0 # NOTE: Reduce coeff in case nuclei distance to mesh boundary is too large.

    placed = []
    for new_point in sorted_volume_grid:
        if len(placed) == max_count:
            break
        # If sampled point is within 2*max_rad of any placed point, skip, otherwise keep point
        current_type = types_to_place[len(placed)] 
        if is_inside(new_point, volume, shrink_value):   
                if is_far_from_points(new_point, placed, 2*max_radius):
                    placed.append(new_point)
                    #print(f"Placing {current_type}, {len(placed)} / {max_count}")
                    points_per_type[current_type].append(new_point)

    result = [(points_per_type[type_], radii[types.index(type_)], type_) for type_ in types]
    return result

    # NOTE: This is the old placing algorithm based on precomuted counts instead of density
    # # Sort counts and radii by radii, starting with maximal radius
    # zipped_data = list(zip(radii, counts, types))
    # sorted_data = sorted(zipped_data, key=lambda x: x[0], reverse=True)
    # bbox = compute_bbox(volume)
    # box_volume = compute_bbox_volume(bbox)
    # for radius, count, type in sorted_data:
    #     # Compute max number of iterations
    #     max_iterations = upper_limit_points(box_volume, radius)

    #     shrink_value = 1.5 * radius if use_strict_boundary else 0 # NOTE: Reduce 1.5 in case nuclei distance to mesh boundary is too large.
    #     points = []
    #     random_points = random_points_in_bbox(bbox, max_iterations)
    #     for new_point in random_points:
    #         if len(points) == count:
    #             break
    #         if is_inside(new_point, volume, shrink_value):   
    #             if is_far_from_points_per_type(new_point, points_per_type, radius):
    #                 if is_far_from_points(new_point, points, 2*radius):
    #                     points.append(new_point)
    #     points_per_type.append((points, radius, type))
    # assert len(counts) == len(points_per_type), "List points_per_type is not the same length as list counts"
    # return points_per_type

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

def nearest_vertex_old(obj, q):
    vertices = [v for v in obj.data.vertices]
    v_coords = [obj.matrix_world @ v.co for v in obj.data.vertices]
    distances = [np.linalg.norm(p-q) for p in v_coords]
    # TODO fix error if necessary
    if distances == []:
        return None, None
    min_index = np.argmin(distances)
    nearest_vertex = vertices[min_index]
    min_distance = distances[min_index]
    return nearest_vertex, min_distance

def nearest_vertex(obj, q):
    if len(obj.data.vertices) == 0:
        return None, None
    vertices = [v for v in obj.data.vertices]
    distances = [(obj.matrix_world @ v.co - q).length for v in vertices]
    min_index = np.argmin(distances)
    return vertices[min_index], distances[min_index]

def is_inside(point, obj, shrink=0):
    '''
    Checks whether a given point lies inside a shrinked version of the given mesh.
    '''
    nearest_vert, _ = nearest_vertex(obj, point)
    if nearest_vert is None:
        return False
    normal = nearest_vert.normal
    closest = obj.matrix_world @ nearest_vert.co - shrink*normal
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