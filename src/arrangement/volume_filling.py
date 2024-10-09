import bpy
import numpy as np
import random

from itertools import combinations
from math import pi, sin, cos
from mathutils import Vector, Matrix

from src.utils.geometry import random_unit_vector


def fill_volume(ratios, density, attributes, volume, seed=None):
    '''
    Idea: Fill volume first with few large cell types that are randomly placed within the bounding box of the volume.
    Then fill the remaining space with small cell types densely.
    '''
    assert len(ratios) == len(attributes), "list of ratios and attributes need to have the same lengths."
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    # TODO: Deal with ratios / counts
    # sum = np.sum(ratios)
    # normalized_ratios = [ratio/sum for ratio in ratios]
    # MAX_COUNT = 600
    # counts = [int(ratio*MAX_COUNT) for ratio in normalized_ratios]
    counts = (0,20,100,20,10) # total: 150


    # Fill sparse and large cell types first
    attribute_list = [attribute for attribute, count in zip(attributes, counts) for _ in range(count)]
    non_fill_attribute_list = [attribute for attribute in attribute_list if not attribute.cell_type.name == "LYM"]
    random.shuffle(attribute_list)
    lattice_points = get_sorted_lattice_points(volume)
    # NOTE: We need more random samples than the number of seeds to place since randomly sampled seeds might intersect. - ck
    # We either stop when all seeds are place or we reach max_samples and no further seeds can be placed.
    random_seeds = []
    max_samples = 3*len(non_fill_attribute_list)
    for sample_idx in range(max_samples):
        pos = random.choice(lattice_points) # sample from bounding box including epithelial volumes so the density is same across images. - ck
        attribute = non_fill_attribute_list[len(random_seeds)]
        if not do_seeds_intersect(random_seeds, (pos, attribute)):
            random_seeds.append((pos, attribute))
        sample_idx += 1
        if len(random_seeds) == len(non_fill_attribute_list):
            break
    # Remove seeds within epithelial volume
    placed_first_seeds = [seed for seed in random_seeds if is_inside(seed[0], volume, seed[1].size)]

    # Fill LYM
    # NOTE: We fill the full bounding box with LYM seeds (smallest cell types) in a tree-like fashion starting from the center of the bbox.
    # Only in the end do we prune seeds which intersect non-fill cell types or lie within the epi volume.
    # NOTE: This is uglyyy & hacky. Generalize this in the future, in case you want to fill remaining volume with different cell types. - ck
    fill_attribute = [attribute for attribute in attributes if attribute.cell_type.name == "LYM"][0]
    MAX_COUNT = 2000 # Tree placing algorithm ends after ahaving placed MAX_COUNT of seeds.
    M = 12 # Number of latitude directions on the sphere
    N = 8 # Number of longitude directions on the sphere
    DENSITY_DISTANCE = 0.0 # Minimum distance between any two seeds
    # Get a set of directions that cover the sphere
    directions = [Vector((cos(2*pi*i/M) * sin(pi*j/N), sin(2*pi*i/M) * sin(pi*j/N), cos(pi*j/N))) for i in range(M) for j in range(1,N)] + [Vector((0,0,1)), Vector((0,0,-1))]
    

    min_corner, max_corner = get_min_max_corners_of_bbox(volume)
    center = (min_corner + max_corner) / 2
    root = (center, fill_attribute)
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
                new_attribute = fill_attribute
                new_position = seed[0] + direction*(seed[1].size + new_attribute.size + DENSITY_DISTANCE)
                if not do_seeds_intersect(placed_seeds, (new_position, new_attribute)) and is_inside_bbox(new_position, min_corner, max_corner):
                    placed_seeds.append((new_position, new_attribute))
                    placed_seeds_by_level[level].append((new_position, new_attribute))
                if len(placed_seeds) >= MAX_COUNT:
                    break
            if len(placed_seeds) >= MAX_COUNT:
                break
        if len(placed_seeds) >= MAX_COUNT:
            break
    print(f"{len(placed_seeds_by_level)} levels used.")

    # Remove seeds within epithelial volume and seeds that intersect sparse & large cell types created in first part. 
    placed_seeds = [seed for seed in placed_seeds if is_inside(seed[0], volume, seed[1].size) and not do_seeds_intersect(placed_first_seeds, seed, min_distance=0)]
    placed_seeds += placed_first_seeds
    random.shuffle(placed_seeds)
    assert density <= 1.0 and density >= 0.0, "Invalid density"
    placed_seeds = placed_seeds[:int(len(placed_seeds)*density)]

    points_per_attribute = []
    for attribute in attributes: 
        points = [seed[0] for seed in placed_seeds if seed[1] == attribute]
        points_per_attribute.append((points, attribute))
    return points_per_attribute


def do_seeds_intersect(seeds, candidate_seed, min_distance=0):
        current_point, current_attribute = candidate_seed
        for seed_point, seed_attribute in seeds:
            if (seed_point - current_point).length <= current_attribute.size + seed_attribute.size + min_distance:
                return True
        return False


def is_inside_bbox(seed_point, min_corner, max_corner):
    return (min_corner[0] <= seed_point[0] <= max_corner[0] and min_corner[1] <= seed_point[1] <= max_corner[1] and min_corner[2] <= seed_point[2] <= max_corner[2])
    

def get_sorted_lattice_points(volume):
    MAX_COUNT = 6000 # Default: 6000 / Maximum number of cells placed in the volume
    DELTA = 0.01 # Default: 0.01 / Lattice spacing, smaller values yield more random placements but increases processing time
    
    min_corner, max_corner = get_min_max_corners_of_bbox(volume)

    lattice_counts = (int((max_corner[0] - min_corner[0]) / DELTA), int((max_corner[1] - min_corner[1]) / DELTA), int((max_corner[2] - min_corner[2]) / DELTA))
    lattice_points = [min_corner + Vector((x*DELTA, y*DELTA, z*DELTA)) for x in range(lattice_counts[0]) for y in range(lattice_counts[1]) for z in range(lattice_counts[2])]
    midpoint = min_corner + Vector((lattice_counts[0]/2, lattice_counts[1]/2, lattice_counts[2]/2)) * DELTA
    random.shuffle(lattice_points)
    lattice_points = lattice_points[:MAX_COUNT]
    return sorted(lattice_points, key=lambda point: (point - midpoint).length)


def get_min_max_corners_of_bbox(volume):
    bounds = volume.bound_box
    world_corners = [volume.matrix_world @ Vector(corner) for corner in bounds]
    min_corner = Vector((min(x for x, y, z in world_corners), min(y for x, y, z in world_corners), min(z for x, y, z in world_corners)))
    max_corner = Vector((max(x for x, y, z in world_corners), max(y for x, y, z in world_corners), max(z for x, y, z in world_corners)))
    return min_corner, max_corner


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


def nearest_vertex(obj, q):
    if len(obj.data.vertices) == 0:
        return None, None
    vertices = [v for v in obj.data.vertices]
    distances = [(obj.matrix_world @ v.co - q).length for v in vertices]
    min_index = np.argmin(distances)
    return vertices[min_index], distances[min_index]