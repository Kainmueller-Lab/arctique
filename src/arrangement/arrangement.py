import random

from collections import namedtuple
from mathutils import Vector

from src.objects.cells import CellAttribute
from src.utils.geometry import *
from src.utils.helper_methods import *
from src.utils.surface_filling import fill_surface
from src.utils.volume_filling import fill_volume, fill_volume_old
from src.utils.voronoi import *
import src.utils.helper_methods as hm


NucleusSeed = namedtuple('NucleusSeed', 'centroid scale direction')
# NOTE: A NucleusSeed object is a 3d-ellipsoid defined by its 3d centroid, its scale
# (3d vector containing the length of axes in x, y and z direction), and its direction
# (3d normalized vector to which the x-axis will point).


############################################################
# Permutation table for Perlin noise
permutation = [151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225,
               140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247,
               120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177,
               33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165,
               71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211,
               133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25,
               63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196,
               135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217,
               226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206,
               59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119, 248,
               152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, 129, 22,
               39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246,
               97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51,
               145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84,
               204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114,
               67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180]

# Duplicate the permutation table to avoid overflow
permutation += permutation

def fade(t):
    # Fade function to smooth the interpolation
    return t * t * t * (t * (t * 6 - 15) + 10)

def lerp(t, a, b):
    # Linear interpolation
    return a + t * (b - a)

def grad(hash, x, y, z):
    # Gradient function to calculate the dot product
    h = hash & 15
    u = x if h < 8 else y
    v = y if h < 4 else (z if h == 12 or h == 14 else x)
    return ((u if (h & 1) == 0 else -u) +
            (v if (h & 2) == 0 else -v))

def noise(x, y, z):
    # Calculate the unit cube that the point falls in
    X = int(np.floor(x)) & 255
    Y = int(np.floor(y)) & 255
    Z = int(np.floor(z)) & 255

    # Find the relative x, y, z of the point in the cube
    x -= np.floor(x)
    y -= np.floor(y)
    z -= np.floor(z)

    # Compute fade curves for each of x, y, z
    u = fade(x)
    v = fade(y)
    w = fade(z)

    # Hash coordinates of the 8 cube corners
    A = permutation[X] + Y
    AA = permutation[A] + Z
    AB = permutation[A + 1] + Z
    B = permutation[X + 1] + Y
    BA = permutation[B] + Z
    BB = permutation[B + 1] + Z

    # Add blended results from 8 corners of the cube
    return lerp(w, lerp(v, lerp(u, grad(permutation[AA], x, y, z),
                                   grad(permutation[BA], x-1, y, z)),
                           lerp(u, grad(permutation[AB], x, y-1, z),
                                   grad(permutation[BB], x-1, y-1, z))),
                   lerp(v, lerp(u, grad(permutation[AA+1], x, y, z-1),
                                   grad(permutation[BA+1], x-1, y, z-1)),
                           lerp(u, grad(permutation[AB+1], x, y-1, z-1),
                                   grad(permutation[BB+1], x-1, y-1, z-1))))

import math
import random

def perlin(x, y, z, amplitude, frequency, octave_count, persistence, lacunarity):
    value = 0.0
    for _ in range(octave_count):
        value += amplitude * noise(x * frequency, y * frequency, z * frequency)
        amplitude *= persistence
        frequency *= lacunarity
    return value


### Deformation Function
def deformation(seed, input_vector, params):
    np.random.seed(seed)
    # Scale input vector to get noise coordinates
    noise_coords = input_vector
    
    # Generate Perlin noise values
    amplitude, frequency, octave_count, persistence, lacunarity = params
    amplitude = 1
    dx = perlin(noise_coords[0], noise_coords[1], noise_coords[2], amplitude, frequency, octave_count, persistence, lacunarity)
    dy = perlin(noise_coords[0] + 100, noise_coords[1] + 100, noise_coords[2] + 100, amplitude, frequency, octave_count, persistence, lacunarity)
    dz = perlin(noise_coords[0] + 200, noise_coords[1] + 200, noise_coords[2] + 200, amplitude, frequency, octave_count, persistence, lacunarity)
    
    # Create deformation vector
    deformation_vector = Vector((dx, dy, dz))
    return deformation_vector.normalized()
#############################################################################


class CellArrangement:
    # Class variable to keep track of the count
    count = 0

    def __init__(self):
        self.objects = [] # Contains all cell objects (i.e. nucleus and cytoplasm)
        self.nuclei = [] # Contains all nucleus objects
        self.cytoplasm = [] # Contains all cytoplasm objects
        self.name = None
        self.id = CellArrangement.count
        CellArrangement.count += 1

        
class VolumeFill(CellArrangement):
    '''
    A volume fill arrangement. Given a volume mesh, fills the volume with randomly placed nuclei of different cell types without intersection.
    A maximal number of nuclei to place must be given, if number is set too high, Nnuclei are filled until not additional nucleus can be placed into the volume.
    Ratios of corresponding types must be given and determine how many nuclei of each type should be placed.
    If strict_boundary is set to true, nuclei objects will be placed only inside the mesh, otherwise only their locations will be inside the mesh.
    '''
    def __init__(self, mesh, density, types, ratios, strict_boundary = True, seed=None, bounding_box=None):
        """
        Initializes a CellArrangement object with the given parameters.
        Fills a volume with randomly placed nuclei of different types without intersection.

        Parameters:
            - mesh: mesh of volume to populate with nuclei
            - density: density of nuclei to place
            - types: list of nuclei type types that should appear
            - ratios: list of ratios of nuclei types to populate
            - strict_boundary: if true will place nuclei fully inside the mesh, if false only centroids will be placed fully inside the mesh
        """
        super().__init__()
        self.seed = seed
        self.name = "VolumeFill"
        self.mesh = mesh
        self.density = density
        self.subdivision_levels = 2
        self.max_count = 1000 # Max number of nuclei to be placed
        self.types = types
        self.attributes = [CellAttribute.from_type(type) for type in self.types]
        self.ratios = ratios
        self.strict_boundary = strict_boundary # If true will place nuclei fully inside the mesh, if false only centroids will be placed fully inside the mesh
        # Get count
        sum = np.sum(ratios)
        normalized_ratios = [ratio/sum for ratio in ratios]
        self.counts = [int(ratio*self.max_count) for ratio in normalized_ratios]
        # Generate points inside mesh with given minimum distance
        # TODO: Generate based on Mahalanobis distance for scaled spheres
        self.points_per_attribute = fill_volume(
            self.counts, self.density, self.attributes, self.mesh, self.strict_boundary, self.seed)
        # remove all points outside of bounding box before adding objects
        # print('points per type', self.points_per_type)
        # if bounding_box is not None:
        #     self._remove_outside_bbox(bounding_box)

    def add(self):
        for locations, attribute in self.points_per_attribute:
            for idx, location in enumerate(locations):
                direction = Vector(random_unit_vector())
                cell_objects = attribute.add_cell_objects(location, direction, apply_subdivide=True)
                cytoplasm = None
                if len(cell_objects) == 2:
                    cytoplasm = cell_objects[1]
                    cytoplasm.name = f"Cytoplasm_Type_{attribute.cell_type}_{idx}"
                    self.objects.append(cytoplasm)
                    self.cytoplasm.append(cytoplasm)
                nucleus = cell_objects[0]
                hm.shade_switch(nucleus, flat=True)
                nucleus.name = f"Nucleus_Type_{attribute.cell_type}_{idx}"
                self.objects.append(nucleus)
                self.nuclei.append(nucleus)

        # TODO: Add deformation class and put code there
        import time
        start = time.time()
        amplitude = 0.25
        wavelength = 0.05
        octave_count = 2
        persistence = 1.0
        lacunarity = 2.0
        parameters = (amplitude, 1 / wavelength, octave_count, persistence, lacunarity)
        for obj in self.objects:
            centroid = Vector(np.mean([obj.matrix_world @ v.co for v in obj.data.vertices], axis=0))
            for v in obj.data.vertices:
                rad = (obj.matrix_world @ v.co - centroid).length
                v.co = v.co + deformation(42, v.co, parameters) * rad * amplitude
            subdivide(obj, 1)
        print(f"Deformed {len(self.objects)} objects in {time.time() - start} seconds")
        #convert2mesh_list(self.objects)
    
    def _remove_outside_bbox(self, bbox):
        deleted_objects = []
        for obj in self.points_per_attribute:
            # construct bounding box
            pass
            
        #     if not bounding_boxes_intersect(obj, bbox):
        #         deleted_objects.append(obj)
        # remove_objects(deleted_objects)
        # for obj in deleted_objects:
        #     self.objects.remove(obj)


class VoronoiFill(CellArrangement):
    def __init__(self, mesh_obj, surface_obj, type, bounding_box=None):
        """
        Initializes a CellArrangement object with the given parameters.
        Takes as input a mesh_obj which needs to consist of an outer and an inner wall mesh, e.g. an annulus.
        Then the volume is cut into quasi-hexagonal compartments whose size depend on the given attributes.
        Insides each compartment will then be placed a icosphere that fits into this compartment.

        Parameters:
            - mesh_obj: object of volume to populate with nuclei
            - attribute: nuclei type attribute that should populate the mesh
            - bounding_box: optional bounding box to delete compartments before placing nuclei
        """
        super().__init__()
        self.name = "VoronoiFill"
        self.mesh_obj = mesh_obj
        self.surface_obj = surface_obj
        self.type = type
        self.bounding_box = bounding_box
        self.attribute = CellAttribute.from_type(self.type)
        self.radius = self.attribute.size
        self.max_count = 500 # NOTE: The volume will be filled up with maximally this number of nuclei
        self.surface_subdivision_levels = 2 # NOTE: Increase this for finer subdvision and more quasi-random placement, but will be slower
        self.add_nuclei() 


    def add_nuclei(self):
        surface_obj = remove_vertical_and_horizontal_faces(self.surface_obj)
        subdivide_object(surface_obj, self.surface_subdivision_levels)

        # Sort for sampling
        outer_data = [(surface_obj.matrix_world @ v.co, v.normal) for v in surface_obj.data.vertices]
        if len(outer_data) == 0:
            return []
        else:
            if len(outer_data[0]) == 0:
                return []
        root = outer_data[0][0]
        sorted_data = self.sort_data_by_root_dist(outer_data, root)
        sorted_vs = [v for v, _ in sorted_data]
        choice, _ = self.sample_centers(sorted_vs, 2*self.radius, self.max_count)
        
        # Create Voronoi regions
        min_coords = [min([v[0] for v in choice]), min([v[1] for v in choice]), min([v[2] for v in choice])]
        max_coords = [max([v[0] for v in choice]), max([v[1] for v in choice]), max([v[2] for v in choice])]
        padding = 1.0
        auxiliary_points = get_lattice_points(min_coords, max_coords, padding)
        vor = Voronoi(choice + auxiliary_points)
        fr_points = finite_region_points(vor)
        ridges = compute_faces_by_seeds(vor, fr_points)

        region_objects = []
        for point_idx in fr_points:
            region = add_region(vertices = vor.vertices,
                               faces = ridges[point_idx],
                               idx = point_idx)
            region_objects.append(region)
        region_objects = intersect_with_object(region_objects, self.mesh_obj)

        print('part 1 done')
        # delete compartments outside of bounding box -> we could put this into a dedicated function
        if self.bounding_box is not None:
            deleted_objects = []
            for obj in region_objects:
                if not bounding_boxes_intersect(obj, self.bounding_box):
                    deleted_objects.append(obj)
            remove_objects(deleted_objects)
            for obj in deleted_objects:
                region_objects.remove(obj)
            print('deleted compartments', len(deleted_objects))

        # Place nuclei
        for idx, obj in enumerate(region_objects):
            remove_loose_vertices(obj) # NOTE: For some strange reason the intersection with "FAST" solver leads to ca 1000 loose vertices per region. - ck
            prism_coords = [obj.matrix_world @ v.co for v in obj.data.vertices]
            # Don't create nuclei objetcs if Voronoi region too small or too large (avoid artifacts)
            if len(prism_coords) < 4 or diameter(prism_coords) > 10*self.radius or diameter(prism_coords) < self.radius:
                continue

            bpy.ops.mesh.primitive_ico_sphere_add(radius=3, location=centroid(prism_coords))
            nucleus = bpy.context.active_object
            shrinkwrap(obj, nucleus, apply=False, viewport=True)
            smoothen_object(nucleus, self.attribute.smooth_factor, self.attribute.smooth_roundness, apply=False, viewport=True)
            subdivide(nucleus, self.attribute.subdivision_levels, apply=False, viewport=True)

            if self.type.name == "GOB":
                nucleus.name = f"Goblet_Type_{self.type.name}_{idx}"
            else:
                nucleus.name = f"Nucleus_Type_{self.type.name}_{idx}"
            self.objects.append(nucleus)
            self.nuclei.append(nucleus)
        convert2mesh_list(self.objects)
        remove_objects(region_objects + [surface_obj])

        print('part 2 done')

        return []

    def split_vertices(self, mesh):
        inner, outer = [], []
        for v in mesh.vertices:
            prod = np.dot(v.co, v.normal)
            if prod > 0:
                outer.append(v)
            else:
                inner.append(v)
        return inner, outer
    
    def sort_data_by_root_dist(self, data, v):
        def distance_to_v(v1, v2):
            return np.linalg.norm(v1 - v2)
        return sorted(data, key=lambda tuple: distance_to_v(tuple[0], v))
    
    def sample_centers(self, verts, dist, max_count):
        sampled_verts = []
        sampled_ids = []
        for idx, v in enumerate(verts):
            can_be_placed = True
            for q in sampled_verts:
                if (v-q).length < dist:
                    can_be_placed = False
                    break
            if can_be_placed:
                sampled_verts.append(v)
                sampled_ids.append(idx)
            if len(sampled_verts) == max_count:
                break
        return sampled_verts, sampled_ids

    def add(self):
        # NOTE: Remove and refactor later
        pass