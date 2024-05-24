import random

from collections import namedtuple
from mathutils import Vector

from src.objects.cells import CellAttribute
from src.utils.geometry import *
from src.utils.helper_methods import *
from src.utils.surface_filling import fill_surface
from src.utils.volume_filling import fill_volume
from src.utils.voronoi import *


NucleusSeed = namedtuple('NucleusSeed', 'centroid scale direction')
# NOTE: A NucleusSeed object is a 3d-ellipsoid defined by its 3d centroid, its scale
# (3d vector containing the length of axes in x, y and z direction), and its direction
# (3d normalized vector to which the x-axis will point).


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

        
# TODO: Add blowup algorithm (Monte Carlo)
class VolumeFill(CellArrangement):
    '''
    A volume fill arrangement. Given a volume mesh, fills the volume with randomly placed nuclei of different cell types without intersection.
    A maximal number of nuclei to place must be given, if number is set too high, Nnuclei are filled until not additional nucleus can be placed into the volume.
    Ratios of corresponding types must be given and determine how many nuclei of each type should be placed.
    If strict_boundary is set to true, nuclei objects will be placed only inside the mesh, otherwise only their locations will be inside the mesh.
    '''
    def __init__(self, mesh, number, types, ratios, strict_boundary = True):
        """
        Initializes a CellArrangement object with the given parameters.
        Fills a volume with randomly placed nuclei of different types without intersection.

        Parameters:
            - mesh: mesh of volume to populate with nuclei
            - number: number of total nuclei to populate
            - types: list of nuclei type types that should appear
            - ratios: list of ratios of nuclei types to populate
            - strict_boundary: if true will place nuclei fully inside the mesh, if false only centroids will be placed fully inside the mesh
        """
        super().__init__()
        self.name = "VolumeFill"
        self.mesh = mesh
        self.subdivision_levels = 2
        self.number = number
        self.types = types
        self.attributes = [CellAttribute.from_type(type) for type in self.types]
        self.ratios = ratios
        self.strict_boundary = strict_boundary # If true will place nuclei fully inside the mesh, if false only centroids will be placed fully inside the mesh
        # Get count
        sum = np.sum(ratios)
        normalized_ratios = [ratio/sum for ratio in ratios]
        self.counts = [int(ratio*number) for ratio in normalized_ratios]
        # Generate points inside mesh with given minimum distance
        # TODO: Generate based on Mahalanobis distance for scaled spheres
        self.points_per_type = fill_volume(self.counts, self.attributes, self.mesh, self.strict_boundary)

    def add(self):
        for locations, radius, type in self.points_per_type:
            attribute = CellAttribute.from_type(type)
            for idx, location in enumerate(locations):
                direction = Vector(random_unit_vector())
                cell_objects = attribute.add_cell_objects(location, direction)
                cytoplasm = None
                if len(cell_objects) == 2:
                    cytoplasm = cell_objects[1]
                    cytoplasm.name = f"Cytoplasm_Type_{type.name}_{idx}"
                    self.objects.append(cytoplasm)
                    self.cytoplasm.append(cytoplasm)
                nucleus = cell_objects[0]
                nucleus.name = f"Nucleus_Type_{type.name}_{idx}"
                self.objects.append(nucleus)
                self.nuclei.append(nucleus)


class VoronoiFill(CellArrangement):
    def __init__(self, mesh_obj, count, type):
        """
        Initializes a CellArrangement object with the given parameters.
        Takes as input a mesh_obj which needs to consist of an outer and an inner wall mesh, e.g. an annulus.
        Then the volume is cut into quasi-hexagonal compartments whose size depend on the given attributes.
        Insides each compartment will then be placed a icosphere that fits into this compartment.

        Parameters:
            - mesh_obj: object of volume to populate with nuclei
            - count: number of total nuclei to populate
            - attribute: nuclei type attribute that should populate the mesh
        """
        super().__init__()
        self.name = "VoronoiFill"
        self.mesh_obj = mesh_obj
        self.count = count
        self.type = type
        self.attribute = CellAttribute.from_type(self.type)
        self.size_coeff = 0.7 # This can be adjusted. It scales the icospheres w.r.t. to the surrounding compartment. If it is set to 1 it maximally fits into the compartment but can lead to overlaps of meshes. - ck
        self.radius = self.attribute.size * self.attribute.scale[1] # We need the secong largest radius as tesselation distance between seed points
        self.subdivision_level = 4 # Reduce level for computation speed, increase for finer tessellation

        # Compute seeds for icospherical nuclei
        self.nuclei_seeds = self.compute_seeds()

    def compute_seeds(self):
        surface_obj = remove_top_and_bottom_faces(self.mesh_obj)
        subdivide_object(surface_obj, self.subdivision_level)
        mesh = surface_obj.data
        _, outer_vs = self.split_vertices(mesh)

        # Sort for sampling
        outer_data = [(surface_obj.matrix_world @ v.co, v.normal) for v in outer_vs]
        root = outer_data[0][0]
        sorted_data = self.sort_data_by_root_dist(outer_data, root)
        sorted_vs = [v for v, _ in sorted_data]
        choice, _ = self.sample_centers(sorted_vs, 2*self.radius, self.count)
        #choice_normals = [sorted_data[idx][1] for idx in choice_ids]
        
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

        # Place icospheres into Voronoi regions
        seeds = []
        for idx, obj in enumerate(region_objects):
            prism_coords = [obj.matrix_world @ v.co for v in obj.data.vertices]
            diam = diameter(prism_coords)
            c = centroid(prism_coords)
            height = np.sqrt(diam*diam - 4*self.radius*self.radius)
            scale = (0.5*self.size_coeff*height, self.size_coeff*self.radius, self.size_coeff*self.radius)
            direction = (c - choice[idx]) / np.linalg.norm((c - choice[idx]))
            if height < 5 * self.radius: # NOTE: Need this artifact check for now since sometimes a very long icosphere is created. - ck
                seeds.append(NucleusSeed(centroid=c, scale=scale, direction=direction))
        remove_objects(region_objects + [surface_obj])
        return seeds

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
        attribute = CellAttribute.from_type(self.type)
        for idx, seed in enumerate(self.nuclei_seeds):
            cell_objects = attribute.add_cell_objects(seed.centroid, seed.direction)
            cytoplasm = None
            if len(cell_objects) == 2:
                cytoplasm = cell_objects[1]
                cytoplasm.name = f"Cytoplasm_Type_{type.name}_{idx}"
                self.objects.append(cytoplasm)
                self.cytoplasm.append(cytoplasm)
            nucleus = cell_objects[0]
            nucleus.scale = tuple(s / attribute.size for s in seed.scale) # NOTE: Need to rescale since for EPI the scale depends on the Voronoi placement and cannot be given at construction. - ck
            nucleus.name = f"Nucleus_Type_{self.type.name}_{idx}"
            self.objects.append(nucleus)
            self.nuclei.append(nucleus)