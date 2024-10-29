import random

from collections import namedtuple
from mathutils import Vector

from src.arrangement.deformation import deform_objects
from src.arrangement.surface_filling import fill_surface
from src.arrangement.volume_filling import fill_volume
from src.arrangement.voronoi import *
from src.objects.cells import CellAttribute
from src.utils.geometry import *
from src.utils.helper_methods import *
import src.utils.helper_methods as hm


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

        
class VolumeFill(CellArrangement):
    '''
    A volume fill arrangement. Given a volume mesh, fills the volume with randomly placed nuclei of different cell types without intersection.
    A maximal number of nuclei to place must be given, if number is set too high, Nnuclei are filled until not additional nucleus can be placed into the volume.
    Ratios of corresponding types must be given and determine how many nuclei of each type should be placed.
    If strict_boundary is set to true, nuclei objects will be placed only inside the mesh, otherwise only their locations will be inside the mesh.
    '''
    def __init__(self, mesh, density, types, ratios, seed=None, rescaling=0):
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
        self.types = types
        self.attributes = [CellAttribute.from_type(type) for type in self.types]
        self.ratios = ratios
        # Generate points inside mesh with given minimum distance
        self.points_per_attribute = fill_volume(
            self.ratios, self.density, self.attributes, self.mesh, self.seed)


    def add(self):
        for locations, attribute in self.points_per_attribute:
            for idx, location in enumerate(locations):
                direction = Vector(random_unit_vector(seed=idx+self.seed))
                cell_objects = attribute.add_cell_objects(location, direction, apply_subdivide=True)
                cytoplasm = None
                if len(cell_objects) == 2:
                    cytoplasm = cell_objects[1]
                    cytoplasm.name = f"Cytoplasm_Type_{attribute.cell_type.name}_{idx}"
                    self.objects.append(cytoplasm)
                    self.cytoplasm.append(cytoplasm)
                nucleus = cell_objects[0]
                hm.shade_switch(nucleus, flat=True)
                nucleus.name = f"Nucleus_Type_{attribute.cell_type.name}_{idx}"
                self.objects.append(nucleus)
                self.nuclei.append(nucleus)
        # NOTE: We currently deform only objects which are not EOS nuclei and FIBs, because those nuclei are based on metaball objects. - ck
        def is_deformable(obj):
            return not (obj.name.startswith('Nucleus_Type_EOS') or obj.name.startswith('Nucleus_Type_FIB'))
        deformable_objects = [obj for obj in self.objects if is_deformable(obj)]
        # NOTE: One could move the general deformation to the end such that also EPIs can be deformed.
        # However I think that those deforms will look bad due to difference in nucleus size.
        # In the future apply a more detailed noise deformation field for different regions. - ck
        deform_objects(deformable_objects)
        convert2mesh_list(self.objects)

    
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
    def __init__(self, mesh_obj, surface_obj, type, bounding_box=None, rescaling=0):
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
        self.radius = self.attribute.size*(1+rescaling)
        self.max_count = 1000 # NOTE: The volume will be filled up with maximally this number of nuclei
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
        placed_nuclei = []
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
            #subdivide(nucleus, self.attribute.subdivision_levels, apply=False, viewport=False)

            if self.type.name == "GOB":
                nucleus.name = f"Goblet_Type_{self.type.name}_{idx}"
            else:
                nucleus.name = f"Nucleus_Type_{self.type.name}_{idx}"
            placed_nuclei.append(nucleus)
            self.objects.append(nucleus)
            self.nuclei.append(nucleus)
        convert2mesh_list(self.objects)
        remove_objects(region_objects + [surface_obj])

        # add subdivision at once
        hm.subdivide_list(placed_nuclei, self.attribute.subdivision_levels)
        convert2mesh_list(placed_nuclei)

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