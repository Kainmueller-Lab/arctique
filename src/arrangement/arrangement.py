import random

from collections import namedtuple
from mathutils import Vector

from src.objects.cells import Cell
from src.utils.geometry import *
from src.utils.helper_methods import *
from src.utils.surface_filling import fill_surface
from src.utils.volume_filling import fill_volume
from src.utils.voronoi import *

class CellArrangement:
    # Class variable to keep track of the count
    count = 0

    def __init__(self):
        self.objects = [] # Contains the nuclei objects
        self.name = None
        self.id = CellArrangement.count
        CellArrangement.count += 1

class CellList(CellArrangement):
    '''
    A cell list arrangement. Given a list of locations and attributes, places cells with given attributes at those locations.
    '''
    def __init__(self, cell_attributes, locations):
        """
        Initializes a new instance of the CellArrangement class.

        Parameters:
            cell_attributes (dict): A dictionary representing the cell attributes.
            locations (list): A list of locations.
        """
        super().__init__()
        self.cell_attributes = cell_attributes
        self.locations = locations
        self.name = "CellList"

    def add(self):
        for location in self.locations:
            # TODO: This is old code, and needs an update if CellLists turn out to be relevant. - ck 
            cell = Cell(location, self.id, self.type, self.cell_attributes)
            self.objects.append(cell)
        for cell in self.objects:
            cell.add()
        
# TODO: Add blowup algorithm (Monte Carlo)
class VolumeFill(CellArrangement):
    '''
    A volume fill arrangement. Given a volume mesh, fills the volume with randomly placed nuclei of different attributes without intersection.
    A maximal number of nuclei to place must be given, if number is set too high, Nnuclei are filled until not additional nucleus can be placed into the volume.
    Ratios of corresponding attributes must be given and determine how many nuclei of each type should be placed.
    If strict_boundary is set to true, nuclei objects will be placed only inside the mesh, otherwise only their locations will be inside the mesh.
    '''
    def __init__(self, mesh, number, attributes, ratios, strict_boundary = True, seed=None):
        """
        Initializes a CellArrangement object with the given parameters.
        Fills a volume with randomly placed nuclei of different attributes without intersection.

        Parameters:
            - mesh: mesh of volume to populate with nuclei
            - number: number of total nuclei to populate
            - attributes: list of nuclei type attributes that should appear
            - ratios: list of ratios of nuclei types to populate
            - strict_boundary: if true will place nuclei fully inside the mesh, if false only centroids will be placed fully inside the mesh
        """
        super().__init__()
        self.seed = seed
        self.name = "VolumeFill"
        self.mesh = mesh
        self.subdivision_levels = 2
        self.number = number
        self.attributes = attributes
        self.ratios = ratios
        self.strict_boundary = strict_boundary # If true will place nuclei fully inside the mesh, if false only centroids will be placed fully inside the mesh
        # Get count
        sum = np.sum(ratios)
        normalized_ratios = [ratio/sum for ratio in ratios]
        self.counts = [int(ratio*number) for ratio in normalized_ratios]
        # Generate points inside mesh with given minimum distance
        # TODO: Generate based on Mahalanobis distance for scaled spheres
        self.points_per_type = fill_volume(
            self.counts, self.attributes, self.mesh, self.strict_boundary, seed=seed)

    def add(self):
        for points, radius, type in self.points_per_type:
            self.add_nuclei(points, radius, type)

    def add_nuclei(self, locations, radius, type):
        attribute = self.get_attribute_by_type(self.attributes, type)

        # Create a small sphere object for each base point
        for idx, location in enumerate(locations):
            bpy.ops.mesh.primitive_ico_sphere_add(radius=radius)
            nucleus = bpy.context.active_object
            deform_mesh(nucleus, attribute)
            nucleus.name = f"Nucleus_Type_{type}_{idx}"
            nucleus.location = location
            nucleus.scale = attribute.scale
            # TODO: Apply a rotation to the volume nuclei such that they follow the tissue flow or smth like that I am no tissuologist how the hell should I know what that is I mean come on is there really someone who could answer that question yes maybe there is I will ask the Big Goose Of Enlightenment, quak.
            nucleus.rotation_euler = [random.uniform(0, 2*math.pi) for _ in range(3)]
            nucleus.modifiers.new(name="Subdivision", type='SUBSURF')
            nucleus.modifiers["Subdivision"].levels = self.subdivision_levels
            bpy.ops.object.modifier_apply({"object": nucleus}, modifier="Subdivision")
            self.objects.append(nucleus)

    def get_attribute_by_type(self, attributes, type):
        res = None
        for attribute in attributes:
            if attribute.cell_type == type:
                res = attribute
                break
        assert res is not None, "Type not found"
        return res

class SurfaceFill(CellArrangement):
    '''
    A surface fill arrangement. Given a surface mesh, fills the surface with randomly placed nuclei of different attributes.
    Also takes into account the maximal number of nuclei to place and their occuring ratios.
    '''
    def __init__(self, mesh, number, attribute, filler_scale):
        """
        Initializes a CellArrangement object with the given parameters.
        Fills a surface with randomly place nuclei of given attribute
        such that the x-axes of the nuclei are aligned with the surface normals.

        Parameters:
            - mesh: mesh of volume to populate with nuclei
            - number: number of total nuclei to populate
            - attributes: list of nuclei type attributes that should appear
            - ratios: list of ratios of nuclei types to populate
        """
        super().__init__()
        self.name = "SurfaceFill"
        self.mesh = mesh
        self.number = number
        self.attribute = attribute
        self.filler_scale = filler_scale
        self.subdivision_levels = 2
        self.main_verts, self.filler_verts, self.mesh_delta = fill_surface(self.mesh, self.number, self.attribute, self.filler_scale)


    def add(self):
        center_points = [v.co for v in self.main_verts]
        center_normals = [v.normal for v in self.main_verts]
        self.add_nuclei(center_points, center_normals, self.attribute)

        filler_points = [v.co for v in self.filler_verts]
        filler_normals = [v.normal for v in self.filler_verts]
        self.add_nuclei(filler_points, filler_normals, self.attribute, name="_small")

    def add_nuclei(self, points, normals, attribute, name=""):
        for idx, (pt, dir) in enumerate(zip(points, normals)):
            radius = attribute.radius if name == "" else attribute.radius * self.filler_scale
            bpy.ops.mesh.primitive_ico_sphere_add(radius=radius)
            nucleus = bpy.context.active_object
            # TODO: add deform. Is it more efficient to deform after adding? Can loop through all objects
            deform_mesh(nucleus, attribute)
            # NOTE: Random small displacement along the normal axis
            nucleus.location = pt + random.uniform(-0.5*self.mesh_delta, self.mesh_delta) * dir
            # Rotate obj such that local x axis is aligned with surface normal
            rotation_matrix = Matrix.Translation(nucleus.location) @ dir.to_track_quat('X').to_matrix().to_4x4()
            nucleus.matrix_world = rotation_matrix
            nucleus.scale = attribute.scale
            nucleus.name = f"Nucleus_Type_{attribute.cell_type}_{idx}{name}"
            nucleus.modifiers.new(name="Subdivision", type='SUBSURF')
            nucleus.modifiers["Subdivision"].levels = self.subdivision_levels
            bpy.ops.object.modifier_apply({"object": nucleus}, modifier="Subdivision")
            self.objects.append(nucleus)

NucleusSeed = namedtuple('NucleusSeed', 'centroid scale direction')
# NOTE: A NucleusSeed object is a 3d-ellipsoid defined by its 3d centroid, its scale
# (3d vector containing the length of axes in x, y and z direction), and its direction
# (3d normalized vector to which the x-axis will point).

class VoronoiFill(CellArrangement):
    def __init__(self, mesh_obj, count, attribute):
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
        self.attribute = attribute
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
        for idx, s in enumerate(self.nuclei_seeds):
            bpy.ops.mesh.primitive_ico_sphere_add()
            nucleus = bpy.context.active_object
            nucleus.name = f"Nucleus_Type_{self.attribute.cell_type}_{idx}"
            set_orientation(nucleus, s.direction)
            nucleus.location = s.centroid
            nucleus.scale = s.scale
            bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)
            deform_mesh(nucleus, self.attribute)
            bpy.ops.object.modifier_add(type='SUBSURF')
            bpy.context.object.modifiers["Subdivision"].levels = 3
            bpy.ops.object.modifier_apply(modifier="Subdivision")

            self.objects.append(nucleus)