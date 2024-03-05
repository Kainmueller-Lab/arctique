import random

from mathutils import Vector

from src.objects.cells import Cell
from src.utils.geometry import *
from src.utils.helper_methods import *
from src.utils.surface_filling import fill_surface
from src.utils.voronoi import *

class CellArrangement:
    # Class variable to keep track of the count
    count = 0

    def __init__(self):
        self.objects = [] # Contains the nuclei objects
        self.auxiliary_objects = [] # These will be hidden in viewport and rendering but can be used for further computing.
        self.name = None
        self.id = CellArrangement.count
        CellArrangement.count += 1

class CellList(CellArrangement):
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
            # TODO: Make sure the ids are not just from the same range in different lists.
            # They need to be unique for each scene. - ck
            cell = Cell(location, self.id, self.type, self.cell_attributes)
            self.objects.append(cell)
        for cell in self.objects:
            cell.add()
        
# TODO: Add blowup algorithm (Monte Carlo)
class VolumeFill(CellArrangement):
    def __init__(self, mesh, number, attributes, ratios, strict_boundary = True):
        """
        Initializes a CellArrangement object with the given parameters.
        Fills a volume with randomly placed nuclei of different attributes without intersection.

        Parameters:
            - mesh: mesh of volume to populate with nuclei
            - number: number of total nuclei to populate
            - attributes: list of nuclei type attributes that should appear
            - ratios: list of ratios of nuclei types to populate
        """
        super().__init__()
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
        self.points_per_type = generate_points_per_type(self.counts, self.attributes, self.mesh, self.strict_boundary)

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
            # TODO: orientation of nuclei aligned to mesh
            nucleus.rotation_euler = [random.uniform(0, 2*math.pi) for _ in range(3)]
            # TODO: Necessary?
            # Make nuclei children of bounding mesh
            #nucleus.parent = self.mesh
            #world_inv = self.mesh.matrix_world.inverted()
            #nucleus.matrix_parent_inverse = world_inv
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
        self.main_verts, self.filler_verts = fill_surface(self.mesh, self.number, self.attribute, self.filler_scale)


    def add(self):
        center_points = [v.co for v in self.main_verts]
        center_normals = [v.normal for v in self.main_verts]
        self.add_nuclei(center_points, center_normals, self.attribute)

        filler_points = [v.co for v in self.filler_verts]
        filler_normals = [v.normal for v in self.filler_verts]
        self.add_nuclei(filler_points, filler_normals, self.attribute, name="_small")

    def add_nuclei(self, points, normals, attribute, name=""):
        for idx, (pt, dir) in enumerate(zip(points, normals)):
            # TODO: 0.75 is hard coded here. Fix that. - ck
            radius = attribute.size if name == "" else attribute.size * self.filler_scale
            bpy.ops.mesh.primitive_uv_sphere_add(radius=radius)
            nucleus = bpy.context.active_object
            # TODO: add deform
            #deform_mesh(nucleus, attribute)
            nucleus.location = pt
            # Rotate obj such that local x axis is aligned with surface normal
            rotation_matrix = Matrix.Translation(nucleus.location) @ dir.to_track_quat('X').to_matrix().to_4x4()
            nucleus.matrix_world = rotation_matrix
            nucleus.scale = attribute.scale
            nucleus.name = f"Surface_Nucleus_Type_{attribute.cell_type}_{idx}{name}"
            self.objects.append(nucleus)