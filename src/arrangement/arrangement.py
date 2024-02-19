import random

from mathutils import Vector

from src.objects.cells import Cell
from src.utils.geometry import *
from src.utils.helper_methods import *
from src.utils.voronoi import *

class CellArrangement:
    # Class variable to keep track of the count
    count = 0

    def __init__(self):
        self.objects = [] # Contains the nuclei objects
        self.auxiliary_objects = [] # These will be hidden in viewport and rendering but can be used for further computing.
        self.id = CellArrangement.count
        CellArrangement.count += 1

    def add(self):
        pass

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
        self.type = "LIST"

    def add(self):
        for location in self.locations:
            # TODO: Make sure the ids are not just from the same range in different lists.
            # They need to be unique for each scene. - ck
            cell = Cell(location, self.id, self.type, self.cell_attributes)
            self.objects.append(cell)
        for cell in self.objects:
            cell.add()
        

# TODO: Extend to different distribution types, e.g., Gaussian distribution, etc.
class CellDistribution(CellArrangement):
    def __init__(self, cell_attributes, num_cells, min_coords, max_coords):
        """
        Initializes a CellArrangement object with the given parameters.

        Parameters:
            cell_attributes (list): A list of cell attributes.
            num_cells (int): The number of cells.
            min_coords (tuple): The minimum coordinates.
            max_coords (tuple): The maximum coordinates.
        """
        super().__init__()
        self.cell_attributes = cell_attributes
        self.num_cells = num_cells
        self.min_coords = min_coords
        self.max_coords = max_coords
        self.type = "DIST"

    def add(self):
        """
        Generates cells and adds them to the list of objects.

        This function generates a specified number of cells and adds them to the list of objects. Each cell is randomly located within the specified minimum and maximum coordinates. The orientation of each cell is also randomly generated. 
        """
        for _ in range(self.num_cells):
            # Sample uniformly distributed location
            location = Vector([
                random.uniform(self.min_coords.x, self.max_coords.x),
                random.uniform(self.min_coords.y, self.max_coords.y),
                random.uniform(self.min_coords.z, self.max_coords.z)
            ])
            orientation = random_unit_vector()
            orientation = Vector([*orientation])
            cell = Cell(location, self.id, self.type, self.cell_attributes, orientation)
            self.objects.append(cell)
        for cell in self.objects:
            cell.add()


class VoronoiDiagram(CellArrangement):
    def __init__(self, distribution_dict):
        """
        Initializes a CellArrangement object with the given parameters.

        Parameters:
            distribution_dict (dict): Maps cell attribute types to a list of corresp. 3D points
        """
        super().__init__()
        self.distribution_dict = distribution_dict
        # TODO: Make these parameters changeable and dependent of cell attribute type
        self.region_scale = 1
        self.nuclei_scale = 0.5
        self.padding_scale = 0.1
        self.type = "VORO"
        self.empty_regions = [] # Contains regions (=meshes) that should not be populated by random nuclei.

    def add(self):
        # Generate seed points and auxiliary points for Voronoi regions
        # Auxiliary boundary points ensure that the base Voronoi regions are bounded.
        # The regions of the auxiliary points won't be bounded.
        all_points = []
        for point_list in self.distribution_dict.values():
            all_points.extend(point_list)
        min_coords = Vector((min(point[0] for point in all_points),
                               min(point[1] for point in all_points),
                               min(point[2] for point in all_points)))
        max_coords = Vector((max(point[0] for point in all_points),
                               max(point[1] for point in all_points),
                               max(point[2] for point in all_points)))
        self.padding_scale = 1 / len(all_points) # TODO: Check if ok, remove next line otherwise
        padding = (max_coords - min_coords) * self.padding_scale
        auxiliary_points = get_lattice_points(min_coords - padding, max_coords + padding)
        #add_point_cloud(auxiliary_points, radius = 0.2) # NOTE: Only for visualization

        # Add Voronoi regions and shape them via intersection
        generate_voronoi_regions(all_points, auxiliary_points, self.region_scale)
        region_objects = get_objects_with("CellObject")
        box_object = add_box(min_coords-padding, max_coords+padding)
        region_objects = intersect_with_object(region_objects, box_object)

        # Add nuclei objects based on the finite Voronoi regions
        all_attributes = []
        # Collect list of attributes per cell seed
        # NOTE: This is quite hacky. Could be done more elegantly. - ck
        for attribute in self.distribution_dict.keys():
            all_attributes.extend([attribute for _ in self.distribution_dict[attribute]])
        assert len(all_attributes)==len(all_points), "Total number of attributes not matching with number of cell seeds."
        self.objects = self.generate_nucleus_objects(region_objects, all_attributes)

        # Remove auxiliary objects
        remove_objects(region_objects + [box_object])
        # Remove nuclei that intersect or lie within the outer hulls
        self.objects = remove_objects_inside_mesh(self.objects, self.empty_regions)
        # print("Created Voronoi arrangement:")
        # for attribute in self.distribution_dict.keys():
        #     print(f"- {len(self.distribution_dict[attribute])} nuclei of type {attribute.cell_type}")

    def generate_nucleus_objects(self, region_objects, attributes):
        nucleus_objects = shrinkwrap(region_objects)
        for index, nucleus_object in enumerate(nucleus_objects):
            attribute = attributes[index]
            size = attribute.size
            scale = attribute.scale
            cell_type = attribute.cell_type
            nucleus_object.name = f"NucleusObject_{index}_{self.type}_{self.count}_Type_{cell_type}"
            # Scale object to typical cell attribute size
            mean_scale = compute_mean_scale(nucleus_object)
            assert mean_scale > 0, "Nucleus object has a x, y or z-diameter of 0."
            nucleus_object.scale = tuple(x * size / mean_scale for x in scale)
        return nucleus_objects
        
    

class EpithelialArrangement(CellArrangement):
    def __init__(self, param_dict):
        """
        Initializes a CellArrangement object with the given parameters.

        Parameters:
            - param_dict (dict): A dictionary containing the following keys:
                - ico_xy_scale (tuple): The scale of the icosphere w.r.t. to the x-y-axes.
                - z_rot_angle (float): The rotation along the z-axis of the crypt in degrees.
                - center_loc (tuple): The center of the crypt cut in world coordinates.
        """
        super().__init__()
        self.ico_xy_scale = param_dict["ico_xy_scale"] # Scale of the icosphere w.r.t. to the x-y-axes.
        self.z_rot_angle = param_dict["z_rot_angle"]  # Rotation along z-axis of crypt in degrees.
        self.center_loc = param_dict["center_loc"] # Center of crypt cut in world coordinates.

        self.nuclei_limit = 90  # If nuclei count is above this limit, abort script, due to long compute time.
        self.ico_z_scale = 2 # The smaller the number, the denser the nuclei. But it also increases the compute time.
        self.ico_scale = (*self.ico_xy_scale, self.ico_z_scale) # Scale of the icosphere
        self.slice_thickness = 0.1  # Reduce this if the computation time is too long.
        self.subdivision_levels = 1 # Set only to 1 or 2. 2 results in denser nuclei but larger compute time.
        self.tissue_cut_ratio = 1  # Determines which part of the slice should be cut by tissue. If this is less than 1, you get cut nuclei objects.
        self.inner_scale_coeff = 0.9 # Determines the size of the outer hull.
        self.outer_scale_coeff = 1.1 # Determines the size of the inner hull.
        self.max_translate = 0.02 # Maximal displacement of nuclei centroids. Set to None if you want no random translation
        self.min_cut_box = [-8, -8, -self.slice_thickness/2] # Minimal coordinates of the slice cut box
        self.max_cut_box = [8, 8, self.slice_thickness/2] # Maximal coordinates of the slice cut box
        self.min_coords = [-8, -8, -2] # Minimal coordinates for auxiliary points.
        self.max_coords = [8, 8, 2] # Maximal coordinates for auxiliary points.
        self.padding = 0  # Padding of the auxiliary points box.
        self.region_scale = 1  # Scale of the Voronoi regions w.r.t. to the seed
        self.nuclei_scale = 1  # Scale of the nucleus object w.r.t. to the Voronoi region

        self.type = "EPIT"
        self.inner_hull = None
        self.outer_hull = None

    def add(self):
        # Add icospheres whose surface will be used for the nuclei distribution
        ico, inner_ico, outer_ico = self.add_spheres()
        # Sample seeds for Voronoi regions on the median icosphere
        cut_box = (self.min_cut_box, self.max_cut_box)
        points = sample_points_on_mesh(ico, self.nuclei_limit, cut_box, self.max_translate)
        # add_point_cloud(points, radius = 0.01) # Render seed points
        # Add auxiliary boundary points to ensure that the base Voronoi regions are bounded. The regions of the auxiliary points won't be.
        # TODO: Refactor get_lattice_points if padding is not needed - ck
        auxiliary_points = get_lattice_points(self.min_coords, self.max_coords)
        #add_point_cloud(auxiliary_points, radius = 0.2)

        # Generate Voronoi regions and add nucleus objects
        generate_voronoi_regions(points, auxiliary_points, self.region_scale)
        region_objects = get_objects_with("CellObject")
        region_objects = intersect_with_object(region_objects, outer_ico)
        region_objects = subtract_object(region_objects, inner_ico)
        box = add_box(self.min_cut_box, self.max_cut_box)
        region_objects = intersect_with_object(region_objects, box)
        nucleus_objects, artifacts = self.generate_nucleus_objects(region_objects)
        box.scale = tuple(s*a for s,a in zip(box.scale, (1,1,self.tissue_cut_ratio)))
        self.objects = intersect_with_object(nucleus_objects, box)
        # Create inner and outer hull
        self.inner_hull, self.outer_hull = intersect_with_object([inner_ico, outer_ico], box)
        self.inner_hull.name = f"HullObject_Inner_{self.type}_{self.count}"
        self.outer_hull.name = f"HullObject_Outer_{self.type}_{self.count}"
        self.auxiliary_objects = [self.inner_hull, self.outer_hull]
        crypt_objects = self.objects + self.auxiliary_objects
        # Transform crypt objects
        rotate_objects(crypt_objects, self.z_rot_angle)
        translate_objects(crypt_objects, self.center_loc)
        # Remove auxiliary objects
        remove_objects(artifacts)
        remove_objects([ico, box])
        remove_objects(region_objects)

    def add_spheres(self):
        # Create icosphere and subdivide it to desired level
        # NOTE: Do not apply higher level subdiv than 2, it crashes blender.
        bpy.ops.mesh.primitive_ico_sphere_add(enter_editmode=False, align='WORLD', location=(0,0,0), scale=self.ico_scale)
        ico = bpy.context.active_object
        ico.modifiers.new(name="Subdivision", type='SUBSURF')
        ico.modifiers["Subdivision"].levels = self.subdivision_levels
        bpy.ops.object.modifier_apply({"object": ico}, modifier="Subdivision")

        # Create inner and outer ico spheres
        bpy.ops.object.duplicate(linked=False)
        inner_ico = bpy.context.active_object
        bpy.ops.object.duplicate(linked=False)
        outer_ico = bpy.context.active_object
        inner_ico.scale = tuple(x*self.inner_scale_coeff for x in ico.scale)
        outer_ico.scale = tuple(x*self.outer_scale_coeff for x in ico.scale)
        return ico, inner_ico, outer_ico

    def generate_nucleus_objects(self, region_objects):
        # Turn each polytope in a mesh representing its nucleus
        nucleus_objects = shrinkwrap(region_objects, self.nuclei_scale)
        for idx, obj in enumerate(nucleus_objects):
            obj.name = f"NucleusObject_{idx}_{self.type}_{self.count}"
        # Remove too large arrtifacts
        artifacts = [] # NOTE: This lists too large regions
        for obj in nucleus_objects:
            # Calculate bounding box
            diameter = np.max([max(b) - min(b) for b in zip(*obj.bound_box)])
            if diameter > min(self.ico_scale): # NOTE: Maybe there is a better threshold. - ck
                artifacts.append(obj)
        nucleus_objects = [obj for obj in nucleus_objects if obj not in artifacts]
        return nucleus_objects, artifacts
    

class VolumeFill(CellArrangement):
    def __init__(self, mesh, number, attributes, ratios):
        """
        Initializes a CellArrangement object with the given parameters.
        Fills a volume with randomly place nuclei of different attributes without intersection.

        Parameters:
            - mesh: mesh of volume to populate with nuclei
            - number: number of total nuclei to populate
            - attributes: list of nuclei type attributes that should appear
            - ratios: list of ratios of nuclei types to populate
        """
        super().__init__()
        self.mesh = mesh
        self.number = number
        self.attributes = attributes
        self.ratios = ratios
        self.radii = None # TODO
        self.types = None # TODO

        # Get counts
        sum = sum(ratios)
        normalized_ratios = [ratio/sum for ratio in ratios]
        self.counts = [int(ratio*number) for ratio in normalized_ratios]

        # Generate points inside mesh with given minimum distance
        self.points_per_type = generate_points_per_type(self.counts, self.radii, self.types, self.mesh)


    def add(self):
        for points, radius, type in self.points_per_type:
            # TODO: Make general add_nuclei method
            add_point_cloud(points, radius, type)

