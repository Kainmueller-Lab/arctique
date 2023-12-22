import random
from mathutils import Vector
from src.objects.cells import Cell
from src.utils.helper_methods import *
from src.utils.geometry import * 

class CellArrangement:
    # Class variable to keep track of the count
    count = 0

    def __init__(self):
        pass

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
        self.cell_attributes = cell_attributes
        self.locations = locations
        self.type = "LIST"
        self.id = CellArrangement.count
        CellArrangement.count += 1
        self.objects = []

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
        self.cell_attributes = cell_attributes
        self.num_cells = num_cells
        self.min_coords = min_coords
        self.max_coords = max_coords
        self.type = "DIST"
        self.id = CellArrangement.count
        CellArrangement.count += 1
        self.objects = []

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
    def __init__(self, distribution_dict, min_coords, max_coords):
        """
        Initializes a CellArrangement object with the given parameters.

        Parameters:
            distribution_dict (dict): Lists number of cells to add per cell attribute type
            min_coords (tuple): The minimum coordinates.
            max_coords (tuple): The maximum coordinates.
        """
        self.distribution_dict = distribution_dict
        self.cell_count = sum(distribution_dict.values())
        self.min_coords = min_coords
        self.max_coords = max_coords
        # TODO: Make these parameters changeable and dependent of cell attribute type
        self.region_scale = 1
        self.nuclei_scale = 0.5
        self.type = "VORO"
        self.id = CellArrangement.count
        CellArrangement.count += 1
        self.objects = []

    def add(self):
        # Generate N random 3D vectors bewteen min and max coords
        points = [list(map(random.uniform, self.min_coords, self.max_coords)) for _ in range(self.cell_count)]
        #add_point_cloud(points, radius = 0.01)

        # Add auxiliary boundary points to ensure that the base Voronoi regions are bounded
        # The regions of the auxiliary points won't be.
        # NOTE: You need to choose one of those three
        #auxiliary_points = get_octogon_points(self.min_coords, self.max_coords, padding=0.5)
        #auxiliary_points = get_cube_points(self.min_coords, self.max_coords, padding=0.5)
        auxiliary_points = get_lattice_points(self.min_coords, self.max_coords)
        #add_point_cloud(auxiliary_points, radius = 0.2)

        # Generate the Voronoi diagram and necessary data
        vor = Voronoi(points + auxiliary_points)
        #print_voronoi_stats(vor)
        fr_points = finite_region_points(vor)
        #print(f"Finite region points: {fr_points}")
        ridges = compute_faces_by_seeds(vor, fr_points)

        for point_idx in fr_points:
            add_region(vertices = vor.vertices,
                            faces = ridges[point_idx],
                            idx = point_idx)
            # Set the 3D cursor to the desired position
            bpy.context.scene.cursor.location = points[point_idx]
            # Set the object's origin to the 3D cursor location
            bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
            # Scale the object (adjust the scale factors as needed)
            scale = self.region_scale
            bpy.ops.transform.resize(value=(scale, scale, scale), orient_type='LOCAL')
            obj = bpy.context.active_object
            obj.select_set(False)

        # Collect region objects in the scene
        region_objects = get_objects_with("CellObject")
        # Create bounding box of distribution
        box_object = add_box(self.min_coords, self.max_coords)
        # Run the intersection function to get polytopes representing cell membranes
        cell_objects = intersect_with_object(region_objects, box_object)
        # Turn each polytope in a mesh representing its nucleus
        self.objects = add_nuclei_from(cell_objects, self.nuclei_scale)
        remove_objects(cell_objects)
