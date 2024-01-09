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
    def __init__(self, distribution_dict):
        """
        Initializes a CellArrangement object with the given parameters.

        Parameters:
            distribution_dict (dict): Maps cell attribute types to a list of corresp. 3D points
        """
        self.distribution_dict = distribution_dict
        # TODO: Make these parameters changeable and dependent of cell attribute type
        self.region_scale = 1
        self.nuclei_scale = 0.5
        self.padding_scale = 0.1 # TODO: Adapt to cellcount and box size
        self.type = "VORO"
        self.id = CellArrangement.count
        CellArrangement.count += 1
        self.objects = []

    def add(self):
        # Collect all nuclei center points
        all_points = []
        for point_list in self.distribution_dict.values():
            all_points.extend(point_list)

        # Compute auxiliary boundary points and padding
        min_coords = Vector((min(point[0] for point in all_points),
                               min(point[1] for point in all_points),
                               min(point[2] for point in all_points)))
        max_coords = Vector((max(point[0] for point in all_points),
                               max(point[1] for point in all_points),
                               max(point[2] for point in all_points)))
        padding = (max_coords - min_coords) * self.padding_scale

        # Add auxiliary boundary points to ensure that the base Voronoi regions are bounded
        # The regions of the auxiliary points won't be bounded.
        # NOTE: You need to choose one of those three
        #auxiliary_points = get_octogon_points(self.min_coords, self.max_coords, padding=0.5)
        #auxiliary_points = get_cube_points(self.min_coords, self.max_coords, padding=0.5)
        auxiliary_points = get_lattice_points(min_coords - padding, max_coords + padding)
        #add_point_cloud(auxiliary_points, radius = 0.2)

        # Generate the Voronoi diagram and necessary data
        vor = Voronoi(all_points + auxiliary_points)
        fr_points = finite_region_points(vor)
        if not len(fr_points)==len(all_points):
            print("Less nuclei than expected have been generated.")
        assert len(fr_points)==len(all_points), "Less nuclei than expected have been generated."
        ridges = compute_faces_by_seeds(vor, fr_points)

        for point_idx in fr_points:
            add_region(vertices = vor.vertices,
                            faces = ridges[point_idx],
                            idx = point_idx)
            # Set the 3D cursor to the desired position
            bpy.context.scene.cursor.location = all_points[point_idx]
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
        box_object = add_box(min_coords-padding, max_coords+padding)
        # Run the intersection function to get polytopes representing cell membranes
        cell_objects = intersect_with_object(region_objects, box_object)
        # Collect list of attributes per cell seed
        # NOTE: This is quite hacky. Could be done more elegantly. - ck
        all_attributes = []
        for attribute in self.distribution_dict.keys():
            all_attributes.extend([attribute for _ in self.distribution_dict[attribute]])
        assert len(all_attributes)==len(all_points), "Total number of attributes not matching with number of cell seeds."
        # Turn each polytope in a mesh representing its nucleus
        self.objects = self.add_nuclei_from(cell_objects, all_attributes)
        remove_objects(cell_objects)
        print("Created Voronoi arrangement:")
        for attribute in self.distribution_dict.keys():
            print(f"- {len(self.distribution_dict[attribute])} nuclei of type {attribute.cell_type}")

    def add_nuclei_from(self, cell_objects, attributes):
        nucleus_objects = []
        for cell_idx, cell_object in enumerate(cell_objects):
            attribute = attributes[cell_idx]
            size = attribute.size
            scale = attribute.scale
            cell_type = attribute.cell_type

            bpy.ops.mesh.primitive_cube_add(enter_editmode=False, align='WORLD', location=cell_object.location, scale=(1, 1, 1))
            nucleus_object = bpy.context.active_object
            index = cell_object.name.split('_')[1]
            nucleus_object.name = f"NucleusObject_{index}_{self.type}_Type_{cell_type}"
            shrinkwrap = nucleus_object.modifiers.new(name="Shrinkwrap Modifier", type='SHRINKWRAP')
            shrinkwrap.target = cell_object
            bpy.ops.object.modifier_apply(modifier="Shrinkwrap Modifier")
            subsurf = nucleus_object.modifiers.new("Subsurface Modifier", type='SUBSURF')
            subsurf.levels = 2
            bpy.ops.object.modifier_apply(modifier="Subsurface Modifier")
            nucleus_object.scale = tuple(x * size for x in scale)
            #nucleus_object.scale = (nuclei_scale,) * 3
            nucleus_objects.append(nucleus_object)
        return nucleus_objects
