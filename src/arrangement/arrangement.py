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
        self.padding_scale = 0.1
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
        # TODO: Check if ok, remove next line otherwise
        self.padding_scale = 1 / len(all_points)
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

    # TODO: Fusion this one with helper method add_nuclei_shaped(). - ck
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
            # Scale object to typical cell attribute size
            mean_scale = self.compute_mean_scale(nucleus_object)
            assert mean_scale > 0, "Nucleus object has a x, y or z-diameter of 0."
            nucleus_object.scale = tuple(x * size / mean_scale for x in scale)
            #nucleus_object.scale = (nuclei_scale,) * 3
            nucleus_objects.append(nucleus_object)
        return nucleus_objects
    
    def compute_mean_scale(self, object):
        '''
        Given a mesh object this method returns the mean scale of the object.
        The mean scale is defined as the geometric mean of the x, y and z-diameter of the mesh.      
        '''
        mesh = object.data
        min_coords = (min(vert.co[0] for vert in mesh.vertices),
                      min(vert.co[1] for vert in mesh.vertices),
                      min(vert.co[2] for vert in mesh.vertices))
        max_coords = (max(vert.co[0] for vert in mesh.vertices),
                      max(vert.co[1] for vert in mesh.vertices),
                      max(vert.co[2] for vert in mesh.vertices))
        diameter = [max - min for max, min in zip(max_coords, min_coords)]
        return (diameter[0]*diameter[1]*diameter[2]) ** (1/3)
    


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
        self.ico_xy_scale = param_dict["ico_xy_scale"] # Scale of the icosphere w.r.t. to the x-y-axes.
        self.z_rot_angle = param_dict["z_rot_angle"]  # Rotation along z-axis of crypt in degrees.
        self.center_loc = param_dict["center_loc"] # Center of crypt cut in world coordinates.

        self.nuclei_limit = 120  # If nuclei count is above this limit, abort script, due to long compute time.
        self.ico_z_scale = 2 # The smaller the number, the denser the nuclei. But it also increases the compute time.
        self.ico_scale = (*self.ico_xy_scale, self.ico_z_scale) # Scale of the icosphere
        self.slice_thickness = 0.1  # Reduce this if the computation time is too long.
        self.subdivision_levels = 1 # Set only to 1 or 2. 2 results in denser nuclei but larger compute time.
        self.tissue_cut_ratio = 1  # Determines which part of the slice should be cut by tissue. If this is less than 1, you get cut nuclei objects.
        self.inner_scale_coeff = 0.9 # Determines the size of the outer hull.
        self.outer_scale_coeff = 1.1 # Determines the size of the inner hull.
        self.random_translate = True # Enable or disable random displacement of nuclei centroids; True leads to more realistic results.
        self.max_translate = 0.02 # Maximal displacement of nuclei centroids.
        self.min_cut_box = [-8, -8, -self.slice_thickness/2] # Minimal coordinates of the slice cut box
        self.max_cut_box = [8, 8, self.slice_thickness/2] # Maximal coordinates of the slice cut box
        self.min_coords = [-8, -8, -2] # Minimal coordinates for auxiliary points.
        self.max_coords = [8, 8, 2] # Maximal coordinates for auxiliary points.
        self.padding = 0  # Padding of the auxiliary points box.
        self.region_scale = 1  # Scale of the Voronoi regions w.r.t. to the seed
        self.nuclei_scale = 1  # Scale of the nucleus object w.r.t. to the Voronoi region

        self.type = "EPIT"
        self.id = CellArrangement.count
        CellArrangement.count += 1
        self.objects = []

    def add(self):
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
        ico.name = "Surface_Med"
        inner_ico.name = "Surface_Inner"
        outer_ico.name = "Surface_Outer"


        # Generate seeds for Voronoi diagram
        bpy.context.view_layer.objects.active = ico
        bpy.ops.object.mode_set(mode='EDIT')
        mesh = bmesh.from_edit_mesh(ico.data)
        bmesh.ops.triangulate(mesh, faces=mesh.faces[:]) # Triangulate # NOTE: Check if this is truly necessary. - ck
        #print(f"Mesh has {len(mesh.faces)} faces and {len(mesh.verts)} vertices.")

        points = []
        for face in mesh.faces:
            centroid = face.calc_center_median()
            points.append(centroid)
            for v in face.verts:
                v_loc = ico.matrix_world @ v.co
                p = 0.5*(centroid + v_loc)
                points.append(p)
        #        q = 0.25*(centroid - v_loc) + centroid # NOTE: Can be used in case the nuclei are not dense enough. - ck
        #        points.append(q)
        for vert in mesh.verts:
            v_loc = ico.matrix_world @ vert.co
            points.append(v_loc)
        bpy.ops.object.mode_set(mode='OBJECT')
        # Retain only points that lie between min and max CUT_BOX
        # NOTE: Comment this line out if you want to get Voronoi cells on the whole surface.
        # WARNING: This could lead to very long computation times.
        points = [p for p in points if all(min_c <= val <= max_c for val, min_c, max_c in zip(p, self.min_cut_box, self.max_cut_box))]
        assert len(points) <= self.nuclei_limit, f"About to render {len(points)} nuclei. Stopped to avoid long compute time.\nYou can reduce the slice thickness to generate less nuclei." 
        if self.random_translate:
            points = [[p[i] + random.uniform(-1, 1)*self.max_translate for i in range(3)] for p in points]
        #add_point_cloud(points, radius = 0.01) # Render seed points


        # Add auxiliary boundary points to ensure that the base Voronoi regions are bounded.
        # The regions of the auxiliary points won't be.
        # TODO: Refactor get_lattice_points if padding is not needed - ck
        auxiliary_points = get_lattice_points(self.min_coords, self.max_coords)
        #add_point_cloud(auxiliary_points, radius = 0.2)

        # Generate the Voronoi diagram
        vor = Voronoi(points + auxiliary_points)
        #print_voronoi_stats(vor)
        fr_points = finite_region_points(vor)
        #print(f"Finite region points: {fr_points}")
        ridges = compute_faces_by_seeds(vor, fr_points)

        for _, point_idx in enumerate(fr_points):
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
        # Run the intersection function to get polytopes representing cell membranes
        region_objects = intersect_with_object(region_objects, outer_ico)
        region_objects = subtract_object(region_objects, inner_ico)
        box = add_box(self.min_cut_box, self.max_cut_box)
        region_objects = intersect_with_object(region_objects, box)

        # Turn each polytope in a mesh representing its nucleus
        nucleus_objects = add_nuclei_shaped(region_objects, self.nuclei_scale)
        # Remove too large arrtifacts
        artifacts = [] # NOTE: This lists too large regions
        for obj in nucleus_objects:
            # Calculate bounding box
            diameter = np.max([max(b) - min(b) for b in zip(*obj.bound_box)])
            if diameter > min(self.ico_scale): # NOTE: Maybe there is a better threshold. - ck
                artifacts.append(obj)
        nucleus_objects = [obj for obj in nucleus_objects if obj not in artifacts]

        # Intersect again
        box.scale = tuple(s*a for s,a in zip(box.scale, (1,1,self.tissue_cut_ratio)))
        nucleus_objects = intersect_with_object(nucleus_objects, box)
        # Create inner and outer hull
        inner_hull, outer_hull = intersect_with_object([inner_ico, outer_ico], box)
        crypt_objects = nucleus_objects + [inner_hull, outer_hull]
        # Transform crypt objects
        rotate_objects(crypt_objects, self.z_rot_angle)
        translate_objects(crypt_objects, self.center_loc)
        # Remove auxiliary objects
        remove_objects(artifacts)
        remove_objects([ico, box])
        remove_objects(region_objects)

