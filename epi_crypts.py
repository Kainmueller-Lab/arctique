import bmesh
import bpy
import numpy as np
import random

from math import radians, sin, cos, pi
from mathutils import Matrix, Vector, geometry, Quaternion
from scipy.spatial import Voronoi # NOTE: You might install it directly into Blender's python using pip. - ck
from scipy.spatial.transform import Rotation

# NOTES about current version: - ck
# - This script models tissue cuts through epithelial crypts.
# - The nuclei are created from a Voronoi diagram whose seeds are placed onto the mesh of a scaled icosphere.
# - The seed points are not chosen randomly but set as the vertices and face centers of the mesh.
# - The resulting Voronoi regions are intersected with an inner and outer icosphere to create small nuclei objects along the surface.
# - In order to reduce the computing time not the whole surface and nuclei are rendered but only a slice of changeable thickness.
# - The resulting objects form a cut though such a crypt with nuclei at its boundary.
# - In addition two hull objects a re created that can be used to define the inside and outside of the crypt.
# - It is possible to manually set the scale, location and rotation in xy-plane of the resulting crypt cut.



################ HELPERS ####################
def clear_scene():
    bpy.context.view_layer.objects.active = bpy.context.scene.objects[0]
    bpy.ops.object.mode_set(mode='OBJECT')
    # Remove all objects
    bpy.ops.object.select_all(action='DESELECT')
    bpy.ops.object.select_by_type(type='MESH')
    bpy.ops.object.delete()
    # Remove all materials
    for mat in bpy.data.materials:
        bpy.data.materials.remove(mat, do_unlink=True)
    # Remove orphaned meshes
    for mesh in bpy.data.meshes:
        bpy.data.meshes.remove(mesh, do_unlink=True)

def print_voronoi_stats(vor):    
    # Vertex locations
    print(f"Vertices: {vor.vertices}")
    # - The Voronoi vertices represent the points where three or more Voronoi cells meet.
    # - Each vertex is a 3D point in the input space.

    # Ridge vertices
    print(f"Ridge vertices: {vor.ridge_vertices}")
    # - The ridge vertices represent the vertices of the Voronoi ridges (faces between Voronoi cells).
    # - Each ridge is represented by a set of vertex indices.
    # - If one of the indices is -1, the ridge extends to infinity.

    # Ridge dictionary
    print(f"Ridge dict: {vor.ridge_dict}")
    # - The ridge dictionary maps pairs of point indices to the corresponding ridge index.
    # - It provides a way to look up the ridge index based on the indices of the points it connects.

    # Points
    print(f"Regions per point: {vor.point_region}")
    # - The point_region array maps each input point to the index of the Voronoi region (cell) it belongs to.
    # - The index corresponds to the index of the Voronoi region in the 'regions' array.

    # Ridge points
    print(f"Ridge points: {vor.ridge_points}")
    # - The ridge_points array provides the indices of the points that form the corners of each Voronoi ridge.
    # - Each row corresponds to a ridge, and the two columns contain the indices of the points.

    # Regions
    print(f"Regions: {vor.regions}")
    # - The regions array contains the vertices of each Voronoi region (polygonal cell).
    # - It includes the indices of the Voronoi vertices that form the boundary of the region.
    # - Some entries may be empty lists, representing unbounded regions.

def add_point_cloud(locations, radius):
    # Create a small sphere object for each base point
    for idx, location in enumerate(locations):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=location)
        sphere = bpy.context.active_object
        sphere.name = f"Sphere_{idx}"

def finite_region_points(vor):
    fr_points = []
    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if -1 not in region and len(region) > 0:
            fr_points.append(point_idx)
    return fr_points


def voronoi_test():
    # Generate 27 integer points in the 0, -1, 1 lattice
    points = [[i, j, k] for i in range(-1, 2) for j in range(-1, 2) for k in range(-1, 2)]
    vor = Voronoi(points)
    fr_points = finite_region_points(vor)
    assert len(fr_points) == 1
    assert fr_points[0] == 13

def add_region(vertices, faces, idx = 0):
    '''
    Creates a 3D polytope mesh in the scene.
    vertices: a list of 3D positions
    faces: a list of lists of vertex indices, one list of vertex indices cocrresponding to a face
    '''
    mesh = bpy.data.meshes.new(name=f"CellMesh_{idx}")
    obj = bpy.data.objects.new(f"CellObject_{idx}", mesh)
    bpy.context.scene.collection.objects.link(obj)
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)

    # Get a BMesh to create the mesh
    bm = bmesh.new()
    # Add vertices
    for vertex in vertices:
        bm.verts.new(vertex)
    bm.verts.new([0, 0, 10]) # Infinite vertex
    bm.verts.ensure_lookup_table() # NOTE: Necessary to allow vertex indexing
    # Add faces
    for face_indices in faces:
        bm.faces.new([bm.verts[i] for i in face_indices])
    # Update the mesh with the new geometry
    bm.to_mesh(mesh)
    bm.free()


def get_cube_points(min_coords, max_coords, padding):
    points = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                point = [max_coords[0]+padding if i else min_coords[0]-padding,
                        max_coords[1]+padding if j else min_coords[1]-padding,
                        max_coords[2]+padding if k else min_coords[2]-padding]
                points.append(point)
    return points

def get_octogon_points(min_coords, max_coords, padding):
    points = []
    cube_center = [0.5 * (i + j) for i, j in zip(min_coords, max_coords)]
    for i in range(3):
        point1, point2 = cube_center.copy(), cube_center.copy()
        point1[i] = max_coords[i] + padding
        points.append(point1)
        point2[i] = min_coords[i] - padding
        points.append(point2)
    return points

def get_lattice_points(min_coords, max_coords, padding=0):
    points = []
    cube_center = Vector([0.5 * (i + j) for i, j in zip(min_coords, max_coords)])
    diameter = Vector([j - i for i, j in zip(min_coords, max_coords)])
    diameter += diameter * padding
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for k in [-1, 0, 1]:
                if not (i==0 and j==0 and k==0):
                    point = cube_center + Vector([a*b for (a,b) in zip(diameter, [i,j,k])])
                    points.append(point)
    return points

def add_box(min_coords, max_coords):
    bpy.ops.mesh.primitive_cube_add(size=1, enter_editmode=False, align='WORLD', location=((max_coords[0]+min_coords[0])/2, (max_coords[1]+min_coords[1])/2, (max_coords[2]+min_coords[2])/2))
    bpy.context.active_object.dimensions = [(max_coords[i] - min_coords[i]) for i in range(3)]
    bpy.context.active_object.location = [(max_coords[i] + min_coords[i]) / 2 for i in range(3)]
    return bpy.context.active_object

def intersect_with_object(target_objects, intersect_object):
    # Iterate through each target object
    for target_object in target_objects:
        boolean = target_object.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
        boolean.operation = 'INTERSECT'
        boolean.object = intersect_object

        # Apply the boolean modifier and remove the cube object
        bpy.ops.object.modifier_apply({"object": target_object}, modifier="Boolean Modifier")
    return target_objects

def subtract_object(target_objects, subtract_object):
    # Iterate through each target object
    for target_object in target_objects:
        boolean = target_object.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
        boolean.operation = 'DIFFERENCE'
        boolean.use_self = True
        boolean.object = subtract_object

        # Apply the boolean modifier and remove the cube object
        bpy.ops.object.modifier_apply({"object": target_object}, modifier="Boolean Modifier")
    return target_objects

def add_nuclei_shaped(cell_objects, nuclei_scale):
    nucleus_objects = []
    for cell_object in cell_objects:
        bpy.ops.mesh.primitive_cube_add(enter_editmode=False, align='WORLD', location=cell_object.location)
        nucleus_object = bpy.context.active_object
        index = cell_object.name.split('_')[1]
        nucleus_object.name = f"NucleusObject_{index}"
        shrinkwrap = nucleus_object.modifiers.new(name="Shrinkwrap Modifier", type='SHRINKWRAP')
        shrinkwrap.target = cell_object
        bpy.ops.object.modifier_apply(modifier="Shrinkwrap Modifier")
        subsurf = nucleus_object.modifiers.new("Subsurface Modifier", type='SUBSURF')
        subsurf.levels = 2
        bpy.ops.object.modifier_apply(modifier="Subsurface Modifier")
        nucleus_object.scale = (nuclei_scale, nuclei_scale, nuclei_scale)
        nucleus_objects.append(nucleus_object)
    return nucleus_objects
        
def remove_objects(object_list):
    for obj in object_list:
        bpy.data.objects.remove(obj, do_unlink=True)

def get_objects_with(string):
    # Create a list to store matching objects
    object_list = []
    # Iterate through all objects
    for obj in bpy.data.objects:
        # Check if the object name starts with the defined pattern
        if obj.name.startswith(string):
            object_list.append(obj)
    return object_list

def compute_faces_by_seeds(vor, seeds):
    ridges = {}
    # Print ridge pairs of finite regions
    for point_idx in seeds:
        # Get boundary ridges (= faces, as list of face vertices)
        ridge_list = []
        for ridge_pair in vor.ridge_dict.keys():
            if point_idx in set(ridge_pair):
                ridge = vor.ridge_dict[ridge_pair]
                #print(f"Boundary ridge pair: {ridge_pair} - {ridge}")
                ridge_list.append(ridge)
        ridges[point_idx] = ridge_list
    print(f"Created the Voronoi region faces")
    return ridges


def is_point_inside_mesh(mesh, point):
    # Get the mesh data
    mesh_data = mesh.data
    # Transform the point to the mesh's local coordinates
    local_point = mesh.matrix_world.inverted() @ point
    # Check if the point is inside the mesh
    is_inside = geometry.intersect_point_triangles(local_point, mesh_data.vertices, mesh_data.polygons)
    return is_inside




def rotate_3d_points_(points, alpha):
    # Convert angle to radians
    alpha_rad = np.radians(alpha)
    # Create a Rotation object for the z-axis rotation
    rotation = Rotation.from_euler('z', alpha_rad)
    # Apply the rotation to all points
    rotated_points_array = rotation.apply(points)
    # Convert back to a list of tuples
    rotated_points = [list(p) for p in rotated_points_array]
    return rotated_points

def rotate_objects(objects, alpha):
    # Convert angle to radians
    alpha_rad = radians(alpha)
    # Select the objects and activate the context
    bpy.context.view_layer.objects.active = objects[0]
    bpy.ops.object.select_all(action='DESELECT')
    for obj in objects:
        obj.select_set(True)
    # Apply the rotation using bpy.ops.transform.rotate
    bpy.ops.transform.rotate(value=alpha_rad, orient_axis='Z', orient_type='GLOBAL',
                             orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                             orient_matrix_type='GLOBAL', constraint_axis=(False, False, True),
                             mirror=False, use_proportional_edit=False, proportional_edit_falloff='SMOOTH',
                             proportional_size=1, use_proportional_connected=False, use_proportional_projected=False,
                             snap=False, snap_elements={'INCREMENT'}, use_snap_project=False, snap_target='CLOSEST',
                             use_snap_self=True, use_snap_edit=True, use_snap_nonedit=True, use_snap_selectable=False)
    bpy.ops.object.select_all(action='DESELECT')
    return  

def translate_objects(objects, location):
    # Select the objects and activate the context
    bpy.context.view_layer.objects.active = objects[0]
    bpy.ops.object.select_all(action='DESELECT')
    for obj in objects:
        obj.select_set(True)
    # Apply translation using bpy.ops.transform.translate
    bpy.ops.transform.translate(value=location, orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                                orient_matrix_type='GLOBAL', constraint_axis=(False, False, True),
                                mirror=False, use_proportional_edit=False, proportional_edit_falloff='SMOOTH',
                                proportional_size=1, use_proportional_connected=False, use_proportional_projected=False,
                                snap=False, snap_elements={'INCREMENT'}, use_snap_project=False, snap_target='CLOSEST',
                                use_snap_self=True, use_snap_edit=True, use_snap_nonedit=True, use_snap_selectable=False)
    # Deselect all objects
    bpy.ops.object.select_all(action='DESELECT')
    return                 


################# PARAMETERS ####################

NUCLEI_LIMIT = 120 # If nuclei count is above this limit stop, due to long compute time
SLICE_THICKNESS = 0.1 # Reduce this if the computation time is too long.
SUBDIVISION_LEVELS = 1
TISSUE_CUT_RATIO = 1 # Determines which part of the slice should be cut by tissue. If this is less than 1 you get cut nuclei objects.

ICO_SCALE = (0.3, 0.1, 1)
INNER_SCALE_COEFF = 0.9
OUTER_SCALE_COEFF = 1.1

RANDOM_TRANSLATE = True
MAX_TRANSLATE = 0.02
Z_ROT_ANGLE = 40 # Rotation along z-axis of crypt in degrees
CENTER_LOC = (0,1,0)


MIN_CUT_BOX = [-8, -8, 0]
MAX_CUT_BOX = [8, 8, SLICE_THICKNESS]

MIN_COORDS = [-8,-8,-2]
MAX_COORDS = [8,8,2]
REGION_SCALE = 1 # Scale of the Voronoi regions w.r.t. to the seed
NUCLEI_SCALE = 1 # Scale of the nucleus object w.r.t. to the Voronoi region
PADDING = 0 # NOTE: Can be removed

# Set seed for reproducibility
#np.random.seed(42)


################# MAIN METHOD ####################

clear_scene()
voronoi_test()

# Create icosphere and subdivide it to desired level
# NOTE: Do not apply higher level subdiv than 2, it crashes blender.
bpy.ops.mesh.primitive_ico_sphere_add(enter_editmode=False, align='WORLD', location=(0,0,0), scale=ICO_SCALE)
ico = bpy.context.active_object
ico.modifiers.new(name="Subdivision", type='SUBSURF')
ico.modifiers["Subdivision"].levels = SUBDIVISION_LEVELS
bpy.ops.object.modifier_apply({"object": ico}, modifier="Subdivision")

# Create inner and outer ico spheres
bpy.ops.object.duplicate(linked=False)
inner_ico = bpy.context.active_object
bpy.ops.object.duplicate(linked=False)
outer_ico = bpy.context.active_object
inner_ico.scale = tuple(x*INNER_SCALE_COEFF for x in ico.scale)
outer_ico.scale = tuple(x*OUTER_SCALE_COEFF for x in ico.scale)
ico.name = "Surface_Med"
inner_ico.name = "Surface_Inner"
outer_ico.name = "Surface_Outer"


# Generate seeds for Voronoi diagram
bpy.context.view_layer.objects.active = ico
bpy.ops.object.mode_set(mode='EDIT')
mesh = bmesh.from_edit_mesh(ico.data)
bmesh.ops.triangulate(mesh, faces=mesh.faces[:]) # Triangulate
face_count = len(mesh.faces)
#print(f"Mesh has {face_count} faces and {len(mesh.verts)} vertices.")

points = []
for face in mesh.faces:
    centroid = face.calc_center_median()
    points.append(centroid)
    for v in face.verts:
        v_loc = ico.matrix_world @ v.co
        p = 0.5*(centroid + v_loc)
        points.append(p)
#        q = 0.25*(centroid - v_loc) + centroid # Can be used in case the nuclei are not dense enough. - ck
#        points.append(q)
for vert in mesh.verts:
    v_loc = ico.matrix_world @ vert.co
    points.append(v_loc)
bpy.ops.object.mode_set(mode='OBJECT')
# Retain only points that lie between min and max CUT_BOX
# NOTE: Comment this line out if you want to get Voronoi cells on the whole surface.
# WARNING: This could lead to very long computation times.
points = [p for p in points if all(min_c <= val <= max_c for val, min_c, max_c in zip(p, MIN_CUT_BOX, MAX_CUT_BOX))]
assert len(points) <= NUCLEI_LIMIT, f"About to render {len(points)} nuclei. Stopped to avoid long compute time.\nYou can reduce the slice thickness to generate less nuclei." 
if RANDOM_TRANSLATE:
    points = [[p[i] + random.uniform(-1, 1)*MAX_TRANSLATE for i in range(3)] for p in points]
#add_point_cloud(points, radius = 0.01) # Render seed points


# Add auxiliary boundary points to ensure that the base Voronoi regions are bounded.
# The regions of the auxiliary points won't be.
auxiliary_points = get_lattice_points(MIN_COORDS, MAX_COORDS, PADDING)
#add_point_cloud(auxiliary_points, radius = 0.2)

# Generate the Voronoi diagram
vor = Voronoi(points + auxiliary_points)
#print_voronoi_stats(vor)
fr_points = finite_region_points(vor)
#print(f"Finite region points: {fr_points}")
ridges = compute_faces_by_seeds(vor, fr_points)

for idx, point_idx in enumerate(fr_points):
    add_region(vertices = vor.vertices,
                       faces = ridges[point_idx],
                       idx = point_idx)
    # Set the 3D cursor to the desired position
    bpy.context.scene.cursor.location = points[point_idx]
    # Set the object's origin to the 3D cursor location
    bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
    # Scale the object (adjust the scale factors as needed)
    scale = REGION_SCALE
    bpy.ops.transform.resize(value=(scale, scale, scale), orient_type='LOCAL')
    obj = bpy.context.active_object
    obj.select_set(False)

# Collect region objects in the scene
region_objects = get_objects_with("CellObject")
# Run the intersection function to get polytopes representing cell membranes
region_objects = intersect_with_object(region_objects, outer_ico)
region_objects = subtract_object(region_objects, inner_ico)
box = add_box(MIN_CUT_BOX, MAX_CUT_BOX)
region_objects = intersect_with_object(region_objects, box)

# Turn each polytope in a mesh representing its nucleus
nucleus_objects = add_nuclei_shaped(region_objects, NUCLEI_SCALE)
# Remove too large arrtifacts
artifacts = [] # NOTE: This lists too large regions
for obj in nucleus_objects:
    # Calculate bounding box
    diameter = np.max([max(b) - min(b) for b in zip(*obj.bound_box)])
    if diameter > min(ICO_SCALE): # NOTE: Maybe there is a better threshold. - ck
        artifacts.append(obj)
nucleus_objects = [obj for obj in nucleus_objects if obj not in artifacts]

# Intersect again
box.scale = tuple(s*a for s,a in zip(box.scale, (1,1,TISSUE_CUT_RATIO)))
nucleus_objects = intersect_with_object(nucleus_objects, box)
# Create inner and outer hull
inner_hull, outer_hull = intersect_with_object([inner_ico, outer_ico], box)
crypt_objects = nucleus_objects + [inner_hull, outer_hull]
# Transform crypt objects
rotate_objects(crypt_objects, Z_ROT_ANGLE)
translate_objects(crypt_objects, CENTER_LOC)
# Remove auxiliary objects
remove_objects(artifacts)
remove_objects([ico, box])
remove_objects(region_objects)