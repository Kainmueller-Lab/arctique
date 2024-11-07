import bmesh
import bpy
import numpy as np
import random

from mathutils import Vector
from scipy.spatial import Voronoi # NOTE: You might install it directly into Blender's python using pip. - ck

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

def finite_region_points(vor):
    fr_points = []
    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if -1 not in region and len(region) > 0:
            fr_points.append(point_idx)
    return fr_points

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
    return ridges

def voronoi_test():
    # Generate 27 integer points in the 0, -1, 1 lattice
    points = [[i, j, k] for i in range(-1, 2) for j in range(-1, 2) for k in range(-1, 2)]
    vor = Voronoi(points)
    fr_points = finite_region_points(vor)
    assert len(fr_points) == 1
    assert fr_points[0] == 13

def generate_voronoi_regions(all_points, auxiliary_points, region_scale):
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
        scale = region_scale
        bpy.ops.transform.resize(value=(scale, scale, scale), orient_type='LOCAL')
        obj = bpy.context.active_object
        obj.select_set(False)

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
    bm.verts.ensure_lookup_table() # NOTE: Necessary to allow vertex indexing
    # Add faces
    for face_indices in faces:
        bm.faces.new([bm.verts[i] for i in face_indices])
    # Update the mesh with the new geometry
    bm.to_mesh(mesh)
    bm.free()
    return obj

def get_cube_points(min_coords, max_coords, padding = 0):
    points = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                point = [max_coords[0]+padding if i else min_coords[0]-padding,
                        max_coords[1]+padding if j else min_coords[1]-padding,
                        max_coords[2]+padding if k else min_coords[2]-padding]
                points.append(point)
    return points

def get_octogon_points(min_coords, max_coords, padding = 0):
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
    diameter = Vector([j - i + 2*padding for i, j in zip(min_coords, max_coords)])
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for k in [-1, 0, 1]:
                if not (i==0 and j==0 and k==0):
                    point = cube_center + Vector([a*b for (a,b) in zip(diameter, [i,j,k])])
                    points.append(point)
    return points

def generate_lattice_parameters(theta, only_one=False):
   '''
   Theta defines a lattice of 4 isosceles parallelograms of area 1.
   Each such parallelogram defines a rotated ellipsoid fitting into it.
   '''
   theta = np.radians(theta)
   length = 1 / np.sin(theta)
   mu = 0.6 # NOTE: This is handpicked and steers the distance between crypts in the lattice. - ck
   x_axis = mu * length * np.cos(theta/2)
   y_axis = mu * length * np.sin(theta/2)
   v = Vector([length,0,0])
   w = Vector([length*np.cos(theta), length*np.sin(theta), 0])
   max_delta_scale = 0.05
   max_delta_angle = 10
   ico_scales = [(x_axis + max_delta_scale*random.uniform(-1,1), y_axis + max_delta_scale*random.uniform(-1,1)) for _ in range(4)]
   angles = [np.degrees(theta)/2 + max_delta_angle*random.uniform(-1,1) for _ in range(4)]
   base_pt = Vector([-0.5, -0.5, 0.5])
   lattice = [Vector([0,0,0]), v, w, v+w]
   centers = [base_pt + vec for vec in lattice]
   if only_one:
      ico_scales = [ico_scales[0]]
      angles = [angles[0]]
      centers = [centers[0]]
   return ico_scales, angles, centers