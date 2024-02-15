import bmesh
import bpy
import random

from scipy.spatial import Voronoi # NOTE: You might install it directly into Blender's python using pip. - ck

from src.utils.helper_methods import *

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

def compute_mean_scale(object):
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

def remove_objects_inside_mesh(objects, empty_regions):
    nuclei_to_remove = []
    nuclei_to_keep = []
    for nucleus_object in objects:
        for region in empty_regions:
            if is_point_inside_mesh(region, nucleus_object.location):
                #print(f"Removing nucleus {nucleus_object.name}")
                nuclei_to_remove.append(nucleus_object)
            else:
                nuclei_to_keep.append(nucleus_object)
    remove_objects(nuclei_to_remove)
    return nuclei_to_keep

def sample_points_on_mesh(obj, nuclei_limit, cut_box, max_translate=None):
    # Generate seeds for Voronoi diagram
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.mode_set(mode='EDIT')
    mesh = bmesh.from_edit_mesh(obj.data)
    bmesh.ops.triangulate(mesh, faces=mesh.faces[:]) # Triangulate # NOTE: Check if this is truly necessary. - ck
    #print(f"Mesh has {len(mesh.faces)} faces and {len(mesh.verts)} vertices.")
    points = []
    for face in mesh.faces:
        centroid = face.calc_center_median()
        points.append(centroid)
        for v in face.verts:
            v_loc = obj.matrix_world @ v.co
            p = 0.5*(centroid + v_loc)
            points.append(p)
    #        q = 0.25*(centroid - v_loc) + centroid # NOTE: Can be used in case the nuclei are not dense enough. - ck
    #        points.append(q)
    for vert in mesh.verts:
        v_loc = obj.matrix_world @ vert.co
        points.append(v_loc)
    bpy.ops.object.mode_set(mode='OBJECT')
    # Retain only points that lie between min and max CUT_BOX
    min_cut_box, max_cut_box = cut_box
    points = [p for p in points if all(min_c <= val <= max_c for val, min_c, max_c in zip(p, min_cut_box, max_cut_box))]
    if len(points) > nuclei_limit:
        print(f"WARNING: Adding only {nuclei_limit} nuclei instead {len(points)} nuclei due to long compute time.\nYou can reduce the slice thickness to generate less nuclei.")
        points = random.sample(points, k=nuclei_limit)
    if max_translate is not None:
        points = [[p[i] + random.uniform(-1, 1)*max_translate for i in range(3)] for p in points]
    return points
    