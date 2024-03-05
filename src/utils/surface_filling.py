import bpy
import bmesh
import numpy as np
import random

from itertools import combinations
from mathutils import Vector, Matrix

def edge_length(mesh, edge):
    v1 = mesh.vertices[edge.vertices[0]].co
    v2 = mesh.vertices[edge.vertices[1]].co
    return (v2 - v1).length

def min_max_edges(mesh):
    lengths = [edge_length(mesh, edge) for edge in mesh.edges]
    res = (min(lengths), max(lengths)) if lengths else (0, 0)
    return res

def triangulate_object(obj):
    me = obj.data
    # Get a BMesh representation
    bm = bmesh.new()
    bm.from_mesh(me)
    bmesh.ops.triangulate(bm, faces=bm.faces[:])
    # V2.79 : bmesh.ops.triangulate(bm, faces=bm.faces[:], quad_method=0, ngon_method=0)
    # Finish up, write the bmesh back to the mesh
    bm.to_mesh(me)
    bm.free()
    
def nearest_point_dist(pt, points):
    distances = [np.linalg.norm(pt, q) for q in points]
    return np.min(distances)

def add_point_cloud(locations, radius):
    # Create a small sphere object for each base point
    for idx, location in enumerate(locations):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=location)
        sphere = bpy.context.active_object
        sphere.name = f"Point_{idx}"
        
def mesh_area(mesh):
    # Initialize total area
    total_area = 0.0
    # Iterate over the faces and compute their areas
    for face in mesh.polygons:
        # Get the vertices of the face
        verts = [mesh.matrix_world @ mesh.vertices[i].co for i in face.vertices]
        assert len(verts) == 3, "Face must be triangle"
        a = verts[1] - verts[0]
        b = verts[2] - verts[0]
        # Compute the area of the face using cross product
        area = a.cross(b).length
        total_area += np.abs(area) / 2.0
    return total_area

def refine_mesh(mesh, delta, smoothness=1):
    bpy.ops.object.mode_set(mode='EDIT')
    for face in mesh.polygons:
        # Get longest side length of face
        vertices = [mesh.vertices[index].co for index in face.vertices]
        diameter = max((v1 - v2).length for v1 in vertices for v2 in vertices)

        # If face is too large subdivide it
        if diameter > delta:
            face.select = True
            level = int(np.ceil(diameter / delta))
            bpy.ops.mesh.subdivide(number_cuts=level, smoothness=smoothness)
            bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')
    return mesh

def sample_centers(verts, dist, max_count):
    sampled_verts = []
    for v in verts:
        is_free = True
        for q in sampled_verts:
            if (v.co-q.co).length < dist:
                is_free = False
                break
        if is_free:
            sampled_verts.append(v)
        if len(sampled_verts) == max_count:
            break
    return sampled_verts
        
def sample_fillers(mesh_verts, center_verts, center_dist, fill_dist, max_count):
    filler_verts = []
    for v in mesh_verts:
        is_free = True
        for q in center_verts:
            if (v.co-q.co).length < 0.5*(center_dist+fill_dist):
                is_free = False
                break
        for q in filler_verts:
            if (v.co-q.co).length < fill_dist:
                is_free = False
                break
        if is_free:
            filler_verts.append(v)
        if len(filler_verts) == max_count:
            break
    return filler_verts
        
        
def add_oriented_points(points, directions, scale, radius):
    for idx, (pt, dir) in enumerate(zip(points, directions)):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=radius)
        sphere = bpy.context.active_object
        sphere.location = pt
        
        rotation_matrix = Matrix.Translation(sphere.location) @ dir.to_track_quat('X').to_matrix().to_4x4()
        # Apply the rotation to the object
        sphere.matrix_world = rotation_matrix
        sphere.scale = scale
        
        sphere.name = f"Point_{idx}"
        

def fill_surface(obj, max_point_count, attribute, filler_scale):
    '''
    First refines the mesh of the object until the edges are sufficiently small.
    Then samples the maximal number of vertices on the refined mesh such that intersection free placement of nuclei is possible.
    '''
    min_dist = 2 * attribute.size * attribute.scale[1] # NOTE: Use medium radius of scale, as max radius is for normal direction. - ck
    # NOTE: That's the best idea so far for creating packed surfaces. Is there better way? - ck
    fill_dist = min_dist * filler_scale # Smaller radius of nucleus to fill the gaps between large ones
    # NOTE: Find optimal value, min_dist/3 looks good? - ck
    mesh_delta = min_dist/3 # Refine mesh until vertices in mesh grid are at most mesh_delta apart.

    triangulate_object(obj)
    mesh = obj.data

    # Refine mesh
    # _, max_edge = min_max_edges(mesh)
    # print(f"Before refinement: {len(mesh.vertices)} vertices, mesh offset {max_edge}")
    mesh = refine_mesh(mesh, mesh_delta)
    # _, max_edge = min_max_edges(mesh)
    # print(f"After refinement: {len(mesh.vertices)} vertices, mesh offset {max_edge}")

    # Sample nuclei centers
    grid_verts = list(mesh.vertices)
    random.shuffle(grid_verts)
    center_verts = sample_centers(grid_verts, min_dist, max_point_count)
    filler_verts = sample_fillers(grid_verts, center_verts, min_dist, fill_dist, max_point_count)
    #print(f"Placed {len(center_verts)} centers and {len(filler_verts)} fillers")
    return center_verts, filler_verts, mesh_delta
