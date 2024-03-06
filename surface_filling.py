import bpy
import bmesh
import numpy as np
import random

from itertools import combinations
from mathutils import Vector, Matrix


def delete_objects():
    # Select all objects 
    for obj in bpy.context.scene.objects:
        obj.hide_viewport = False
    bpy.ops.object.select_all(action='SELECT')
    # Delete the selected objects
    bpy.ops.object.delete()
    # Remove orphaned meshes
    for mesh in bpy.data.meshes:
        bpy.data.meshes.remove(mesh, do_unlink=True)

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
        verts = [obj.matrix_world @ mesh.vertices[i].co for i in face.vertices]
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
        vertices = [obj.data.vertices[index].co for index in face.vertices]
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
        
    

#############  MAIN

# IDEA: Given a mesh, triangulate it, subdivide the faces until all edges are below a threshold, and then collapse too small edges

#bpy.ops.object.mode_set(mode='OBJECT')
delete_objects()

MIN_DIST = 0.12
FILL_DIST = 0.09
SCALE = (2,1,1)

MAX_POINT_COUNT = 200

# TODO: Find optimal value, maybe MIN_DIST/3?
mesh_delta = MIN_DIST/3 # Maximal edge length in refined mesh


bpy.ops.mesh.primitive_ico_sphere_add()
obj = bpy.context.active_object
triangulate_object(obj)
obj.scale = (1, 0.5, 0.5)
bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)
mesh = obj.data

area = mesh_area(mesh)
radius = MIN_DIST/2
print(f"\nArea of Mesh: {area}")
print(f"Number of faces: {len(mesh.polygons)}")
print(f"Max number of cells of radius {radius} for dense packing: {0.6*area/(np.pi*radius*radius)}")

# Refine mesh
min_edge, max_edge = min_max_edges(mesh)
print(f"Before refinement: {len(mesh.vertices)} vertices, mesh offset {max_edge}")
mesh = refine_mesh(mesh, mesh_delta)
min_edge, max_edge = min_max_edges(mesh)
print(f"After refinement: {len(mesh.vertices)} vertices, mesh offset {max_edge}")

# Sample nuclei centers
grid_points = [v.co for v in mesh.vertices]
grid_normals = [v.normal for v in mesh.vertices]

grid_verts = list(mesh.vertices)
random.shuffle(grid_verts)
center_verts = sample_centers(grid_verts, MIN_DIST, MAX_POINT_COUNT)
filler_verts = sample_fillers(grid_verts, center_verts, MIN_DIST, FILL_DIST, MAX_POINT_COUNT)


print(f"Placed {len(center_verts)} centers and {len(filler_verts)} fillers")

center_points = [v.co for v in center_verts]
center_normals = [v.normal for v in center_verts]
add_oriented_points(center_points, center_normals, SCALE, MIN_DIST/2)

filler_points = [v.co for v in filler_verts]
filler_normals = [v.normal for v in filler_verts]
add_oriented_points(filler_points, filler_normals, SCALE, FILL_DIST/2)


#add_point_cloud(center_points, MIN_DIST/2)
#add_point_cloud(filler_points, FILL_DIST/2)
##mesh.hide_viewport = True