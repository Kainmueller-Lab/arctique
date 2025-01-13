import bpy
import bmesh
import numpy as np
import random

from itertools import combinations
from mathutils import Vector, Matrix

def triangulate_object(obj):
    """
    Triangulates a mesh object. This is necessary for the mesh to be renderable.

    Parameters:
        obj (bpy.types.Object): The object to triangulate.
    """
    me = obj.data
    # Get a BMesh representation
    bm = bmesh.new()
    bm.from_mesh(me)
    bmesh.ops.triangulate(bm, faces=bm.faces[:])
    # V2.79 : bmesh.ops.triangulate(bm, faces=bm.faces[:], quad_method=0, ngon_method=0)
    # Finish up, write the bmesh back to the mesh
    bm.to_mesh(me)
    bm.free()

def refine_mesh(mesh, delta, smoothness=1):
    """
    Refines a mesh by subdividing its faces based on a given maximum diameter.

    Parameters:
        mesh (bpy.types.Mesh): The mesh object to refine.
        delta (float): The maximum allowable diameter for any face. Faces larger than this will be subdivided.
        smoothness (float, optional): The smoothness factor for subdivision. Default is 1.

    Returns:
        bpy.types.Mesh: The refined mesh object with subdivided faces.
    """
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
    """
    Samples a subset of vertices from a list, ensuring a minimum distance between any two sampled vertices.

    Parameters:
        verts (list): A list of vertex objects with a 'co' attribute representing their coordinates.
        dist (float): The minimum required distance between any two sampled vertices.
        max_count (int): The maximum number of vertices to sample.

    Returns:
        list: The list of sampled vertex objects.
    """

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
    """
    Samples a subset of vertices to serve as fillers, ensuring a minimum distance from both center vertices and other fillers.

    Parameters:
        mesh_verts (list): A list of vertex objects with a 'co' attribute representing their coordinates.
        center_verts (list): A list of center vertex objects with a 'co' attribute representing their coordinates.
        center_dist (float): The distance used to calculate the separation threshold from center vertices.
        fill_dist (float): The minimum required distance between any two filler vertices.
        max_count (int): The maximum number of filler vertices to sample.

    Returns:
        list: The list of sampled filler vertex objects.
    """

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
        

def fill_surface(obj, max_point_count, attribute, filler_scale):
    """
    Fills the surface of a mesh with cell nuclei.

    Parameters:
        obj (bpy.types.Object): The object to fill with cell nuclei.
        max_point_count (int): The maximum number of cell nuclei to place on the surface.
        attribute (CellAttribute): The cell attribute to use for filling the surface.
        filler_scale (float): The ratio of the filler nucleus size to the center nucleus size.

    Returns:
        tuple: A tuple containing:
            - center_verts (list): The list of sampled center vertex objects.
            - filler_verts (list): The list of sampled filler vertex objects.
            - mesh_delta (float): The distance used to refine the mesh.
    """
    min_dist = 2 * attribute.radius * attribute.scale[1] # NOTE: Use medium radius of scale, as max radius is for normal direction. - ck
    # NOTE: That's the best idea so far for creating packed surfaces. Is there better way? - ck
    fill_dist = min_dist * filler_scale # Smaller radius of nucleus to fill the gaps between large ones
    # NOTE: Find optimal value, min_dist/3 looks good? - ck
    mesh_delta = min_dist/3 # Refine mesh until vertices in mesh grid are at most mesh_delta apart.

    # Refine mesh
    triangulate_object(obj)
    mesh = obj.data
    mesh = refine_mesh(mesh, mesh_delta)

    # Sample nuclei centers
    grid_verts = list(mesh.vertices)
    random.shuffle(grid_verts)
    center_verts = sample_centers(grid_verts, min_dist, max_point_count)
    filler_verts = sample_fillers(grid_verts, center_verts, min_dist, fill_dist, max_point_count)
    return center_verts, filler_verts, mesh_delta
