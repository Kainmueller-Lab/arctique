import bmesh
import bpy
import numpy as np
from mathutils import Matrix, Vector


def move_selection(offset_vector):
    selection = bpy.context.selected_objects
    for obj in selection:
        obj.location += offset_vector
        
        
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


def add_point_cloud(locations, radius):
    # Create a small sphere object for each base point
    for idx, location in enumerate(locations):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=location)
        sphere = bpy.context.active_object
        sphere.name = f"Point_{idx}"


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

def get_lattice_points(min_coords, max_coords):
    points = []
    cube_center = Vector([0.5 * (i + j) for i, j in zip(min_coords, max_coords)])
    diameter = Vector([j - i for i, j in zip(min_coords, max_coords)])
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

def intersect_with_object(target_objects, box_object):
    # Iterate through each target object
    for target_object in target_objects:
        boolean = target_object.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
        boolean.operation = 'INTERSECT'
        boolean.object = box_object

        # Apply the boolean modifier and remove the cube object
        bpy.ops.object.modifier_apply({"object": target_object}, modifier="Boolean Modifier")
    bpy.data.objects.remove(box_object, do_unlink=True)
    return target_objects

def add_nuclei_from(cell_objects, nuclei_scale):
    nucleus_objects = []
    for cell_object in cell_objects:
        bpy.ops.mesh.primitive_cube_add(enter_editmode=False, align='WORLD', location=cell_object.location, scale=(1, 1, 1))
        nucleus_object = bpy.context.active_object
        index = cell_object.name.split('_')[1]
        nucleus_object.name = f"NucleusObject_{index}"
        shrinkwrap = nucleus_object.modifiers.new(name="Shrinkwrap Modifier", type='SHRINKWRAP')
        shrinkwrap.target = cell_object
        bpy.ops.object.modifier_apply(modifier="Shrinkwrap Modifier")
        subsurf = nucleus_object.modifiers.new("Subsurface Modifier", type='SUBSURF')
        subsurf.levels = 2
        bpy.ops.object.modifier_apply(modifier="Subsurface Modifier")
        nucleus_object.scale = (nuclei_scale,) * 3
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



def lerp(a, b, t):
    return a*(1-t) + b*t
    

# Function to set the local x-axis orientation of an object along a direction vector
def set_orientation(obj, direction_vector):
    # Calculate the rotation matrix to align the object  with the direction vector
    rotation_matrix = Matrix.Translation(obj.location) @ direction_vector.to_track_quat('X').to_matrix().to_4x4()
    # Apply the rotation to the object
    obj.matrix_world = rotation_matrix


def compute_normal(curve_function, t): # Only in x-y-plane
    epsilon = 1e-4  # Small value to avoid division by zero
    # Compute the tangent vector
    tangent_vector = (curve_function(t + epsilon) - curve_function(t - epsilon))
    # Normalize the tangent vector
    tangent_vector /= np.linalg.norm(tangent_vector)
    assert tangent_vector.z == 0

    # Compute the normal vector (perpendicular to the tangent)
    normal_vector = Vector([tangent_vector.y, -tangent_vector.x, 0])
    return normal_vector

def random_unit_vector():
    # Generate random values for x, y, and z
    x, y, z = np.random.uniform(-1, 1, 3)

    # Create a vector
    vector = np.array([x, y, z])

    # Normalize the vector to make it a unit vector
    unit_vector = vector / np.linalg.norm(vector)

    return unit_vector
