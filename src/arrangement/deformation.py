import numpy as np
import time
import random

from mathutils import Vector
import mathutils.noise as noise_blender
import bpy

from src.utils.geometry import subdivide

# Permutation table for Perlin noise
permutation = [151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225,
               140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247,
               120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177,
               33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165,
               71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211,
               133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25,
               63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196,
               135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217,
               226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206,
               59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119, 248,
               152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, 129, 22,
               39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246,
               97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51,
               145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84,
               204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114,
               67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180]

# Duplicate the permutation table to avoid overflow
permutation += permutation

def fade(t):
    # Fade function to smooth the interpolation
    return t * t * t * (t * (t * 6 - 15) + 10)

def lerp(t, a, b):
    # Linear interpolation
    return a + t * (b - a)

def grad(hash, x, y, z):
    # Gradient function to calculate the dot product
    h = hash & 15
    u = x if h < 8 else y
    v = y if h < 4 else (z if h == 12 or h == 14 else x)
    return ((u if (h & 1) == 0 else -u) +
            (v if (h & 2) == 0 else -v))

def noise(x, y, z):
    # Calculate the unit cube that the point falls in
    X = int(np.floor(x)) & 255
    Y = int(np.floor(y)) & 255
    Z = int(np.floor(z)) & 255

    # Find the relative x, y, z of the point in the cube
    x -= np.floor(x)
    y -= np.floor(y)
    z -= np.floor(z)

    # Compute fade curves for each of x, y, z
    u = fade(x)
    v = fade(y)
    w = fade(z)

    # Hash coordinates of the 8 cube corners
    A = permutation[X] + Y
    AA = permutation[A] + Z
    AB = permutation[A + 1] + Z
    B = permutation[X + 1] + Y
    BA = permutation[B] + Z
    BB = permutation[B + 1] + Z

    # Add blended results from 8 corners of the cube
    return lerp(w, lerp(v, lerp(u, grad(permutation[AA], x, y, z),
                                   grad(permutation[BA], x-1, y, z)),
                           lerp(u, grad(permutation[AB], x, y-1, z),
                                   grad(permutation[BB], x-1, y-1, z))),
                   lerp(v, lerp(u, grad(permutation[AA+1], x, y, z-1),
                                   grad(permutation[BA+1], x-1, y, z-1)),
                           lerp(u, grad(permutation[AB+1], x, y-1, z-1),
                                   grad(permutation[BB+1], x-1, y-1, z-1))))


def perlin(x, y, z, amplitude, frequency, octave_count, persistence, lacunarity):
    value = 0.0
    for _ in range(octave_count):
        value += amplitude * noise(x * frequency, y * frequency, z * frequency)
        amplitude *= persistence
        frequency *= lacunarity
    return value

def deformation(input_vector, params):
    amplitude, frequency, octave_count, persistence, lacunarity = params
    dx = perlin(input_vector[0], input_vector[1], input_vector[2], amplitude, frequency, octave_count, persistence, lacunarity)
    dy = perlin(input_vector[0] + 100, input_vector[1] + 100, input_vector[2] + 100, amplitude, frequency, octave_count, persistence, lacunarity)
    dz = perlin(input_vector[0] + 200, input_vector[1] + 200, input_vector[2] + 200, amplitude, frequency, octave_count, persistence, lacunarity)

    # Create deformation vector
    # NOTE: We normalize the deformation vector to a unit vector. Otherwise the deformations could be to drastic. 
    deformation_vector = Vector((dx, dy, dz)).normalized() 
    return deformation_vector 

def deform_objects(obj_list, deformation_strength=0.25): # TODO add tuple / properties and strength
    """
    Deforms a list of objects using Perlin noise parameters.
    Parameters:
        obj_list (list): List of objects to be deformed.
        deformation_strength (float, optional): Strength of the deformation. Default is 0.25.
        Notes:
            - The amplitude is currently set to 1, and the final deformation strength should not exceed 1/4 of the object radius.
            - There is a potential issue where deformation of touching large and small cells could lead to intersection. This needs to be checked.
            Perlin Noise Parameters:
            - amplitude: Max strength of deformation.
            - wavelength: Average distance between two maximal deformations in the same direction.
            - octave_count: Number of octaves in the noise function.
            - persistence: Persistence of the noise function.
            - lacunarity: Lacunarity of the noise function.
        The function calculates the centroid of each object and deforms each vertex based on its distance from the centroid and the Perlin noise parameters. The objects are then subdivided once.
        Timing:
        The function prints the time taken to deform all objects.
        TODO:
            - Add tuple/properties and strength.
            - Verify if the current approach to setting amplitude and deformation strength is appropriate.
            - Check for potential intersections when deforming touching large and small cells.
    """
    # Set Perlin noise parameters
    # NOTE: Currently amplitude is set to 1 and the final deformation strength should not exceed 1/4 of the object radius.
    # TODO: Is this a good idea? Check this. -ck
    # It could happen that a deformation of touching large cell and small cell leads to intersection. Check this!
    amplitude = 1 # Max strength of deformation
    wavelength = 0.05 # Average distance between two maximal deformations in same direction
    octave_count = 2
    persistence = 1.0
    lacunarity = 2.0
    parameters = (amplitude, 1 / wavelength, octave_count, persistence, lacunarity)

    start = time.time()
    for obj in obj_list:
        centroid = Vector(np.mean([obj.matrix_world @ v.co for v in obj.data.vertices], axis=0))
        for v in obj.data.vertices:
            rad = (obj.matrix_world @ v.co - centroid).length
            v.co = v.co + deformation(v.co, parameters) * rad * deformation_strength
        subdivide(obj, 1)
    print(f"Deformed {len(obj_list)} objects in {time.time() - start} seconds")


# def deform_objects(obj_dict):
#     '''
#     Deforms a dictionary of objects using Perlin noise parameters.
#     Parameters:
#         obj_dict (dict): Nested dictionary of objects to be deformed including the noise parameters.
#         {object_name: {object: object, deformations: {name: {strength: float, params: **kwargs}}}
#     '''
#     for obj_name, d in obj_dict.items():
#         obj = d['object']
#         deformations = d['deformations']
#         for name, deformation in deformations.items():
#             deformation_strength = deformation['strength']
#             deformation_params = deformation['params']
#             if name == 'perlin':
#                 deform_objects_perlin([obj], deformation_strength, deformation_params)


def elastic_deform(obj_name, deformation_strength=0.025, noise_scale=3, seed=3, quad_falloff=0):
    """
    Apply elastic deformation to a mesh object in Blender.
    This function deforms the vertices of a specified mesh object using a noise-based algorithm.
    The deformation is controlled by the strength, scale, and seed parameters.
    Parameters:
        obj_name (str): The name of the object to deform.
        deformation_strength (float, optional): The strength of the deformation. Default is 0.025.
        noise_scale (float, optional): The scale of the noise used for deformation. Default is 3.
        seed (int, optional): The seed for the random number generator to ensure reproducibility. Default is 3.
    Returns:
        None
    """
    # ensure that the object exists
    obj = bpy.data.objects.get(obj_name)
    if obj is None:
        print(f"Object '{obj_name}' not found!")
        return
    
    # ensure the object is a mesh
    if obj.type != 'MESH':
        print(f"Object '{obj_name}' is not a mesh!")
        return
    
    # get mesh
    bpy.ops.object.mode_set(mode='OBJECT')
    mesh = obj.data
    
    # apply random deformation to vertices
    random.seed(seed)
    seed_offset_x = Vector((random.randint(0, 1000), random.randint(0, 1000), random.randint(0, 1000)))
    seed_offset_y = Vector((random.randint(0, 1000), random.randint(0, 1000), random.randint(0, 1000)))
    for vertex in mesh.vertices:
        world_space_coord = obj.matrix_world @ vertex.co
        shift_x = noise_blender.noise(world_space_coord * noise_scale + seed_offset_x) * deformation_strength
        shift_y = noise_blender.noise(world_space_coord * noise_scale + seed_offset_y) * deformation_strength
        shift_x = quad_falloff * (shift_x * shift_x) + (1-quad_falloff) * shift_x
        shift_y = quad_falloff * (shift_y * shift_y) + (1-quad_falloff) * shift_y
        world_space_coord += Vector((shift_x, shift_y, 0))
        vertex.co = obj.matrix_world.inverted() @ world_space_coord
    mesh.update()