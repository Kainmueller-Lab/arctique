import numpy as np
import time

from mathutils import Vector

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

def deform_objects(obj_list, deformation_strength=0.25):
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