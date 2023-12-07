import bpy
import random
from mathutils import Vector

# IMPORT SOURCES
import sys
import os

dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir )

import src.helper_methods as hm

# this next part forces a reload in case you edit the source after you first start the blender session
import imp
imp.reload(hm)


class RandomDeformedCell:
    def __init__(self, size=1.0, scale=(1,1,1), deformation_strength=0.1, orientation=None):
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.orientation = orientation

        # Create a cube
        bpy.ops.mesh.primitive_cube_add(size=self.size, location=(0, 0, 0))
        self.cell_object = bpy.context.active_object
        
        # Apply two levels of subdivision surface
        modifier = self.cell_object.modifiers.new("Subsurface Modifier", "SUBSURF")
        modifier.levels = 2
        
        # Orient
        if orientation:
            hm.set_orientation(self.cell_object, orientation)

        # Scale 
        self.cell_object.scale = self.scale
        
        # Deform the vertices randomly
        self.deform_vertices()
        

    def deform_vertices(self):
        # Iterate through each vertex and deform its position
        for vertex in self.cell_object.data.vertices:
            original_position = vertex.co.copy()
            deformation_vector = Vector([
                random.uniform(-self.deformation_strength, self.deformation_strength),
                random.uniform(-self.deformation_strength, self.deformation_strength),
                random.uniform(-self.deformation_strength, self.deformation_strength)
            ])*Vector(self.scale)
            vertex.co = original_position + deformation_vector
            
            
class RandomCellField:
    def __init__(self, num_cells, min_coords, max_coords):
        self.num_cells = num_cells
        self.min_coords = min_coords
        self.max_coords = max_coords
        
    def generate_cells(self, cell_size, cell_scale, deformation_strength):
        for _ in range(self.num_cells):
            cell = RandomDeformedCell(cell_size, cell_scale, deformation_strength)
            cell.cell_object.location = Vector([
                random.uniform(self.min_coords.x, self.max_coords.x),
                random.uniform(self.min_coords.y, self.max_coords.y),
                random.uniform(self.min_coords.z, self.max_coords.z)
            ])
            
            
class CellCurve:
    def __init__(self, num_cells, curve_function, curve_noise, curve_interval=(0,1), scale_weight=None):
        self.num_cells = num_cells
        self.curve_function = curve_function
        self.curve_noise = curve_noise
        self.curve_interval = curve_interval
        self.scale_weight = scale_weight
        self.increment = (curve_interval[1]-curve_interval[0]) / num_cells
        
    def generate_cells(self, cell_size, cell_scale, deformation_strength):
        t = self.curve_interval[0]
        for idx in range(self.num_cells):
            # Change orientation of cells such that x-axis points to curve normal
            normal = hm.compute_normal(self.curve_function, t)
            # Scale cells along the curve according to scale_weight
            current_scale = tuple([self.scale_weight(t) * s for s in cell_scale]) if self.scale_weight!=None else cell_scale
            cell = RandomDeformedCell(cell_size, current_scale, deformation_strength, normal)
            cell.cell_object.location = Vector([
                self.curve_function(t).x + random.uniform(-self.curve_noise.x, self.curve_noise.x),
                self.curve_function(t).y + random.uniform(-self.curve_noise.y, self.curve_noise.y),
                # TODO: Change this once we do not delete cells with zloc > 0
                self.curve_function(t).z + random.uniform(0, self.curve_noise.z) 
            ])
            # TODO: Handle increment weight better such that it can be tweaked by parameters
            increment_weight = self.scale_weight(t)*self.scale_weight(t) if self.scale_weight else 1
            current_increment = self.increment * increment_weight
            t += current_increment

        