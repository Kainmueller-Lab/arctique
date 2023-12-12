import bpy 
import random

from typing import List, Optional

from mathutils import Vector

from rHEnder.shading.shading import Material
from rHEnder.utils.helper_methods import set_orientation


class CellAttribute():
    def __init__(self):
        self.cell_type = None
        self.size = None
        self.scale = None
        self.deformation_strength = None
        self.attribute_name = None


class CellAttributeA(CellAttribute):
    def __init__(self, cell_type = "A", size = 0.04, scale = (1,1,1), deformation_strength = 0.007, attribute_name = "Cell Type A"):
        self.cell_type = cell_type
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name

class CellAttributeB(CellAttribute):
    def __init__(self, cell_type = "B", size = 0.07, scale = (2.3,1,1), deformation_strength = 0.007, attribute_name = "Cell Type B"):
        self.cell_type = cell_type
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name

class Cell:
    def __init__(self, cell_id: int, cell_attributes: CellAttribute, location: Vector, cell_name: str, orientation: Optional[Vector] = None):
        self.cell_id = cell_id
        self.cell_attributes = cell_attributes
        self.location = location
        self.orientation = orientation
        self.cell_name = cell_name
        self.semantic_id = None
        self.material = None
        self.cell_object = None

    def set_material(self, material: Material):
        self.material = material

    def set_material(self, semantic_id: int):
        self.semantic_id = semantic_id

    def add(self):
        # Create a cube
        bpy.ops.mesh.primitive_cube_add(size=self.cell_attributes.size, location=self.location)
        self.cell_object = bpy.context.active_object
        
        # Apply two levels of subdivision surface
        modifier = self.cell_object.modifiers.new("Subsurface Modifier", "SUBSURF")
        modifier.levels = 2
        
        # Orient
        if self.orientation:
            set_orientation(self.cell_object, self.orientation)

        # Scale 
        self.cell_object.scale = self.cell_attributes.scale
        
        # Deform the vertices randomly
        self.deform_mesh(self.cell_attributes.deformation_strength)

    def deform_mesh(self, deformation_strength):
        # Iterate through each vertex and deform its position
        for vertex in self.cell_object.data.vertices:
            original_position = vertex.co.copy()
            deformation_vector = Vector([
                random.uniform(-1, 1),
                random.uniform(-1, 1),
                random.uniform(-1, 1)
            ])*Vector(self.cell_attributes.scale)*self.cell_attributes.deformation_strength
            vertex.co = original_position + deformation_vector