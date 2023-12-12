import bpy 
import random

from typing import List, Optional

from mathutils import Vector

from src.shading.shading import Material
from src.utils.helper_methods import set_orientation

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
    def __init__(self, idx: int, location: Vector, arrangement_id: int, arrangement_type: str, cell_attributes: CellAttribute, orientation: Optional[Vector] = None, index: int = 1):
        self.cell_attributes = cell_attributes
        self.cell_id = idx # TODO: Add arrangement id as id. - ck
        self.location = location
        self.orientation = orientation
        # TODO: Improve cell naming for better overview in blender. - ck
        self.cell_name = f"Cell_{idx}_Arr_{arrangement_type}_{arrangement_id}_Type_{self.cell_attributes.cell_type}"
        self.arrangement_type = None
        self.index = index
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
        self.cell_object.name = self.cell_name
        
        pass_index = bpy.props.IntProperty(name="Pass Index", subtype='UNSIGNED')
        bpy.context.object.pass_index = self.index
        
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