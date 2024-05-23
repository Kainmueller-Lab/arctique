import bpy 
import bmesh
import random
import math

from enum import Enum
from mathutils import Vector
from typing import Optional

from src.shading.shading import Material
from src.utils.geometry import set_orientation

class CellType(Enum):
    PLA = 1
    LYM = 2
    EOS = 2
    FIB = 3
    EPI = 4
    GOB = 5


class CellAttribute():
    def __init__(self):
        self.cell_type = None
        self.size = None
        self.scale = None
        self.deformation_strength = None
        self.attribute_name = None
        self.max_bending_strength = None

# NOTE: The size and scale convention for CellAttributes is as follows:
# - size: is the maximal ellipsoid radius of the nucleus 
# - scale: is the ellipsoid scale of the nucleus
# The scale of an attribute should always be normalized such the maximal scale value is 1. 
# This ensures that the maximal diameter of a cell nucleus is twice its size.

class PLA(CellAttribute):
    def __init__(self, 
                 size = 0.06,
                 scale = (1,1,1),
                 deformation_strength = 0.8,
                 attribute_name = "Plasma Cell",
                 max_bending_strength = 0.2):
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class LYM(CellAttribute):
    def __init__(self, 
                 size = 0.03,
                 scale = (1,1,1),
                 deformation_strength = 0.8,
                 attribute_name = "Lymphocyte",
                 max_bending_strength = 0.2):
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class EOS(CellAttribute):
    def __init__(self, 
                 size = 0.06,
                 scale = (1,0.3,0.3),
                 deformation_strength = 0.8,
                 attribute_name = "Eosinophile",
                 max_bending_strength = 0.2):
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class FIB(CellAttribute):
    def __init__(self, 
                 size = 0.08,
                 scale = (1,0.2,0.2),
                 deformation_strength = 0.8,
                 attribute_name = "Fibroblast",
                 max_bending_strength = 0.2):
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class EPI(CellAttribute):
    def __init__(self, 
                 size = 0.03,
                 scale = (1,1,1),
                 deformation_strength = 0.8,
                 attribute_name = "Epithelial Cell",
                 max_bending_strength = 0.2):
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class GOB(CellAttribute):
    def __init__(self,
                 size = 0.03,
                 scale = (1,1,1),
                 deformation_strength = 0.8,
                 attribute_name = "Goblet Cell",
                 max_bending_strength = 0.2):
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength



class CellAttributeA(CellAttribute):
    def __init__(self, cell_type = "A",
                 size = 0.03,
                 scale = (1,1,1),
                 deformation_strength = 0.8,
                 attribute_name = "Cell Type A",
                 max_bending_strength = 0.2):
        self.cell_type = cell_type
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class CellAttributeB(CellAttribute):
    def __init__(self, cell_type = "B", size = 0.08, scale = (1, 0.5, 0.4), deformation_strength = 0.6, attribute_name = "Cell Type B", max_bending_strength = 0.3):
        self.cell_type = cell_type
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class CellAttributeC(CellAttribute):
    def __init__(self, cell_type = "C", size = 0.04, scale = (1, 1, 0.5), deformation_strength = 0.8, attribute_name = "Cell Type C", max_bending_strength = 0.3):
        self.cell_type = cell_type
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class CellAttributeEpi(CellAttribute):
    # TODO: Make cell type a fixed value
    def __init__(self, cell_type = "Epi", size = 0.1, scale = (1, 0.6, 0.6), deformation_strength = 0.2, attribute_name = "Epithelial cell", max_bending_strength = 0.3):
        self.cell_type = cell_type
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class MixAttribute:
    def __init__(self, true_attribute: CellAttribute, mixing_attribute: CellAttribute, mix: float):
        '''
        Produces a cell attribute that is a combination of two cell attributes in terms of numerical values.
        Mix is the mixing strength of the true attribute and is a number between 0 and 1.
        Mix 0 produces the true attribute, mix 1 produces the mixing attribute.
        '''
        self.true_attribute = true_attribute
        self.mixing_attribute = mixing_attribute
        self.mix = mix

        self.cell_type = f"{true_attribute.cell_type}_Mix_{mix}_{mixing_attribute.cell_type}"
        self.size = (1-mix)*true_attribute.size + mix*mixing_attribute.size
        self.scale = tuple((1-mix)*true_attribute.scale[i] + mix*mixing_attribute.scale[i] for i in range(3))
        self.deformation_strength = (1-mix)*true_attribute.deformation_strength + mix*mixing_attribute.deformation_strength
        self.attribute_name = f"{1-mix} {true_attribute.attribute_name} + {mix} {mixing_attribute.attribute_name}"
        self.max_bending_strength = (1-mix)*true_attribute.max_bending_strength + mix*mixing_attribute.max_bending_strength

class Cell:
    cell_count = 0

    def __init__(self, location: Vector, arrangement_id: int, arrangement_type: str, cell_attributes: CellAttribute, orientation: Optional[Vector] = None):
        self.cell_attributes = cell_attributes
        self.cell_id = Cell.cell_count
        Cell.cell_count += 1
        self.location = location
        self.orientation = orientation
        self.cell_name = f"Cell_{self.cell_id}_Arr_{arrangement_type}_{arrangement_id}_Type_{self.cell_attributes.cell_type}"
        self.pass_index = self.cell_id
        self.arrangement_type = None
        self.semantic_id = None
        self.material = None
        self.cell_object = None

    def set_material(self, material: Material):
        self.material = material

    def set_material(self, semantic_id: int):
        self.semantic_id = semantic_id

    def add(self):
        """
        Adds a sphere object to the scene with the specified radius, location, and scale.
        Sets the pass index of the sphere object to the specified pass index.
        Orients the sphere object according to the specified orientation.
        Applies two levels of subdivision surface to the sphere object.
        Deforms the mesh of the sphere object randomly.
        Bends the sphere object along the z axis.
        """
        # Create a sphere
        bpy.ops.mesh.primitive_ico_sphere_add(radius=self.cell_attributes.size, location=self.location, scale=self.cell_attributes.scale)
        self.cell_object = bpy.context.active_object
        self.cell_object.name = self.cell_name
        
        # NOTE: Is this line necessary? variable pass_index is not used. - ck
        pass_index = bpy.props.IntProperty(name="Pass Index", subtype='UNSIGNED')
        bpy.context.object.pass_index = self.pass_index
        
        # Orient
        if self.orientation:
            set_orientation(self.cell_object, self.orientation)

        # Apply two levels of subdivision surface
        modifier = self.cell_object.modifiers.new("Subsurface Modifier", "SUBSURF")
        modifier.levels = 2    

        # NOTE: Some options here for generating diverse looking cell shapes. - ck
        # 1) Use deform_mesh_old() and bend_mesh() afterwards. This is the most stable approach but the deformations don't look good.
        # 2) Use deform_mesh() and bend_mesh() afterwards. Not stable, Blender crashes after generating a random number of objects. But the deformations look good.
        # 3) Use a different approach using Voronoi diagrams. TBD.

        # Deform the sphere mesh randomly
        self.deform_mesh()
        # Bend mesh along the z axis
        self.bend_mesh()



    def deform_mesh(self):
        """
    	Deforms the mesh by randomly translating a subset of vertices using proportional edit.
    	That is, neighboring vertices are also translated proportionally

    	Parameters:
    	- None
    	
    	Return:
    	- None
    	
    	Internal Variables:
    	- TRANSLATION_RANGE: The maximum range of translation for each vertex.
    	- TRANSFORM_COUNT: The number of transformations to apply.
    	- PROPORTIONAL_SIZE: The size of the proportional edit range.
    	"""
        # TODO: Put these variables as members to CellAttributes class. - ck
        # NOTE: One can fine tune these values for better results. So far this is good enough. - ck
        TRANSLATION_RANGE = 0.05
        TRANSFORM_COUNT = 3
        PROPORTIONAL_SIZE = 0.2     

        # apply subsurface modifier
        bpy.ops.object.modifier_apply({"object": bpy.context.object}, modifier="Subsurface Modifier")
        # extract mesh
        mesh = self.cell_object.data
        # deselect all faces
        mesh.polygons.foreach_set("select", (False,) * len(mesh.polygons))
        # deselect all edges
        mesh.edges.foreach_set("select", (False,) * len(mesh.edges))
        # deselect all vertices
        mesh.vertices.foreach_set("select", (False,) * len(mesh.vertices))
        # translate random mesh vertices using proportional edit
        for _ in range(TRANSFORM_COUNT):
            transform = Vector([random.uniform(-1, 1),
                        random.uniform(-1, 1),
                        random.uniform(-1, 1)])*TRANSLATION_RANGE
            v = mesh.vertices[random.randint(0, len(mesh.vertices) - 1)]
            v.select = True
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.transform.translate(value=transform, 
                                    constraint_axis=(False, False, False),
                                    orient_type='GLOBAL',
                                    mirror=False, 
                                    use_proportional_edit = True,
                                    use_proportional_connected =True,
                                    proportional_edit_falloff='SMOOTH',
                                    proportional_size=PROPORTIONAL_SIZE)
            bpy.ops.object.mode_set(mode='OBJECT')

    
    def deform_mesh_old(self):
        """
        Deforms the mesh by adding a random displacement to each vertex position.
        This function iterates through each vertex of the mesh and deforms its position
        by adding a random displacement vector. The displacement vector is calculated
        by multiplying a random vector with values between -1 and 1 by the scale and
        deformation strength attributes of the cell. The original position of each 
        vertex is stored and then updated by adding the deformation vector.
        """
        # Iterate through each vertex and deform its position
        for vertex in self.cell_object.data.vertices:
            original_position = vertex.co.copy()
            deformation_vector = Vector([
                random.uniform(-1, 1),
                random.uniform(-1, 1),
                random.uniform(-1, 1)
            ])*Vector(self.cell_attributes.scale)*self.cell_attributes.deformation_strength
            vertex.co = original_position + deformation_vector


    def bend_mesh(self):
        """
        Bend mesh along Z axis.

        Parameters:
        - None

        Returns:
        - None
        """
        # Bend mesh along Z axis
        modifier = self.cell_object.modifiers.new("Simple Deform Modifier", "SIMPLE_DEFORM")
        modifier_index = len(self.cell_object.modifiers) - 1  # Index of the last added modifier
        modifier = self.cell_object.modifiers[modifier_index]
        modifier.deform_method = 'BEND'

        # Create an empty
        bpy.ops.object.empty_add(type='ARROWS', align='WORLD', location=self.cell_object.location, scale=(1, 1, 1))
        empty = bpy.context.active_object

        # Set the origin and deform axis
        modifier.origin = empty
        modifier.deform_axis = 'Z'
        bending_strength = random.uniform(-1,1)*self.cell_attributes.max_bending_strength
        modifier.angle = 2*math.pi*bending_strength

        # Remove empty object
        bpy.data.objects.remove(empty, do_unlink=True)