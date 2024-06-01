import bpy 
import bmesh
import random
import math

from enum import Enum
from mathutils import Vector
from typing import Optional

from src.shading.materials import Material
from src.utils.geometry import *

class CellType(Enum):
    PLA = 1
    LYM = 2
    EOS = 3
    FIB = 4
    EPI = 5
    GOB = 6
    MIX = 7

TYPE_MIXING = 0.3


class CellAttribute():
    def __init__(self):
        self.cell_type = None
        self.size = None
        self.attribute_name = None
        self.subdivision_levels = 1

    def add_cell_objects(self, location, direction):
        '''
        Returns a 1- or 2-element list of blender objects representing the nucleus and optionally cytoplasm.
        The cytoplasm is optional and depends on the cell type.
        Currently only PLA and EOS cells have cytoplasm.
        '''
        pass

    def from_type(cell_type):
        assert cell_type.name in cell_type_to_attribute.keys(), f"Unknown cell type {cell_type.name}."
        return cell_type_to_attribute[cell_type.name]


# NOTE: The size and scale convention for CellAttributes is as follows:
# - size: is the maximal ellipsoid radius of the nucleus 
# - scale: is the ellipsoid scale of the nucleus
# The scale of an attribute should always be normalized such the maximal scale value is 1. 
# This ensures that the maximal diameter of a cell nucleus is twice its size.
class PLA(CellAttribute):
    def __init__(self):
        super().__init__()
        self.cell_type = CellType.PLA
        self.size = 0.075 # 0.9
        self.nucleus_size = 0.05
        self.scale = (1,0.8,0.7) # (1, 0.6, 0.5)
        self.deformation_strength = 0.4 # 0.7
        self.attribute_name = "Plasma Cell"
        self.max_bending_strength = 0.2
    
    def add_cell_objects(self, location, direction):
        # Add cytoplasm
        bpy.ops.mesh.primitive_ico_sphere_add(radius=self.size)
        cytoplasm = bpy.context.active_object
        deform_mesh(cytoplasm, self)
        subdivide(cytoplasm, self.subdivision_levels)
        set_orientation(cytoplasm, direction)
        cytoplasm.location = location
        cytoplasm.scale = self.scale

        # Add nucleus
        bpy.ops.mesh.primitive_ico_sphere_add(radius=self.nucleus_size)
        nucleus = bpy.context.active_object
        deform_mesh(nucleus, self)
        subdivide(nucleus, self.subdivision_levels)
        set_orientation(nucleus, direction)
        nucleus.location = location
        nucleus.scale = self.scale
        return [nucleus, cytoplasm] # NOTE: Important: Nucleus needs to be first in list. - ck

class LYM(CellAttribute):
    def __init__(self):
        super().__init__()
        self.cell_type = CellType.LYM
        self.size = 0.04 # 0.04,
        self.nucleus_size = 0.025
        self.scale = (1, 0.9, 0.8) #(1, 0.9, 0.8)
        self.deformation_strength = 0.2 #0.7
        self.attribute_name = "Lymphocyte"
        self.max_bending_strength = 0.2
    
    def add_cell_objects(self, location, direction):
        bpy.ops.mesh.primitive_ico_sphere_add(radius=self.size)
        nucleus = bpy.context.active_object
        deform_mesh(nucleus, self)
        subdivide(nucleus, self.subdivision_levels)
        set_orientation(nucleus, direction)
        nucleus.location = location
        nucleus.scale = self.scale
        return [nucleus]

class EOS(CellAttribute):
    def __init__(self):
        super().__init__()
        self.cell_type = CellType.EOS
        self.size = 0.056
        self.nucleus_size = 0.025
        self.scale = (1,1,1)
        self.deformation_strength = 0.2
        self.attribute_name = "Eosinophil"
        self.max_bending_strength = 0.2

    def add_cell_objects(self, location, direction):
        # Add cytoplasm
        bpy.ops.mesh.primitive_ico_sphere_add(radius=self.size)
        cytoplasm = bpy.context.active_object
        deform_mesh(cytoplasm, self)
        subdivide(cytoplasm, self.subdivision_levels)
        set_orientation(cytoplasm, direction)
        cytoplasm.location = location
        cytoplasm.scale = self.scale

        # Add nucleus
        coeff = 0.55 # Controls the displacement between to metaballs, 0: no distance
        delta = 2.0*coeff*self.nucleus_size
        rad1 = self.size*(1-coeff) + self.deformation_strength*self.size*random.uniform(-1,1)
        rad2 = 2.0*self.size*(1-coeff) - rad1
        assert rad1 > 0 and rad2 > 0, "Negative radius."
        center = self.size - rad1
        bpy.ops.object.metaball_add(type='BALL', radius=rad1, enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
        bpy.ops.object.metaball_add(type='BALL', radius=rad2, enter_editmode=False, align='WORLD', location=(delta, 0, 0), scale=(1, 1, 1))
        bpy.ops.object.convert(target='MESH')
        bpy.context.scene.cursor.location = (center, 0, 0)
        bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
        nucleus = bpy.context.active_object
        set_orientation(nucleus, direction)
        nucleus.location = location
        nucleus.scale = self.scale
        return [nucleus, cytoplasm]

class FIB(CellAttribute):
    def __init__(self):
        super().__init__()
        self.cell_type = CellType.FIB
        self.size = 0.1
        self.nucleus_size = 0.1
        self.scale = (1,0.6,0.5)
        self.deformation_strength = 0.2
        self.attribute_name = "Fibroblast"
        self.max_bending_strength = 0.7

    def add_cell_objects(self, location, direction):
        coeff = 0.3 # Controls the displacement between to metaballs, 0: no distance
        delta = 2.0*coeff*self.size
        rad1 = self.size*(1-coeff) + self.deformation_strength*self.size*random.uniform(-1,1)
        rad2 = 2.0*self.size*(1-coeff) - rad1
        assert rad1 > 0 and rad2 > 0, "Negative radius."
        center = self.size - rad1
        bpy.ops.object.metaball_add(type='ELLIPSOID', radius=rad1, enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
        bpy.ops.object.metaball_add(type='ELLIPSOID', radius=rad2, enter_editmode=False, align='WORLD', location=(delta, 0, 0), scale=(1, 1, 1))
        bpy.ops.object.convert(target='MESH')
        bpy.context.scene.cursor.location = (center, 0, 0)
        bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
        nucleus = bpy.context.active_object
        bend_mesh(nucleus, self.max_bending_strength)
        set_orientation(nucleus, direction)
        nucleus.location = location
        nucleus.scale = self.scale
        return [nucleus]

class EPI(CellAttribute):
    def __init__(self):
        super().__init__()
        self.cell_type = CellType.EPI
        self.size = 0.04
        self.attribute_name = "Epithelial Cell"
        self.smooth_factor = 2 # (float 1) Controls the roundness of the object. 0 = identical to surrounding mesh, the higher this number, the rounder the mesh.
        self.smooth_roundness = 2 # (int, 2) Controls the size of objects, the higher this number, the smaller the mesh. Should be at least 2.

    def add_cell_objects(self, location, direction):
        pass

class GOB(CellAttribute):
    def __init__(self):
        super().__init__()
        self.cell_type = CellType.GOB
        self.size = 0.1
        self.attribute_name = "Goblet Cell"
        self.smooth_factor = 1.5 # (float, 1) Controls the roundness of the object. 0 = identical to surrounding mesh, the higher this number, the rounder the mesh.
        self.smooth_roundness = 2 # (int, 2) Controls the size of objects, the higher this number, the smaller the mesh. Should be at least 2.
        self.subdivision_levels = 2 

    def add_cell_objects(self, location, direction):
        pass



class MixAttribute(CellAttribute):
    def __init__(self, true_attribute: CellAttribute, mixing_attribute: CellAttribute, mix: float):
        '''
        Produces a cell attribute that is a combination of two cell attributes in terms of numerical values.
        Mix is the mixing strength of the true attribute and is a number between 0 and 1.
        Mix 0 produces the true attribute, mix 1 produces the mixing attribute.
        '''
        super().__init__()
        self.true_attribute = true_attribute
        self.mixing_attribute = mixing_attribute
        self.mix = mix

        # TODO: Add lerp function
        self.cell_type = CellType.MIX
        self.size = lerp(true_attribute.size, mixing_attribute.size, mix)
        self.nucleus_size = lerp(true_attribute.nucleus_size, mixing_attribute.nucleus_size, mix)
        self.scale = tuple(lerp(true_attribute.scale[i], mixing_attribute.scale[i], mix) for i in range(3))
        self.deformation_strength = lerp(true_attribute.deformation_strength, mixing_attribute.deformation_strength, mix)
        self.attribute_name = f"{1-mix} {true_attribute.attribute_name} + {mix} {mixing_attribute.attribute_name}"
        self.max_bending_strength = lerp(true_attribute.max_bending_strength, mixing_attribute.max_bending_strength, mix)

    def lerp(a, b, t):
        return a*(1-t) + b*t
    
    def add_cell_objects(self, location, direction):
            bpy.ops.mesh.primitive_ico_sphere_add(radius=self.nucleus_size)
            nucleus = bpy.context.active_object
            deform_mesh(nucleus, self)
            subdivide(nucleus, self.subdivision_levels)
            set_orientation(nucleus, direction)
            nucleus.location = location
            nucleus.scale = self.scale
            # TODO: Should a PLA-LYM mix type have cytoplasm? If yes, need to implement it.
            return [nucleus]

cell_type_to_attribute = {
    CellType.PLA.name : PLA(),
    CellType.LYM.name : LYM(),
    CellType.EOS.name : EOS(),
    CellType.FIB.name : FIB(),
    CellType.EPI.name : EPI(),
    CellType.GOB.name : GOB(),
    CellType.MIX.name : MixAttribute(PLA(), LYM(), TYPE_MIXING),
}

def initialize_mixing_attribute(mixing_coefficient):
    '''
    Initializes the mixing attribute with the given mixing coefficient.
    '''
    cell_type_to_attribute[CellType.MIX.name] = MixAttribute(PLA(), LYM(), mixing_coefficient)
    