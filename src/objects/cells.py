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
    EOS = 3
    FIB = 4
    EPI = 5
    GOB = 6


class CellAttribute():
    def __init__(self):
        self.cell_type = None
        self.size = None
        self.scale = None
        self.deformation_strength = None
        self.attribute_name = None
        self.max_bending_strength = None

    def from_type(cell_type):
        assert cell_type.name in cell_type_to_attribute.keys(), f"Unknown cell type {cell_type.name}."
        return cell_type_to_attribute[cell_type.name]


# NOTE: The size and scale convention for CellAttributes is as follows:
# - size: is the maximal ellipsoid radius of the nucleus 
# - scale: is the ellipsoid scale of the nucleus
# The scale of an attribute should always be normalized such the maximal scale value is 1. 
# This ensures that the maximal diameter of a cell nucleus is twice its size.
class PLA(CellAttribute):
    def __init__(self, 
                 cell_type = CellType.PLA,
                 size = 0.06,
                 scale = (1,1,1),
                 deformation_strength = 0.8,
                 attribute_name = "Plasma Cell",
                 max_bending_strength = 0.2):
        super().__init__()
        self.cell_type = cell_type
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class LYM(CellAttribute):
    def __init__(self, 
                 cell_type = CellType.LYM,
                 size = 0.03,
                 scale = (1,1,1),
                 deformation_strength = 0.8,
                 attribute_name = "Lymphocyte",
                 max_bending_strength = 0.2):
        self.cell_type = cell_type
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class EOS(CellAttribute):
    def __init__(self, 
                 cell_type = CellType.EOS,
                 size = 0.06,
                 scale = (1,0.3,0.3),
                 deformation_strength = 0.8,
                 attribute_name = "Eosinophile",
                 max_bending_strength = 0.2):
        self.cell_type = cell_type
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class FIB(CellAttribute):
    def __init__(self, 
                 cell_type = CellType.FIB,
                 size = 0.08,
                 scale = (1,0.2,0.2),
                 deformation_strength = 0.8,
                 attribute_name = "Fibroblast",
                 max_bending_strength = 0.2):
        self.cell_type = cell_type
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class EPI(CellAttribute):
    def __init__(self, 
                 cell_type = CellType.EPI,
                 size = 0.1,
                 scale = (1, 0.6, 0.6),
                 deformation_strength = 0.2,
                 attribute_name = "Epithelial Cell",
                 max_bending_strength = 0.3):
        self.cell_type = cell_type
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class GOB(CellAttribute):
    def __init__(self,
                 cell_type = CellType.GOB,
                 size = 0.03,
                 scale = (1,1,1),
                 deformation_strength = 0.8,
                 attribute_name = "Goblet Cell",
                 max_bending_strength = 0.2):
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

cell_type_to_attribute = {
    CellType.PLA.name : PLA(),
    CellType.FIB.name : FIB(),
    CellType.EOS.name : EOS(),
    CellType.LYM.name : LYM(),
    CellType.EPI.name : EPI(),
    CellType.GOB.name : GOB()
}