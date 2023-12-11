
import abc

from src.shading.shading import Material


class CellAttribute(metaclass=abc.ABCMeta):
    def __init__(self):
        self.cell_type = None
        self.cell_id = None
        self.cell_name = None
        self.shape = None
        self.orientation = None
        self.elongation = None



class CellAttributeEndothelia(CellAttribute):
    def __init__(self,cell_id,  cell_name, cell_type = "Endothelia", shape="stretched", orientation="vector", elongation=0.1, material=Material()):
        self.cell_type = cell_type
        self.cell_id = cell_id
        self.cell_name = cell_name
        self.shape = shape
        self.orientation = orientation
        self.elongation = elongation
        self.material = material



class Cell:
    def __init__(self, cell_attributes: CellAttribute, semantic_id: int, cell_name: str, material: Material):
        self.cell_attributes = cell_attributes
        self.semantic_id = semantic_id
        self.cell_name = cell_name
        self.material = material

    def generate(self):
        pass


