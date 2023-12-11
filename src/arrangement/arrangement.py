import random

from mathutils import Vector

from src.objects.cells import Cell

class CellList:
    def __init__(self, cell_attributes, locations, name = None):
        self.cell_attributes = cell_attributes
        self.locations = locations
        self.name = name
        self.objects = []

    def generate_cells(self):
        for idx, location in enumerate(self.locations):
            # TOOD: Make sure the ids are not just from tha same range in different lists. They need to be unique for each scene. - ck
            cell = Cell(idx, self.cell_attributes, location, self.name)
            self.objects.append(cell)

    def add_to_scene(self):
        for cell in self.objects:
            cell.add()


# TODO: Extend to different ditribution types, e.g. Gaussian distribution etc.
class CellDistribution:
    def __init__(self, cell_attributes, num_cells, min_coords, max_coords, name = None):
        self.cell_attributes = cell_attributes
        self.num_cells = num_cells
        self.min_coords = min_coords
        self.max_coords = max_coords
        self.name = name
        self.objects = []
        
    def generate_cells(self):
        for idx in range(self.num_cells):
            # Sample uniformly distributed location
            location = Vector([
                random.uniform(self.min_coords.x, self.max_coords.x),
                random.uniform(self.min_coords.y, self.max_coords.y),
                random.uniform(self.min_coords.z, self.max_coords.z)
            ])
            cell = Cell(idx, self.cell_attributes, location, self.name)
            self.objects.append(cell)

    def add_to_scene(self):
        for cell in self.objects:
            cell.add()