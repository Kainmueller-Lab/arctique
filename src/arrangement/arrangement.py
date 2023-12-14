import random
from mathutils import Vector
from src.objects.cells import Cell

class CellArrangement:
    # Class variable to keep track of the count
    count = 0

    def __init__(self):
        pass

    def generate_cells(self):
        pass

    def add(self):
        for cell in self.objects:
            cell.add()

class CellList(CellArrangement):
    def __init__(self, cell_attributes, locations):
        self.cell_attributes = cell_attributes
        self.locations = locations
        self.type = "LIST"
        self.id = CellArrangement.count
        CellArrangement.count += 1
        self.objects = []

    def generate_cells(self):
        for location in self.locations:
            # TODO: Make sure the ids are not just from the same range in different lists.
            # They need to be unique for each scene. - ck
            cell = Cell(location, self.id, self.type, self.cell_attributes)
            self.objects.append(cell)
        
# TODO: Extend to different distribution types, e.g., Gaussian distribution, etc.
class CellDistribution(CellArrangement):
    def __init__(self, cell_attributes, num_cells, min_coords, max_coords):
        self.cell_attributes = cell_attributes
        self.num_cells = num_cells
        self.min_coords = min_coords
        self.max_coords = max_coords
        self.type = "DIST"
        self.id = CellArrangement.count
        CellArrangement.count += 1
        self.objects = []

    def generate_cells(self):
        for _ in range(self.num_cells):
            # Sample uniformly distributed location
            location = Vector([
                random.uniform(self.min_coords.x, self.max_coords.x),
                random.uniform(self.min_coords.y, self.max_coords.y),
                random.uniform(self.min_coords.z, self.max_coords.z)
            ])
            cell = Cell(location, self.id, self.type, self.cell_attributes)
            self.objects.append(cell)
