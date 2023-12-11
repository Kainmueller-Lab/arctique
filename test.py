from mathutils import Vector

from src.arrangement.arrangement import CellList, CellDistribution
from src.objects.cells import CellAttributeA, CellAttributeB
from src.objects.cells import Cell

from src.objects.tissue import Tissue
from src.shading.shading import Material

from src.scene import BioMedicalScene, Camera, LightSource

# Create random cell field

cell_distribution = CellDistribution(
    cell_attributes = CellAttributeA(),
    num_cells = 20,
    min_coords = Vector([-2, -2, -2]),
    max_coords = Vector([2, 2, 2])
)
cell_distribution.generate_cells()

for cell in cell_distribution.objects:
    cell.add()

locations = [Vector([-2, -2, -2]), Vector([-2, -2, 2]), Vector([-2, 2, -2]), Vector([-2, 2, 2]), Vector([2, -2, -2]), Vector([2, -2, 2]), Vector([2, 2, -2]), Vector([2, 2, 2])]
cell_list = CellList(
    cell_attributes = CellAttributeB(),
    locations = locations
)
cell_list.generate_cells()

for cell in cell_list.objects:
    cell.add()


# # create the necessary objects
# my_end_cell_attr = CellAttributeEndothelia(cell_id=1, cell_name="mycell")
# my_end_cell = Cell(my_end_cell_attr, cell_name="myCell", semantic_id=1, material=Material())

# my_tissue = Tissue(Material())

# my_light_source = LightSource()
# my_camera = Camera()

# # create scene
# my_scene = BioMedicalScene(my_light_source, my_camera)

# # render
# my_scene.render()
