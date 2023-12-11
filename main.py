# a usecase how to use this library


from src.arrangement.arrangement import CellArrangement
from src.objects.cells import CellAttributeEndothelia
from src.objects.cells import Cell

from src.objects.tissue import Tissue
from src.shading.shading import Material

from src.scene import BioMedicalScene, Camera, LightSource

# create the necessary objects
my_end_cell_attr = CellAttributeEndothelia(cell_id=1, cell_name="mycell")
my_end_cell = Cell(my_end_cell_attr, cell_name="myCell", semantic_id=1, material=Material())

my_tissue = Tissue(Material())

my_light_source = LightSource()
my_camera = Camera()

# create scene
my_scene = BioMedicalScene(my_light_source, my_camera)

# render
my_scene.render()
