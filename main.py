from mathutils import Vector

# IMPORT SOURCES
import rHEnder.arrangement.arrangement as arr
import rHEnder.objects.cells as cells
import rHEnder.objects.tissue as tissue
import rHEnder.shading.shading as shading
import rHEnder.scene as scene

###################  PARAMETER  #####################
# args_camera = {'pos'}


###################  MAIN  METHOD  #####################
# create the necessary objects
scene.BioMedicalScene.clear()
    
# add microscope objects
my_tissue = tissue.Tissue(shading.Material())
my_light_source = scene.LightSource()
my_camera = scene.Camera()

# create scene
my_scene = scene.BioMedicalScene(my_light_source, my_camera)

# define cell arrangements
cell_distribution = arr.CellDistribution(
    cell_attributes = cells.CellAttributeA(),
    num_cells = 200,
    min_coords = Vector([-1, -1, 0.4]),
    max_coords = Vector([1, 1, 0.6])
)
locations = [Vector([-0.5, -0.5, 0.5]), Vector([0.5, -0.5, 0.5]), Vector([-0.5, 0.5, 0.5]), Vector([0.5, 0.5, 0.5])]
cell_list = arr.CellList(
    cell_attributes = cells.CellAttributeB(),
    locations = locations
)

# add cell arrangements to scene
my_scene.add_arangement(cell_distribution)
my_scene.add_arangement(cell_list)


# render scene
my_scene.render()
