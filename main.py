# IMPORT SOURCES
dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir )

from src.arrangement.arrangement import CellArrangement
from src.objects.cells import CellAttributeEndothelia
from src.objects.cells import Cell
from src.objects.tissue import Tissue
from src.shading.shading import Material
from src.scene import BioMedicalScene, Camera, LightSource

# this next part forces a reload in case you edit the source after you first start the blender session
import imp
imp.reload(CellArrangement)
imp.reload(CellAttributeEndothelia)
imp.reload(Cell)
imp.reload(Tissue)
imp.reload(Material)
imp.reload(BioMedicalScene)
imp.reload(Camera)
imp.reloadLightSource)

###################  PARAMETER  #####################
# args_camera = {'pos'}


###################  MAIN  METHOD  #####################
# create the necessary objects

# create cell arrangements
# cells_a = CellArrangement()
# cell_list = cells_a.cells

# add microscope objects
my_tissue = Tissue(Material())
my_light_source = LightSource()
my_camera = Camera()

# create scene
my_scene = BioMedicalScene(my_light_source, my_camera)#, cell_list)

# render scene
my_scene.render()



# scn = bpy.context.scene
# cam1 = bpy.data.cameras.new("Camera 1")
# cam1.lens = 18

# # create the first camera object
# cam_obj1 = bpy.data.objects.new("Camera 1", cam1)
# cam_obj1.location = (0, 0, 2)
# cam_obj1.rotation_euler = (0, 0, 0)
# scn.collection.objects.link(cam_obj1)