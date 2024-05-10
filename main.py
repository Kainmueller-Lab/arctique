import bpy
import sys
import os


# IMPORT SOURCES
dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir )

import src.arrangement.arrangement as arr 
import src.objects.cells as cells
import src.objects.tissue as tissue
import src.shading.shading as shading
import src.scene as scene
import src.utils as utils
import src.utils.geometry as geom
import src.utils.surface_filling as sf

# this next part forces a reload in case you edit the source after you first start the blender session
#import imp
import importlib as imp # imp module is deprecated since python 3.12
imp.reload(arr)
imp.reload(cells)
imp.reload(tissue)
imp.reload(shading)
imp.reload(scene)
imp.reload(utils)
imp.reload(geom)
imp.reload(sf)

###################  PARAMETER  #####################
# args_camera = {'pos'} # no change just test

# NOTE: It might be good to use realistic units here.
# Tissue usually 5-10 mu thick, Epithelial cells usally 9-17 mu thick, nuclei apparently also 5-10 mu thick.
TISSUE_THICKNESS = 0.2
TISSUE_SIZE = 2
TISSUE_LOCATION = (0, 0, 0.5)
TISSUE_PADDING = 0.04


###################  MAIN  METHOD  #####################
# create the necessary objects
scene.BioMedicalScene.clear()
    
# add microscope objects
my_materials = shading.Material()

my_tissue = tissue.Tissue(my_materials.tissue_staining, thickness=TISSUE_THICKNESS, size=TISSUE_SIZE, location=TISSUE_LOCATION) # thickness and location of tissue should encapsulate min and max z-coordinates of cells 
my_light_source = scene.LightSource(material=my_materials.light_source)
my_camera = scene.Camera()

# create scene
my_scene = scene.BioMedicalScene(my_light_source, my_camera)

# add bounding volumes
# NOTE: MIX_VOL bounds the volume for the mixed cell types
# NOTE: EPI_VOL bounds the volume for the epithelial cell types.
MIX_VOL, EPI_VOL = utils.geometry.add_dummy_volumes(my_tissue, TISSUE_PADDING)

# add mix volume filling
MIX_COUNT = 240
# NOTE: Create nuclei of type A which are mixed with nuclei of type C with a factor of 0.3.
# A mix factor of 0 produces the pure true attribute, mix factor 1 produces the pure mixing attribute.
ATTRIBUTES = [cells.MixAttribute(cells.CellAttributeA(), cells.CellAttributeB(), 0.3), cells.CellAttributeA(), cells.CellAttributeB(), cells.CellAttributeC()]
RATIOS = [0.2, 0.2, 0.2, 0.4]
volume_fill = arr.VolumeFill(MIX_VOL, MIX_COUNT, ATTRIBUTES, RATIOS, strict_boundary=True)
my_scene.add_arrangement(volume_fill) # NOTE: 240 nuclei take about 15 s

add epi volume filling
EPI_COUNT = 200
EPI_ATTRIBUTE = cells.CellAttributeEpi(size=0.1, scale=(1, 0.5, 0.5))
crypt_fill = arr.VoronoiFill(EPI_VOL, EPI_COUNT, EPI_ATTRIBUTE)
my_scene.add_arrangement(crypt_fill) # NOTE: 200 nuclei take about 30 s

# Add tissue
my_scene.add_tissue(tissue=my_tissue.tissue)
my_scene.cut_cells()
my_scene.add_staining(material=my_materials.nuclei_staining)

# render scene
RENDER_PATH = 'C:/Users/cwinklm/Documents/Alpacathon/rendered_HE/renders/'
#RENDER_PATH = 'renders/'

my_scene.render(filepath = RENDER_PATH,  # where to save renders
               scene = True, # if true scene is rendered
               single_masks = True, # if true singel cell masks are rendered
               semantic_mask = True, # if true semantic mask is generated
               instance_mask = True, # if true instance mask is generated
               depth_mask = True, # if true depth mask is generated
               obj3d = True, # if true scene is saved as 3d object
               output_shape = (500, 500), # dimensions of output
               max_samples = 10) # number of samples for rendering. Fewer samples will render more quickly. Default is 1024








