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
import src.objects.macro_structures as macro
import src.objects.tissue_architecture as arch

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
imp.reload(macro)
imp.reload(arch)

###################  PARAMETER  #####################
# args_camera = {'pos'} # no change just test

# NOTE: It might be good to use realistic units here.
# Tissue usually 5-10 mu thick, Epithelial cells usally 9-17 mu thick, nuclei apparently also 5-10 mu thick.
TISSUE_THICKNESS = 0.2
TISSUE_SIZE = 2
TISSUE_LOCATION = (0, 0, 0.5)
TISSUE_PADDING = 0.5


###################  MAIN  METHOD  #####################
# create the necessary objects
scene.BioMedicalScene.clear()
    
# add microscope objects
my_materials = shading.Material()
my_tissue = tissue.Tissue(my_materials.muscosa, thickness=TISSUE_THICKNESS, size=TISSUE_SIZE, location=TISSUE_LOCATION) # thickness and location of tissue should encapsulate min and max z-coordinates of cells 
my_light_source = scene.LightSource(material=my_materials.light_source)
my_camera = scene.Camera()

# create scene
my_scene = scene.BioMedicalScene(my_light_source, my_camera)

# Define volume and surface objects
# NOTE: In the end the volume and surface objects should come from the epithelial tissue macrostructure. - ck
vol_scale = (0.6, 0.3, 1)
surf_scale = (0.8, 0.5, 1)

# Add tissue
my_scene.add_tissue(tissue=my_tissue.tissue)

# add macrostructure
tissue_arch = arch.TissueArch()
tissue_arch.random_crop(my_tissue.tissue)
SURF_OBJ, VOL_OBJ = tissue_arch.get_architecture()
my_scene.bound_architecture(volumes=[VOL_OBJ], surfaces=[SURF_OBJ])

# VOL_OBJ, SURF_OBJ = utils.geometry.add_dummy_objects(my_tissue, TISSUE_PADDING, vol_scale, surf_scale)

# NOTE: For some very weird reason you need to create the surface filling before the volume filling.
# Otherwise the surface filling won't work and it won't even refine the mesh. :/ - ck
# TODO: Fix that
# Add surface filling
# SURF_NUMBER = 20
# SURF_ATTRIBUTE = cells.CellAttributeEpi(size=0.1, scale=(1, 0.5, 0.5))
# FILLER_SCALE = 0.8 # Scale of the size of smaller filler nuclei w.r.t to the original nuclei size
#surface_fill = arr.SurfaceFill(SURF_OBJ, SURF_NUMBER, SURF_ATTRIBUTE, FILLER_SCALE)

# Add volume filling
NUMBER = 20
ATTRIBUTES = [cells.CellAttributeA(), cells.CellAttributeB(), cells.CellAttributeC()]
RATIOS = [0.05, 0.15, 0.8]
volume_fill = arr.VolumeFill(VOL_OBJ, NUMBER, ATTRIBUTES, RATIOS, strict_boundary=False)
my_scene.add_arrangement(volume_fill)
#my_scene.add_arrangement(surface_fill)

# FINAL CUTTING
my_scene.cut_cells()
my_scene.cut_tissue()
my_scene.add_tissue_staining(materials=[my_materials.muscosa])
my_scene.add_staining(material=my_materials.nuclei_staining)

# render scene
#RENDER_PATH = 'renders/'

#my_scene.render(filepath = RENDER_PATH,  # where to save renders
#               scene = True, # if true scene is rendered
#               single_masks = True, # if true singel cell masks are rendered
#               semantic_mask = True, # if true semantic mask is generated
#               instance_mask = True, # if true instance mask is generated
#               depth_mask = True, # if true depth mask is generated
#               obj3d = True, # if true scene is saved as 3d object
#               output_shape = (500, 500), # dimensions of output
#               max_samples = 10) # number of samples for rendering. Fewer samples will render more quickly. Default is 1024