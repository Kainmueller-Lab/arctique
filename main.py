import bpy
import random
import sys
import time
import os

TYPE_MIXING = 0.3

# IMPORT SOURCES
dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir )

import src.arrangement.arrangement as arr 
import src.objects.cells as cells
import src.objects.tissue as tissue
import src.shading.materials as materials
import src.scene as scene
import src.utils as utils
import src.utils.geometry as geom
import src.utils.surface_filling as sf
import src.utils.volume_filling as vf
import src.objects.tissue_architecture as arch
import src.utils.helper_methods as hm

from src.objects.cells import CellType

# this next part forces a reload in case you edit the source after you first start the blender session
#import imp
import importlib as imp # imp module is deprecated since python 3.12
imp.reload(arr)
imp.reload(cells)
imp.reload(tissue)
imp.reload(materials)
imp.reload(scene)
imp.reload(utils)
imp.reload(geom)
imp.reload(sf)
imp.reload(vf)
imp.reload(hm)

###################  PARAMETER  #####################
# args_camera = {'pos'} # no change just test

# NOTE: It might be good to use realistic units here.
# Tissue usually 5-10 mu thick, Epithelial cells usally 9-17 mu thick, nuclei apparently also 5-10 mu thick.
TISSUE_THICKNESS = 0.05
TISSUE_SIZE = 1.28
TISSUE_LOCATION = (0, 0, 0.5)
TISSUE_PADDING = 0.5
SEED = random.randint(0, 10000)

MIX_COUNT = 100
RATIOS = [0.1, 0.3, 0.4, 0.1, 0.1]
MIX_TYPES = [CellType.MIX, CellType.PLA, CellType.LYM, CellType.EOS, CellType.FIB] # NOTE: MIX is a PLA nucleus that has a mixed shape interpolated to LYM. - ck
# NOTE: TYPE_MIXING is an unused dummy variable. You have to set it in the cells.py file.
# Once we use config files this should easily be solvable. - ck
TYPE_MIXING = 0.3 # Mixing strength of PLA nuclei (between 0 and 1, 0 is a pure PLA shape, 1 is a pure LYM shape)


###################  MAIN  METHOD  #####################
# create the necessary objects
scene.BioMedicalScene.clear()
    
 # 1) initialize microscope objects and add to scene
my_materials = materials.Material(seed=SEED)
my_tissue = tissue.Tissue(
    my_materials.muscosa, thickness=TISSUE_THICKNESS,
    size=TISSUE_SIZE, location=TISSUE_LOCATION)
my_light_source = scene.LightSource(material=my_materials.light_source)
my_camera = scene.Camera()
my_scene = scene.BioMedicalScene(my_light_source, my_camera)
my_scene.add_tissue(tissue=my_tissue.tissue)

# 2) create macrostructures in tissue block, rotate and scale them and cut them
tissue_arch = arch.TissueArch(seed=SEED)
tissue_arch.random_crop(my_tissue.tissue)
macro_structure = tissue_arch.get_architecture()
crypt, epi_volume, crypt_vol_2, mucosa = macro_structure
my_scene.bound_architecture(
        volumes=[epi_volume, crypt_vol_2, mucosa], surfaces=[crypt],
        padding=TISSUE_PADDING)
vol_goblet = hm.copy_object(crypt_vol_2, 'vol_goblet')
extended_stroma = hm.copy_object(mucosa, 'extended_stroma')
hm.add_boolean_modifier(extended_stroma, epi_volume, operation='UNION', name='add epi to stroma', apply=True)
hm.add_boolean_modifier(vol_goblet, extended_stroma, name='Remove inner volume', apply=True)

# 3) populate scene with nuclei/cells
# add mix volume filling
start = time.time()
#volume_fill = arr.VolumeFill(mucosa, MIX_COUNT, MIX_TYPES, RATIOS, strict_boundary=True, seed=SEED)
end1 = time.time()
print(f"Volume filling took {end1 - start} s")
#my_scene.add_arrangement(volume_fill) # NOTE: 240 nuclei take about 20 s
end2 = time.time()
print(f"Volume adding took {time.time() - end1} s")

# add epi volume filling
epi_fill = arr.VoronoiFill(epi_volume, mucosa, CellType.EPI)
goblet_fill = arr.VoronoiFill(vol_goblet, extended_stroma, CellType.GOB)
end3 = time.time()
print(f"Voronoi filling took {end3 - end2} s")
my_scene.add_arrangement(epi_fill) # NOTE: 200 nuclei take about 40 s
my_scene.add_arrangement(goblet_fill)
end4 = time.time()
print(f"Voronoi adding took {end4 - end3} s")

# 4) cut objects and add staining
#my_scene.cut_cells()
#my_scene.cut_tissue()
# my_scene.add_tissue_staining(materials=[my_materials.muscosa, my_materials.crypt_staining])
# my_scene.add_staining(material=my_materials.nuclei_mask)
# my_scene.add_staining(material=my_materials.nuclei_staining)

# Hide non cell objects
my_scene.hide_non_cell_objects()

# # render scene
# RENDER_PATH = 'C:/Users/cwinklm/Documents/Alpacathon/rendered_HE/renders2d_test/'
# RENDER_PATH = 'renders/'

# my_scene.render(filepath = RENDER_PATH,  # where to save renders
#               scene = True, # if true scene is rendered
#               single_masks = True, # if true singel cell masks are rendered
#               semantic_mask = True, # if true semantic mask is generated
#               instance_mask = True, # if true instance mask is generated
#               depth_mask = True, # if true depth mask is generated
#               obj3d = True, # if true scene is saved as 3d object
#               output_shape = (500, 500), # dimensions of output
#               max_samples = 10) # number of samples for rendering. Fewer samples will render more quickly. Default is 1024



# my_scene.render3d(filepath = RENDER_PATH,  # where to save renders
#                scene = True, # if true scene is rendered
#                n_slices = 10,
#                semantic_mask = True, # if true semantic mask is generated
#                instance_mask = True, # if true instance mask is generated
#                depth_mask = True, # if true depth mask is generated
#                obj3d = True, # if true scene is saved as 3d object
#                output_shape = (500, 500), # dimensions of output
#                max_samples = 10) # number of samples for rendering. Fewer samples will render more quickly. Default is 1024