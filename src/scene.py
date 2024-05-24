import bpy
import json
import os
from PIL import Image

import src.arrangement.arrangement as arr
import src.utils.geometry as geo
import src.utils.helper_methods as hm
import src.utils.plot_helpers as ph
import src.utils.camera_utils as cu

#import imp
import importlib as imp # imp module is deprecated since python 3.12
from pathlib import Path
import time
import numpy as np
import random
imp.reload(arr)
imp.reload(geo)
imp.reload(hm)
imp.reload(ph)

def fn_print_time_when_render_done(dummy):
    
    print("----- the time is: ", time.time())

class Camera:
    def __init__(
            self, name='camera 1', pos=(0, 0, 1.5622), rot=(0, 0, 0), size=2,
            focus_pos=(0, 0, 0.62), fstop=0.8, sensor_width=36):
        self.scene = bpy.context.scene
        
        # add camera to scene
        cam = bpy.data.cameras.new(name)
        cam.lens = 18
        self.cam_obj = bpy.data.objects.new(name, cam)
        self.cam_obj.location = pos
        self.cam_obj.rotation_euler = rot
        self.cam_obj.data.type = 'ORTHO'
        self.cam_obj.data.ortho_scale = size
        self.cam_obj.data.sensor_width = sensor_width
        self.cam_obj.data.dof.aperture_fstop = fstop

        # add focus plane to scene
        bpy.ops.mesh.primitive_plane_add(size=2.5, location=focus_pos)
        self.focus = bpy.context.active_object
        self.focus.name = 'focus'
        self.focus.hide_viewport = True
        self.focus.hide_render = True

        # set camera focus
        self.cam_obj.data.dof.use_dof = True
        self.cam_obj.data.dof.focus_object = self.focus

        # add camera for mask rendering
        cam2 = bpy.data.cameras.new(name)
        cam2.lens = 18
        self.cam_obj_mask = bpy.data.objects.new(name + '_mask', cam2)
        self.cam_obj_mask.location = pos
        self.cam_obj_mask.rotation_euler = rot
        self.cam_obj_mask.data.type = 'ORTHO'
        self.cam_obj_mask.data.ortho_scale = size

    def switch_to_mask_camera(self, scene):
        scene.camera = self.cam_obj_mask

    def switch_to_main_camera(self, scene):
        scene.camera = self.cam_obj


class LightSource:
    def __init__(self, material, name='lightsource'):
        # create mesh
        bpy.ops.mesh.primitive_plane_add(size=2, location=(0, 0, -0.2))
        self.light_source = bpy.context.active_object

        # add shading
        bpy.context.active_object.data.materials.append(material)
        bpy.context.active_object.active_material = material


class BioMedicalScene:
    def __init__(self, light_source: LightSource, camera: Camera, sample_name: int = None):
        self.light_source = light_source
        self.camera = camera
        self.sample_name = sample_name
        self.arrangements = []
        self.cell_objects = []
        self.nuclei_objects = []
        self.cytoplasm_objetcs = []
        self.scene = bpy.context.scene
        self.scene.camera = self.camera.cam_obj
        self._clear_compositor()

    @staticmethod
    def clear():
        hm.clear_scene()
        # hm.delete_objects()

    def _clear_compositor(self):
        if self.scene.node_tree is not None:
            bpy.context.scene.node_tree.nodes.clear()

    def add_tissue(self, tissue):
        self.tissue = tissue
        bpy.context.view_layer.objects.active = bpy.context.scene.objects['tissue']
        self.tissue.hide_viewport = True
        self.tissue.hide_render = True

        # add tissue empty for final slicing
        self.tissue_empty = tissue.copy()
        self.tissue_empty.data = tissue.data.copy()
        bpy.context.collection.objects.link(self.tissue_empty)
        self.tissue_empty.name = 'tissue_empty'
        self.tissue_empty.hide_viewport = True
        self.tissue_empty.hide_render = True
        self.tissue_empty.scale.z = 0.999

        # add bounding box for omitting objects outside of tissue
        self.tissue_bound = tissue.copy()
        self.tissue_bound.data = tissue.data.copy()
        bpy.context.collection.objects.link(self.tissue_bound)
        self.tissue_bound.name = 'tissue_bound'
        self.tissue_bound.hide_viewport = True
        self.tissue_bound.hide_render = True
        self.tissue_bound.scale.x = 1.4
        self.tissue_bound.scale.y = 1.4
        self.tissue_bound.scale.z = 1.4

    def bound_architecture(self, volumes=[], surfaces=[], padding=0.1):
        self.tissue_bound.dimensions.z = self.tissue.dimensions.z + padding
        self.volumes = volumes
        self.surfaces = surfaces
        for v in self.volumes:
            boolean = v.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
            boolean.operation = 'INTERSECT'
            boolean.object = self.tissue_bound
            bpy.context.view_layer.objects.active = v
            bpy.ops.object.modifier_apply(modifier=boolean.name)
        for s in self.surfaces:
            boolean = s.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
            boolean.operation = 'INTERSECT'
            boolean.object = self.tissue_bound
            bpy.context.view_layer.objects.active = s
            bpy.ops.object.modifier_apply(modifier=boolean.name)

    def cut_tissue(self):
        for v in self.volumes:
            boolean = v.modifiers.new(name="Boolean Modifier 2", type='BOOLEAN')
            boolean.operation = 'INTERSECT'
            boolean.object = self.tissue
            #bpy.ops.object.modifier_apply(modifier=boolean.name)

    def cut_cells(self, boolean_object=None):
        for cell in self.cell_objects:
            boolean = cell.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
            boolean.operation = 'INTERSECT'
            if boolean_object is None:
                boolean.object = self.tissue_empty
            else:
                boolean.object = boolean_object

    def uncut_cells(self): 
        """ reverses the action of cut_cells"""
        for cell in self.cell_objects:
            mod = cell.modifiers['Boolean Modifier']
            cell.modifiers.remove(mod)

    def cut_cells_in_tissue(self):
        for cell in self.cell_objects:
            boolean = self.tissue.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
            boolean.operation = 'DIFFERENCE'
            boolean.object = cell

    def add_tissue_staining(self, materials):
        '''
        adds materials to each macroscopic tissue volume
        by matching the material to the volume by name
        '''
        for m in materials:
            for v in self.volumes:
                if m.name in v.name:
                    v.data.materials.append(m)
                    v.active_material = m
    
    def add_staining(self, material):
        for cell in self.cell_objects:
            cell.data.materials.append(material)
            cell.active_material = material
    
    def add_arrangement(self, cell_arrangement: arr.CellArrangement):
        self.arrangements.append(cell_arrangement)
        cell_arrangement.add()
        print(f"Added arrangement {cell_arrangement.name} with {len(cell_arrangement.objects)} objects.")
        self.cell_objects = self.cell_objects + cell_arrangement.objects
        self.nuclei_objects = self.nuclei_objects + cell_arrangement.nuclei
        self.cytoplasm_objetcs = self.cytoplasm_objetcs + cell_arrangement.cytoplasm

    def rename_nuclei(self):
        for idx, cell in enumerate(self.cell_objects):
            parts = cell.name.split("_")
            cell.name = f"{parts[0]}_{idx}_{parts[1]}_{parts[2]}"

    def hide_everything(self): 
        '''hide all objects in the scene'''
        # hide tissue 
        self.tissue.hide_viewport = True
        self.tissue.hide_render = True
        # hide light source
        self.light_source.light_source.hide_viewport = True
        self.light_source.light_source.hide_render = True
        # hide cells 
        for cell in self.cell_objects: 
            cell.hide_viewport = True
            cell.hide_render = True

        self.hide_non_cell_objects()

    def unhide_everything(self): 
        '''unhide all objects in the scene'''
        # unhide tissue
        self.tissue.hide_viewport = False
        self.tissue.hide_render = False
        # unhide light source
        self.light_source.light_source.hide_viewport = False
        self.light_source.light_source.hide_render = False
        # unhide cells 
        for cell in self.cell_objects: 
            cell.hide_viewport = False
            cell.hide_render = False
        # unhide volumes and surfaces
        for v in self.volumes:
            v.hide_viewport = False
            v.hide_render = False
        for s in self.surfaces:
            s.hide_viewport = False
            s.hide_render = False

    def hide_non_cell_objects(self):
        '''hide all objects except cells in the scene'''
        for obj in self.scene.objects:
            if not obj.name.startswith(('Nucleus', 'Cytoplasm')):
                obj.hide_viewport = True
                obj.hide_render = True

    def setup_scene_render_mask(self, output_shape = (500, 500)): 
        '''specify settings for the rendering of the individual cell masks '''
        self.scene.render.resolution_x = output_shape[0]
        self.scene.render.resolution_y = output_shape[0]

        self.scene.render.engine = "BLENDER_WORKBENCH"
        self.scene.display.shading.light = "FLAT"

        self.scene.display.shading.background_type = "WORLD"
        self.scene.display.shading.single_color = (255, 255, 255)
        self.scene.display.render_aa = "OFF"

        self.scene.render.film_transparent = True # this makes the background transparent and the alpha chanell of output will have only two pixel values
        self.scene.render.image_settings.color_mode = "RGBA"

    def setup_scene_render_default(self, output_shape = (500, 500), max_samples = 1024): 
        '''specify settings for the rendering of the full scene'''
        self.scene.render.resolution_x = output_shape[0]
        self.scene.render.resolution_y = output_shape[0]
        self.scene.render.engine = "CYCLES"
        self.scene.cycles.samples = max_samples
        
        filepath = self.filepath
        if not os.path.exists(filepath):
            os.makedirs(filepath)
        
        self.scene.render.filepath = filepath + f"/images/img_{self.sample_name}.png"

    def export_scene(self): 
        '''Generates a png image of the complete scene using the specifications defined in setup_scene_render_default'''
        bpy.ops.render.render('EXEC_DEFAULT', write_still=True)
    
    def export_masks(self): 
        ''''
        create a binary mask for each individuals cell object
        '''
        self.hide_everything()
        
        filepath = self.filepath
        if not os.path.exists(filepath):
            os.makedirs(filepath)
        self.new_filepath = filepath + f"/masks/instance_individual/{self.sample_name}/"
        
        for cell_info_dict in self.cell_info: 
            cell_object = bpy.data.objects[cell_info_dict["Cellname"]] # get cell object
            cell_object.hide_viewport = False # unhide cell from viewport
            cell_object.hide_render = False # unhide cell from render
            self.scene.render.filepath = self.new_filepath + cell_info_dict["Cellname"] + '.png' #self.filepath + mask_nam
            bpy.ops.render.render('EXEC_DEFAULT', write_still=True) # render single cell mask
            cell_object.hide_viewport = True # hide cell from viewport
            cell_object.hide_render = True # hide cell from render

        self.unhide_everything()
        #ph.reduce_single_masks(self.filepath, [cell_info_dict["Filename"] for cell_info_dict in self.cell_info])# reduce RGBA image to only alpha channel
        
    def create_cell_info(self):    
        ''''
        TODO change to metadata
        create a list of dictionaries wich contains for each cell its type and ID 
        ''' 
        self.cell_info = []
        unique_type_counter = 0
        unique_type_dict = {}
        
        masks_path = self.filepath + f"/masks/instance_individual/{self.sample_name}/"
        metadata_path = self.filepath + f'/metadata'
        if not os.path.exists(metadata_path):
            os.makedirs(metadata_path)
            
        for idx, cell in enumerate(self.cell_objects): 
            cell_id = idx + 1
            cell_name = cell.name

            _, cell_type = hm.get_info_from_cell_name(cell_name)
            
            if cell_type not in unique_type_dict.keys(): 
                unique_type_counter +=1 
                unique_type_dict[cell_type] = unique_type_counter
                
            mask_name = f"{cell_name}.png"
            cell_filename = masks_path + '/' + mask_name
            cell_info_tuple = {"ID": cell_id, "Type": cell_type, "Filename": cell_filename, "Cellname":cell_name, "ID_Type": unique_type_dict[cell_type]}
            self.cell_info.append(cell_info_tuple)

        with open(Path(metadata_path).joinpath(f'metadata_{self.sample_name}.json'), 'w') as f:
            cell_info_dict = {i:info for i, info in enumerate(self.cell_info)}
            json.dump(cell_info_dict, f)

    def define_palette(self, type=""):
        '''
        Assign unique color to each inidividual cell or to each cell type
        '''
        if type == "semantic": 
            unique_cell_types = set([c["Type"] for c in self.cell_info]) # identify unique cell types
            print(unique_cell_types)
            cell_type_dict = {uct : (i+1) for i, uct in enumerate(unique_cell_types)} # assign unique id to each cell type
            print(cell_type_dict)
            palette = ph.make_color_palette(len(cell_type_dict.keys())) # create color palette with one color per cell type 
            print(palette)

        if type == "instance":
            self.instance_mask_names = []
            unique_cell_IDs = [c["ID"] for c in self.cell_info] # get individual cell ids 
            cell_ID_dict = {cid : (i+1) for i, cid in enumerate(unique_cell_IDs)} # assign unique id to each cell type
            palette = ph.make_color_palette(len(cell_ID_dict.keys())) # create color palette with one color per cell type 

        return palette
    
    def setup_scene_render_full_masks(self, output_shape = (500, 500), max_samples = 1024, mask_type = 'semantic'): 
        self.tissue.hide_viewport = True
        self.tissue.hide_render = True
        # hide light source
        self.light_source.light_source.hide_viewport = True
        self.light_source.light_source.hide_render = True

        self.scene.render.resolution_x = output_shape[0]
        self.scene.render.resolution_y = output_shape[1]
        self.scene.render.engine = "CYCLES"
        self.scene.cycles.samples = max_samples
                    
    def setup_node_tree_full_masks(self, filepath = None): 
                
        self.scene.render.use_compositing = True
        self.scene.use_nodes = True
        tree = self.scene.node_tree 
        nodes = tree.nodes
        links = tree.links

        for node in nodes:
            nodes.remove(node)

        render_layer_node = nodes.new('CompositorNodeRLayers')

        math_node = nodes.new("CompositorNodeMath")
        math_node.operation="DIVIDE"
        math_node.inputs[1].default_value = 255

        self.output_node = nodes.new("CompositorNodeOutputFile")
        if filepath is not None: 
            self.output_node.base_path = filepath
            if not os.path.exists(self.output_node.base_path):
                os.makedirs(self.output_node.base_path)
        else:
            self.output_node.base_path = self.filepath + f'/masks/{self.mask_type}'
        self.semantic_path = self.output_node.base_path

        self.output_node.file_slots[0].path = f"tmp_{self.sample_name}"
        self.output_node.file_slots[0].use_node_format = False
        self.output_node.file_slots[0].format.color_mode ="BW"
        self.output_node.file_slots[0].format.color_depth="16"

        links.new(render_layer_node.outputs["IndexOB"], math_node.inputs[0])       # links.new(norm_node.outputs[0], viewer_node.inputs[1])
        links.new(math_node.outputs[0], self.output_node.inputs[0]) 
     
    def set_object_pass_idx(self, flag): 

        for cell_info_dict in self.cell_info: 
            cell_object = bpy.data.objects[cell_info_dict["Cellname"]] # get cell object

            if flag == "instance": 
                cell_object.pass_index = cell_info_dict["ID"]
            if flag == "semantic":  
                cell_object.pass_index = cell_info_dict["ID_Type"]

    def export_full_mask(self, type: str = None): 
        self.hide_non_cell_objects()
        self.set_object_pass_idx(type)
        self.scene.view_layers["ViewLayer"].use_pass_object_index = True
        self.mask_type = type

        self.setup_node_tree_full_masks()
        self.scene.render.filepath = str(Path(self.semantic_path).joinpath("empty.png"))
        print(self.scene.render.filepath)

        bpy.ops.render.render('EXEC_DEFAULT', write_still=True) # render single cell mask

        # the exported image is a BW-png with the pixel values corresponding to the cell types/cell instances
        # typically most pixel values are close to 0 and not well visible in a typical image viewer. Theerfore 
        # we assign colors to to pixels to ensure visibility
        print(self.scene.render.filepath)
        palette = self.define_palette(type=type)
        path_temp = self.semantic_path + f"/tmp_{self.sample_name}0001.png"
        with Image.open(path_temp) as im:
            colored_instance_mask = np.array(im)
        if type == "semantic":
            self.semantic_ids = np.unique(im)
        if type == "instance":
            self.instance_ids = np.unique(im)
        colored_instance_mask = ph.put_palette(colored_instance_mask, palette)
        colored_instance_mask = Image.fromarray(colored_instance_mask.astype(np.uint8))
        if not os.path.exists(self.semantic_path+'_indexing'):
            os.makedirs(self.semantic_path+'_indexing')
        new_path = self.semantic_path+'_indexing'+ f"/{self.sample_name}.png"
        if os.path.exists(new_path):
            os.remove(new_path)
        Path(path_temp).rename(new_path)
        colored_instance_mask.save(str(Path(self.semantic_path).joinpath(f"{self.sample_name}.png")))
        self._clear_compositor()
    
    def enable_depth_output_render_setup(self):
        self.tissue.hide_viewport = True
        self.tissue.hide_render = True
        self.light_source.light_source.hide_viewport = True
        self.light_source.light_source.hide_render = True

        self.scene.render.resolution_x = 500
        self.scene.render.resolution_y = 500
        self.scene.render.engine = "CYCLES"
        self.scene.cycles.samples = 10
    
    def enable_depth_graph_node(self):
        self.scene.render.use_compositing = True
        self.scene.use_nodes = True
        tree = self.scene.node_tree
        self.links = tree.links
        self.scene.view_layers[0].use_pass_z = True 
        
        for n in tree.nodes:
            tree.nodes.remove(n)
        self.rl = tree.nodes.new('CompositorNodeRLayers')
        
        self.vl = tree.nodes.new('CompositorNodeViewer') # The viewer can come in handy for inspecting the results in the GUI
        self.vl.use_alpha = True        
        
        self.links.new(self.rl.outputs[0], self.vl.inputs[0])  # link Image to Viewer Image RGB
        self.links.new(self.rl.outputs['Depth'], self.vl.inputs[1])  # link Render Z to Viewer Image Alpha
    
    def enable_depth_output(self):
        bpy.app.handlers.render_complete.append(fn_print_time_when_render_done)

        self.enable_depth_output_render_setup()
        self.enable_depth_graph_node()
        
        output_path = self.filepath + f'/masks/depth'
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        self.scene.render.filepath = output_path + f'/{self.sample_name}.png'
        self.export_scene()
        
        bpy.app.handlers.render_complete.remove(fn_print_time_when_render_done)
                                
        # dmap = np.flip(np.array(pixels[:]).reshape((500, 500, 4)), axis=0)
        pixels = np.array(bpy.data.images['Viewer Node'].pixels)
        resolution = int(np.sqrt(len(pixels)/4))
        image_with_depth = pixels.reshape(resolution,resolution,4)  #reshaping into image array 4 channel (rgbz)
        
        np.save(f"{output_path}/{self.sample_name}", image_with_depth)
                
        K = cu.get_cam_intrinsic()
        new_depth = self.depth_of_point_to_depth(image_with_depth[:,:,3], K)
        LIMIT_DEPTH = 6e4
        limit_mask = new_depth < LIMIT_DEPTH
        new_depth = new_depth * limit_mask
        
        # self._clear_compositor()
    
    def depth_of_point_to_depth(self, depth_of_point, K):
        """
        # https://blender.stackexchange.com/questions/130970/cycles-generates-distorted-depth
        CYCLES's depth is distance to camera original point
        EEVEE's dpeth is XYZ's Z, same to opencv
        """
        h, w = depth_of_point.shape
        xyzs_1m = (
            np.pad(
                np.mgrid[:h, :w][::-1].reshape(2, -1), ((0, 1), (0, 0)), constant_values=1
            ).T
            @ np.linalg.inv(K).T
        )
        xyzs = (
            xyzs_1m
            * depth_of_point.flatten()[:, None]
            / np.linalg.norm(xyzs_1m, axis=1, keepdims=True)
        )
        depth = xyzs[:, 2].reshape(h, w)
        return depth

    def export_obj3d(self): 
        '''Export the entrie scene as a 3d object'''
        bpy.ops.export_scene.obj(filepath=f"{self.filepath}/my_scene.obj")

    def render(self, 
               filepath: str, 
               scene: bool = True, 
               single_masks: bool = True, 
               semantic_mask: bool = False, 
               instance_mask: bool = False,
               depth_mask: bool = False, 
               obj3d: bool = True,
               output_shape = (512, 512), 
               max_samples = 10,
               n_slices = 10, 
               slice_thickness = None):
        '''
        filepath: the folder where all outputs will be stored
        scene: if true a png of the scene will be generated
        single_masks: if true individual masks for each cell will be rendered
        semantic_mask: if true the individual masks will be combined to a single mask where pixel values distinguish between cell types
        instance_mask: if true the individual masks will be combined to a single mask each cell has a different pixel value
        depth_mask: if true a mask will be generated where pixel values correspond to depth values 
        obj3d: if true the entire scene will be exported in .OBJ format
        output_shape = width and height of output pngs
        max_samples: number of samples for rendering. Fewer samples will render more quickly
        '''
        self.filepath = filepath
        self.create_cell_info()
        
        bpy.app.handlers.render_complete.append(fn_print_time_when_render_done)

        if scene: 
            self.setup_scene_render_default(output_shape=output_shape, max_samples=max_samples)
            self.export_scene()
        
        # switch to non focus camera and switch material of nuclei
        self.camera.switch_to_mask_camera(self.scene)
        self.add_staining(bpy.data.materials.get("nuclei_mask"))

        if semantic_mask: 
            self.setup_scene_render_full_masks(output_shape=output_shape, max_samples=1)
            self.export_full_mask(type = "semantic")
            
        if instance_mask: 
            self.setup_scene_render_full_masks(output_shape=output_shape, max_samples=1)
            self.export_full_mask(type = "instance")
            
        if single_masks:# or semantic_mask or instance_mask:
            self.scan_through_tissue(
                filepath=filepath, n_slices=n_slices, slice_thickness=slice_thickness,
                output_shape=output_shape)

            # self.setup_scene_render_mask(output_shape=output_shape)
            # self.export_masks()

        if obj3d: 
            self.setup_scene_render_default(output_shape=output_shape, max_samples=max_samples)
            self.export_obj3d()

        bpy.app.handlers.render_complete.remove(fn_print_time_when_render_done)
        print("rendering completed")

    def scan_through_tissue(
            self, filepath: str, n_slices=10, slice_thickness=None, output_shape=(500, 500)): 
        '''
            This function will scan through the 3d volume and render 2d images at each slice
            which are then combined to tiff stacks.
            Args:
                filepath: the folder where all outputs will be stored
                n_slices: how many slices should be gnereated
                slice_thickness: thickness of a slice, if not provided it calcualted from n_slices 
                output_shape: tuple, resolution of output 
        '''
        self.filepath = filepath

        # fetch info from pre existing tissue
        tissue_thickness = self.tissue_empty.dimensions.z
        tissue_location = self.tissue_empty.location.z

        tissue_thickness = tissue_thickness

        # prepare slices
        if slice_thickness is None: 
            slice_thickness = tissue_thickness / n_slices
        start = tissue_location + (tissue_thickness - slice_thickness) / 2
        end = tissue_location - (tissue_thickness - slice_thickness) / 2
        slices = np.arange(start, end-end/10**4, -slice_thickness)
        self.tissue_empty.dimensions.z = slice_thickness
        print(slices)

        # setup scene for rendering
        self.hide_non_cell_objects()
        bpy.app.handlers.render_complete.append(fn_print_time_when_render_done)
        type = "instance"
        filepath = self.filepath + f'/masks/instance_3d/{self.sample_name}/'
        self.setup_scene_render_full_masks(output_shape=output_shape, max_samples=1)
        self.set_object_pass_idx(type)
        self.scene.view_layers["ViewLayer"].use_pass_object_index = True
        self.mask_type = type
        self.setup_node_tree_full_masks(filepath=filepath)
        self.scene.render.filepath = str(Path(filepath).joinpath("empty.png"))
        instance_mask_3d = []
        instance_mask_3d_indexing = []
        print(self.scene.render.filepath)
        
        for idx, loc in enumerate(slices): # TODO turn off rendering scene
            self.tissue_empty.location.z = loc
            bpy.ops.render.render('EXEC_DEFAULT', write_still=True) # render single cell mask
            print(self.scene.render.filepath)
            name = self.output_node.file_slots[0].path
            
            # colorize mask
            palette = self.define_palette(type=type)
            path_temp = filepath + f"/{name}0001.png"
            with Image.open(path_temp) as im:
                indexing = np.array(im)
                colored_instance_mask = np.copy(indexing)
            colored_instance_mask = ph.put_palette(colored_instance_mask, palette, ids=self.instance_ids)
            new_filepath = self.filepath + f'/masks/instance_3d_indexing/{self.sample_name}/'
            if not os.path.exists(new_filepath):
                os.makedirs(new_filepath)
            new_path = new_filepath + f"/slice_{self.sample_name}_{idx}.png"
            if os.path.exists(new_path):
                os.remove(new_path)
            Path(path_temp).rename(new_path)
            os.remove(self.scene.render.filepath)
            instance_mask_3d.append(colored_instance_mask)
            instance_mask_3d_indexing.append(indexing)
            colored_instance_mask = Image.fromarray(colored_instance_mask.astype(np.uint8))
            colored_instance_mask.save(str(Path(filepath).joinpath(f"slice_{self.sample_name}_{idx}.png")))
            
        # combine to numpy stack
        instance_mask_3d = np.stack(instance_mask_3d)
        instance_mask_3d_indexing = np.stack(instance_mask_3d_indexing)
        np.save(str(Path(new_filepath).joinpath(f"{self.sample_name}.npy")), instance_mask_3d_indexing)
        np.save(str(Path(filepath).joinpath(f"stack_{self.sample_name}.npy")), instance_mask_3d)

        # restore settings
        self._clear_compositor()
        self.tissue_empty.location.z = tissue_location
        self.tissue_empty.dimensions.z = tissue_thickness

    def add_dummy_objects(self, tissue, padding, vol_scale, surf_scale):
        # Create temporarily padded tissue
        tissue.tissue.scale = tuple(1+padding for _ in range(3))

        bpy.ops.mesh.primitive_cylinder_add() # Example bounding torus mesh
        cylinder = bpy.context.active_object
        #vol_obj.location = Vector(tissue.tissue.location) + Vector((0, 0, 0.5))
        cylinder.scale = vol_scale
        # NOTE: Necessary to transform the vertices of the mesh according to scale
        # It should be used when the object is created, but maybe there's a better place in the methds for it. ck
        bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
        # Intersect with tissue
        bpy.ops.mesh.primitive_cube_add(location=tissue.location)
        box = bpy.context.active_object
        box.scale = (1.1, 1.1, tissue.thickness/tissue.size)
        bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)
        vol_obj = geo.subtract_object([box], cylinder)[0]
        geo.remove_top_and_bottom_faces(vol_obj)
        geo.remove_objects([cylinder])
        vol_obj.name = "Volume"

        bpy.ops.mesh.primitive_cylinder_add()
        surf_obj = bpy.context.active_object
        #surf_obj.location = Vector(tissue.tissue.location) + Vector((0, 0, 0.5))
        surf_obj.scale = surf_scale
        # NOTE: Necessary to transform the vertices of the mesh according to scale
        # It should be used when the object is created, but maybe there's a better place in the methds for it. ck
        bpy.ops.object.transform_apply(location=True, rotation=True, scale=True) 
        # Intersect with tissue
        geo.intersect_with_object([surf_obj], tissue.tissue)
        geo.remove_top_and_bottom_faces(surf_obj)
        surf_obj.name = "Surface"
        # Rescale tissue to original scale
        tissue.tissue.scale = (1,1,1)
        return vol_obj, surf_obj

