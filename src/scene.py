import bpy
import json
import os

import src.arrangement.arrangement as arr
import src.utils.helper_methods as hm
import src.utils.plot_helpers as ph
#import imp
import importlib as imp # imp module is deprecated since python 3.12
from pathlib import Path
import time
import numpy as np
imp.reload(arr)
imp.reload(hm)
imp.reload(ph)

def fn_print_time_when_render_done(dummy):
    
    print("----- the time is: ", time.time())

class Camera:
    def __init__(self, name='camera 1', pos = (0, 0, 2), rot = (0, 0, 0), size = (2, 2)):
        self.scene = bpy.context.scene
        cam = bpy.data.cameras.new(name)
        cam.lens = 18
        
        self.cam_obj = bpy.data.objects.new(name, cam)
        self.cam_obj.location = pos
        self.cam_obj.rotation_euler = rot
        #self.scene.collection.objects.link(self.cam_obj)


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
        self.scene = bpy.context.scene
        self.scene.camera = self.camera.cam_obj
        self._clear_compositor()

    @staticmethod
    def clear():
        hm.delete_objects()

    def _clear_compositor(self):
        if self.scene.node_tree is not None:
            bpy.context.scene.node_tree.nodes.clear()

    def add_tissue(self, tissue):
        self.tissue = tissue
        bpy.context.view_layer.objects.active = bpy.context.scene.objects['tissue']
        self.tissue_empty = tissue.copy()
        self.tissue_empty.data = tissue.data.copy()
        bpy.context.collection.objects.link(self.tissue_empty)
        self.tissue_empty.name = 'tissue_empty'
        self.tissue_empty.hide_viewport = True
        self.tissue_empty.hide_render = True
        self.tissue_empty.scale.z = 1.001
    
    def cut_cells(self):
        for cell in self.cell_objects:
            boolean = cell.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
            boolean.operation = 'INTERSECT'
            boolean.object = self.tissue_empty

    def uncut_cells(self): 
        """ reverses the action of cut_cells"""
        for cell in self.cell_objects:
            mod = cell.modifiers['Boolean Modifier']
            cell.modifiers.remove(mod)

    def cut_tissue(self):
        for cell in self.cell_objects:
            boolean = self.tissue.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
            boolean.operation = 'DIFFERENCE'
            boolean.object = cell

    def add_staining(self, material):
        for cell in self.cell_objects:
            cell.data.materials.append(material)
            cell.active_material = material
    
    def add_arrangement(self, cell_arrangement: arr.CellArrangement):
        self.arrangements.append(cell_arrangement)
        cell_arrangement.add()
        print(f"Added arrangement {cell_arrangement.name} with {len(cell_arrangement.objects)} objects.")
        self.cell_objects = self.cell_objects + cell_arrangement.objects

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
        
        filepath = self.filepath + f"/train"
        if not os.path.exists(filepath):
            os.makedirs(filepath)
        
        self.scene.render.filepath = filepath + f"/{self.sample_name}/images/{self.sample_name}.png"

    def export_scene(self): 
        '''Generates a png image of the complete scene using the specifications defined in setup_scene_render_default'''
        bpy.ops.render.render('EXEC_DEFAULT', write_still=True)
    
    def export_masks(self): 
        ''''
        create a binary mask for each individuals cell object
        '''
        self.hide_everything()
        
        filepath = self.filepath + f"/train"
        if not os.path.exists(filepath):
            os.makedirs(filepath)
        self.new_filepath = filepath + f"/{self.sample_name}/masks/"
        
        for cell_info_dict in self.cell_info: 
            cell_object = bpy.data.objects[cell_info_dict["Cellname"]] # get cell object
            cell_object.hide_viewport = False # unhide cell from viewport
            cell_object.hide_render = False # unhide cell from render
            self.scene.render.filepath = self.new_filepath + cell_info_dict["Cellname"] + '.png' #self.filepath + mask_nam
            bpy.ops.render.render('EXEC_DEFAULT', write_still=True) # render single cell mask
            cell_object.hide_viewport = True # hide cell from viewport
            cell_object.hide_render = True # hide cell from render

        self.unhide_everything()
        
        ph.reduce_single_masks(self.filepath, [cell_info_dict["Filename"] for cell_info_dict in self.cell_info])# reduce RGBA image to only alpha channel
        
    def create_cell_info(self):    
        ''''
        create a list of dictionaries wich contains for each cell its type and ID 
        ''' 
        self.cell_info = []
        
        masks_path = self.filepath + f'/train/{self.sample_name}/masks'
        if not os.path.exists(masks_path):
            os.makedirs(masks_path)
            
        for idx, cell in enumerate(self.cell_objects): 
            cell_id = idx
            cell_name = cell.name 

            cell_type = hm.get_type_from_cell_name(cell_name)
            mask_name = f"{cell_name}.png"
            cell_filename = masks_path + '/' + mask_name
            cell_info_tuple = {"ID": cell_id, "Type": cell_type, "Filename": cell_filename, "Cellname":cell_name}
            self.cell_info.append(cell_info_tuple)

        with open(Path(masks_path).joinpath('data.json'), 'w') as f:
            cell_info_dict = {i:info for i, info in enumerate(self.cell_info)}
            json.dump(cell_info_dict, f)

    def define_palette(self, type=""):
        '''
        Assign unique color to each inidividual cell or to each cell type
        '''
        if type == "semantic": 
            unique_cell_types = set([c["Type"] for c in self.cell_info]) # identify unique cell types 
            cell_type_dict = {uct : (i+1) for i, uct in enumerate(unique_cell_types)} # assign unique id to each cell type
            palette = ph.make_color_palette(len(cell_type_dict.keys())) # create color palette with one color per cell type 

        if type == "instance":
            self.instance_mask_names = []
            unique_cell_IDs = [c["ID"] for c in self.cell_info] # get individual cell ids 
            cell_ID_dict = {cid : (i+1) for i, cid in enumerate(unique_cell_IDs)} # assign unique id to each cell type
            palette = ph.make_color_palette(len(cell_ID_dict.keys())) # create color palette with one color per cell type 

        return palette

    def combine_masks_semantic(self, file_name:int=0, palette=None): 
        ''' Combine Masks for individual cells into a semantic masks assinging specific colors to each cell type '''
        ph.build_semantic_mask(self.filepath, self.cell_info, file_name=file_name, palette=palette)

    def combine_masks_instance(self, file_name:int=0, palette=None): 
        ''' Combine Masks for individual cells into an instance mask assinging specific colors to each individual cell '''
        ph.build_instance_mask(self.filepath, self.cell_info, file_name=file_name, palette=palette)

    def remove_single_masks(self): 
        '''Delete png files of singel cell 2d masks'''
        cell_mask_filenames = [info_tuple["Filename"] for info_tuple in self.cell_info]
        ph.remove_single_masks(cell_mask_filenames)

    def export_depth(self): 
        pass
        # """Obtains depth map from Blender render.
        # :return: The depth map of the rendered camera view as a numpy array of size (H,W).
        # """
                
        # self.tissue.hide_viewport = True
        # self.tissue.hide_render = True
        # # hide light source
        # self.light_source.light_source.hide_viewport = True
        # self.light_source.light_source.hide_render = True

        # self.scene.render.resolution_x = 500
        # self.scene.render.resolution_y = 500
        # self.scene.render.engine = "CYCLES"
        # self.scene.cycles.samples = 100      

        # self.scene.render.use_compositing = True
        # self.scene.use_nodes = True
        # self.scene.view_layers[0].use_pass_z = True
        # tree = self.scene.node_tree 
        # nodes = tree.nodes
        # links = tree.links

        # for node in nodes:
        #     nodes.remove(node)

        # render_layer_node = nodes.new('CompositorNodeRLayers')
        # viewer_node = nodes.new('CompositorNodeViewer')
        # norm_node = nodes.new("CompositorNodeNormalize")
        # viewer_node.use_alpha = True

        # links.new(render_layer_node.outputs['Depth'], norm_node.inputs[0])  # link Render Z to Viewer Image Alpha
        # links.new(norm_node.outputs[0], viewer_node.inputs[1])

        # bpy.ops.render.render(write_still=True) 
        # pixels = bpy.data.images['Viewer Node'].pixels
        
        # dmap = np.flip(np.array(pixels[:]).reshape((500, 500, 4)), axis=0)
        # np.save(f"{self.filepath}/depth", dmap)
    
    def enable_depth_output_render_setup(self):
        self.tissue.hide_viewport = True
        self.tissue.hide_render = True
        self.light_source.light_source.hide_viewport = True
        self.light_source.light_source.hide_render = True

        self.scene.render.resolution_x = 500
        self.scene.render.resolution_y = 500
        self.scene.render.engine = "CYCLES"
        self.scene.cycles.samples = 100 
    
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
        self.enable_depth_output_render_setup()
        self.enable_depth_graph_node()
        
        output_path = self.filepath + f'/train_combined_masks/depth'
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        self.scene.render.filepath = output_path + f'/{self.sample_name}.png'
        self.export_scene()
        
        pixels = bpy.data.images['Viewer Node'].pixels
        
        # dmap = np.flip(np.array(pixels[:]).reshape((500, 500, 4)), axis=0)
        pixels = np.array(bpy.data.images['Viewer Node'].pixels)
        resolution = int(np.sqrt(len(pixels)/4))
        image_with_depth = pixels.reshape(resolution,resolution,4)  #reshaping into image array 4 channel (rgbz)
        
        np.save(f"{output_path}/{self.sample_name}", image_with_depth)

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
               output_shape = (500, 500), 
               max_samples = 10):
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
            
        if single_masks or semantic_mask or instance_mask:
            self.setup_scene_render_mask(output_shape=output_shape)
            self.export_masks()

        if semantic_mask: 
            semantic_palette = self.define_palette(type="semantic")
            self.combine_masks_semantic(palette=semantic_palette, file_name=self.sample_name)
            
            if not single_masks: 
                self.remove_single_masks()
            
        if instance_mask: 
            instance_palette = self.define_palette(type="instance")
            self.combine_masks_instance(palette=instance_palette, file_name=self.sample_name)

            if not single_masks: 
                self.remove_single_masks()
                
        if depth_mask: 
            self.enable_depth_output()
            # self.export_depth()

        # if depth_mask: 
        #     self.setup_scene_render_default(output_shape=output_shape, max_samples=max_samples)
        #     self.export_depth()

        if obj3d: 
            self.setup_scene_render_default(output_shape=output_shape, max_samples=max_samples)
            self.export_obj3d()

        bpy.app.handlers.render_complete.remove(fn_print_time_when_render_done)
        print("rendering completed")

    def hide_auxiliary_objects(self):
        for arr in self.arrangements:
            for obj in arr.auxiliary_objects:
                obj.hide_render = True
                obj.hide_viewport = True

    def render3d(self, 
               filepath: str, 
               scene: bool = False, 
               semantic_mask: bool = False, 
               instance_mask: bool = False,
               depth_mask: bool = False, 
               obj3d: bool = True,
               output_shape = (500, 500), 
               n_slices = 10, 
               slice_thickness = None, 
               min_z = 0.4, 
               max_z = 0.6, 
               semantic_mask_label="",
               instance_mask_label="", 
               max_samples = 100 ):
        
        '''
        filepath: the folder where all outputs will be stored
        scene: if true a png of the scene will be generated
        semantic_mask: if true will create a set of semantic masks where pixel values distinguish between cell types
        instance_mask: if true will create a set of semantic masks where each cell has a different pixel value
        depth_mask: if true a mask will be generated where pixel values correspond to depth values 
        obj3d: if true the entire scene will be exported in .OBJ format
        output_shape = width and height of output pngs
        n_slices: how many slices should be gnereated. will determine the slice thickness is no specific value is provided (i.e. slice_thickness=None). will be ignored if  slice_thickness!=None. 
        slice_thickness: thickness of a slice for scanning through the 3d volume. If not provided it will be calcualted from n_slices. 
        min_z, max_z: where to start/stop scnning
        semantic_mask_label = name to give to each 2d semantic mask (default is simple "semantic_mask"), masks for differenr slices will be identifies by index
        instance_mask_label = name to give to each 2d semantic mask (default is simple "instance_mask"), masks for differenr slices will be identifies by index
        max_samples: number of samples for rendering. Fewer samples will render more quickly
        '''
        # TODO automatically read min_z, max_z from scene properties 
        self.filepath = filepath

        bpy.app.handlers.render_complete.append(fn_print_time_when_render_done)

        if scene: 
            self.setup_scene_render_default(output_shape=output_shape, max_samples=max_samples)
            self.export_scene()

        if semantic_mask or instance_mask:
            self.scan_through_tissue(filepath=filepath, n_slices=n_slices, slice_thickness=slice_thickness, min_z=min_z, max_z=max_z, 
                                     output_shape=output_shape, semantic_mask=semantic_mask, semantic_mask_label=semantic_mask_label, instance_mask=instance_mask, instance_mask_label=instance_mask_label)

        if depth_mask: 
            self.export_depth()

        if obj3d: 
            self.export_obj3d()

        bpy.app.handlers.render_complete.remove(fn_print_time_when_render_done)
        print("rendering completed")


    def scan_through_tissue(self, filepath: str, n_slices=10, slice_thickness=None, min_z = 0.4, max_z =0.6, 
                            output_shape=(500, 500), semantic_mask=True, semantic_mask_label="", instance_mask=True, instance_mask_label=""): 
        '''
            filepath: the folder where all outputs will be stored
            n_slices: how many slices should be gnereated. will determine the slice thickness is no specific value is provided (i.e. slice_thickness=None). will be ignored if  slice_thickness!=None. 
            slice_thickness: thickness of a slice for scanning through the 3d volume. If not provided it will be calcualted from n_slices. 
            min_z, max_z: where to start/stop scnning
            output_shape: tuple, resolution of output 
            semantic_mask: if true will create a set of semantic masks where pixel values distinguish between cell types
            semantic_mask_label = name to give to each 2d semantic mask (default is simple "semantic_mask"), masks for differenr slices will be identifies by index
            instance_mask: if true will create a set of semantic masks where each cell has a different pixel value
            instance_mask_label = name to give to each 2d semantic mask (default is simple "instance_mask"), masks for differenr slices will be identifies by index
        '''
         
        self.filepath = filepath

        if slice_thickness is None: 
            # calculate slice thickness from desired number of slices 
            #TODO throw error when neither thickness nor nslices are given or they are contradictory
            slices = np.linspace(min_z, max_z, n_slices)
            slice_thickness = slices[1] - slices[0]
        else: 
            # calculate number of slices from slice thickness 
            slices = np.arange(min_z, max_z, slice_thickness)
            n_slices = len(slices)

        # scale the moving slice to desired thickness 
        self.tissue_empty.scale.z = slice_thickness 

        self.create_cell_info()

        if semantic_mask: 
            self.semantic_mask_names = []
            semantic_palette = self.define_palette(type="semantic")
        if instance_mask: 
            self.instance_mask_names = []
            instance_palette = self.define_palette(type="instance")

        
        bpy.app.handlers.render_complete.append(fn_print_time_when_render_done)

        for idx, loc in enumerate(slices): # move through scene 
            
            self.tissue_empty.location.z = loc # position moving slice 
            self.cut_cells()# make cell(-parts) invisible outside moving tissue
            self.setup_scene_render_mask(output_shape=output_shape) # render 2d scene at the intersection of scene and moving tissue
            self.export_masks() # export individual cell masks 

            if semantic_mask: 
                # combine individual cells to a 2d semnatic mask
                semantic_mask_name = f"{semantic_mask_label if len(semantic_mask_label)!=0 else 'semantic_mask'}_{idx}"
                self.combine_masks_semantic(file_name=semantic_mask_name, palette=semantic_palette)
                self.semantic_mask_names.append(Path(self.filepath).joinpath(f"{semantic_mask_name}.png"))

            if instance_mask: 
                # combine individual cells to a 2d instance mask
                instance_mask_name = f"{instance_mask_label if len(instance_mask_label)!=0 else 'instance_mask'}_{idx}"
                self.combine_masks_instance(file_name = instance_mask_name, palette=instance_palette)
                self.instance_mask_names.append(Path(self.filepath).joinpath(f"{instance_mask_name}.png"))

            # delete 2d masks of singel cells 
            self.remove_single_masks()

        bpy.app.handlers.render_complete.remove(fn_print_time_when_render_done)
        print("mask rendering completed")
        
        print("combining masks to gif")
        ph.build_gif(self.semantic_mask_names, Path(self.filepath).joinpath("semantic_mask.gif"))
        ph.build_gif(self.instance_mask_names, Path(self.filepath).joinpath("instance_mask.gif"))
        print("done combining masks to gif")
