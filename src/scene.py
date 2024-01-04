import bpy
import src.arrangement.arrangement as arr
import src.utils.helper_methods as hm
import src.utils.plot_helpers as ph
#import imp
import importlib as imp # imp module is deprecated since python 3.12
import time
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
    def __init__(self, light_source: LightSource, camera: Camera):
        self.light_source = light_source
        self.camera = camera
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
            boolean = cell.cell_object.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
            boolean.operation = 'INTERSECT'
            boolean.object = self.tissue_empty

    def uncut_cells(self): 
        """ reverses the action of cut_cells"""
        for cell in self.cell_objects:
            mod = cell.cell_object.modifiers['Boolean Modifier']
            cell.cell_object.modifiers.remove(mod)

    def cut_tissue(self):
        for cell in self.cell_objects:
            boolean = self.tissue.modifiers.new(name="Boolean Modifier", type='BOOLEAN')
            boolean.operation = 'DIFFERENCE'
            boolean.object = cell.cell_object

    def add_staining(self, material):
        for cell in self.cell_objects:
            cell.cell_object.data.materials.append(material)
            cell.cell_object.active_material = material
    
    def add_arangement(self, cell_arrangement: arr.CellArrangement):
        cell_arrangement.generate_cells()
        self.cell_objects = self.cell_objects + cell_arrangement.objects
        cell_arrangement.add()
        self.objects_list = cell_arrangement.objects

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
            cell.cell_object.hide_viewport = True
            cell.cell_object.hide_render = True

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
            cell.cell_object.hide_viewport = False
            cell.cell_object.hide_render = False

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
        self.scene.render.filepath = self.filepath + "scene.png"

    def export_scene(self): 
        '''Generates a png image of the complete scene using the specifications defined in setup_scene_render_default'''
        bpy.ops.render.render('EXEC_DEFAULT', write_still=True)
    
    def export_masks(self): 

        self.hide_everything()
        self.cell_info = []
        for cell in self.cell_objects: 
            # show only one cell 
            cell.cell_object.hide_viewport = False
            cell.cell_object.hide_render = False
            # save mask for one cell 
            mask_name = f"{cell.cell_name}.png"
            self.scene.render.filepath = self.filepath + mask_name
            
            cell_id = cell.cell_id   # cell.cell_object.cell_id
            cell_type = cell.cell_attributes.cell_type
            cell_filename = self.filepath + mask_name
            cell_info_tuple = (cell_id, cell_type, cell_filename)
            self.cell_info.append(cell_info_tuple)

            bpy.ops.render.render('EXEC_DEFAULT', write_still=True)
            # hide cell again
            cell.cell_object.hide_viewport = True
            cell.cell_object.hide_render = True

        self.unhide_everything()
        ph.reduce_single_masks(self.filepath, [info_tuple[2] for info_tuple in self.cell_info])


    def combine_masks_semantic(self, file_name="semantic_mask"): 
        ph.build_semantic_mask(self.filepath, self.cell_info, file_name=file_name)

    def combine_masks_instance(self): 
        cell_mask_filenames = [info_tuple[2] for info_tuple in self.cell_info]
        
        ph.build_instance_mask(self.filepath, cell_mask_filenames)

    def remove_single_masks(self): 
        cell_mask_filenames = [info_tuple[2] for info_tuple in self.cell_info]
        ph.remove_single_masks(cell_mask_filenames)

    def export_depth(self): 
        pass

    def export_obj3d(self): 
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
        max_samples: number of samples for rendering. Fewer samples will render more quickly
        '''

        self.filepath = filepath

        bpy.app.handlers.render_complete.append(fn_print_time_when_render_done)

        if scene: 
            self.setup_scene_render_default(output_shape=output_shape, max_samples=max_samples)
            self.export_scene()

        if single_masks or semantic_mask or instance_mask:
            self.setup_scene_render_mask(output_shape=output_shape)
            self.export_masks()

            if semantic_mask: 
                self.combine_masks_semantic()
            if instance_mask: 
                self.combine_masks_instance()

            if not single_masks: 
                self.remove_single_masks()

        if depth_mask: 
            self.export_depth()

        if obj3d: 
            self.export_obj3d()

        bpy.app.handlers.render_complete.remove(fn_print_time_when_render_done)
        print("rendering completed")



    def render3d(self, 
               filepath: str, 
               scene: bool = False, 
               semantic_mask: bool = False, 
               instance_mask: bool = False,
               depth_mask: bool = False, 
               obj3d: bool = True,
               output_shape = (500, 500), 
               semantic_mask_name = "Semantic_mask_1.png"):


        self.filepath = filepath

        bpy.app.handlers.render_complete.append(fn_print_time_when_render_done)

        if scene: 
            pass

        if semantic_mask or instance_mask:
            self.setup_scene_render_mask(output_shape=output_shape)
            self.export_masks()

            if semantic_mask: 
                self.combine_masks_semantic(filename = semantic_mask_name)
            if instance_mask: 
                self.combine_masks_instance()

            self.remove_single_masks()

        if depth_mask: 
            self.export_depth()

        if obj3d: 
            self.export_obj3d()

        bpy.app.handlers.render_complete.remove(fn_print_time_when_render_done)
        print("rendering completed")



        
        
