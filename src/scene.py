import bpy
import src.arrangement.arrangement as arr
import src.utils.helper_methods as hm
import imp
imp.reload(arr)
imp.reload(hm)

class Camera:
    def __init__(self, name='camera 1', pos = (0, 0, 2), rot = (0, 0, 0), size = (2, 2)):
        self.scene = bpy.context.scene
        cam = bpy.data.cameras.new(name)
        cam.lens = 18

        # create camera object
        cam_obj = bpy.data.objects.new(name, cam)
        cam_obj.location = pos
        cam_obj.rotation_euler = rot
        self.scene.collection.objects.link(cam_obj)

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
        self.scene = bpy.context.scene#camera.scene 
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

    def add_pass_index(self):
        for cell in self.cell_objects:
            pass
    
    def render_engine(self):
        self.scene.render.engine = 'CYCLES' #or 'BLENDER_EEVEE'
        self.scene.render.filepath = self.filepath + "cells.png" #render full picture
        self.scene.render.image_settings.file_format = 'PNG'
    
    def add_output_file(self):
        self.output_file = self.tree.nodes.new("CompositorNodeOutputFile")
        self.output_file.base_path = self.filepath
        self.output_file.format.file_format = 'PNG'
    
    def mask_objects(self, index):
        maskid = self.tree.nodes.new('CompositorNodeIDMask')
        maskid.index = index
        self.scene.node_tree.links.new(self.render_nodes.outputs['IndexOB'], maskid.inputs[0])
        self.output_file.file_slots.new(f"mask_output_{index}")
        self.scene.node_tree.links.new(maskid.outputs['Alpha'], self.output_file.inputs[1])
    
    def render(self, filepath: str, mask: str = False, index: int = 0):
        self.filepath = filepath
        self.render_engine()
        
        self.scene.use_nodes = True
        self.tree = self.scene.node_tree
        self.render_nodes = bpy.context.scene.node_tree.nodes.new(type='CompositorNodeRLayers')
        bpy.context.scene.view_layers["ViewLayer"].use_pass_object_index = True
        self.add_output_file()
        self.scene.node_tree.links.new(self.render_nodes.outputs['Image'], self.output_file.inputs[0])
        
        bpy.ops.render.render('INVOKE_DEFAULT', write_still=True)
        
        if mask:
            self.mask_objects(index)
            bpy.ops.render.render('INVOKE_DEFAULT', write_still=True)


        
        
