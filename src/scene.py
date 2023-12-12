import bpy
import src.utils.helper_methods as hm
import imp
imp.reload(hm)

class Camera:
    def __init__(self, name='camera 1', pos = (0, 0, 2), rot = (0, 0, 0), size = (2, 2)):
        scn = bpy.context.scene
        cam = bpy.data.cameras.new(name)
        cam.lens = 18

        # create camera object
        cam_obj = bpy.data.objects.new(name, cam)
        cam_obj.location = pos
        cam_obj.rotation_euler = rot
        scn.collection.objects.link(cam_obj)


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

    @staticmethod
    def clear():
        hm.delete_objects()
    
    def add_arangement(self):
        pass
    
    def render(self):
        pass
