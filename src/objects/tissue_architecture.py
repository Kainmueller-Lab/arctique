import bpy
import src.objects.macro_structures as macro
import numpy as np


class TissueArch():
    def __init__(self, seed=0):
        # 1) build epithelial layer
        # 1a) build crypt
        self.seed = seed
        np.random.seed(self.seed)
        self.crypt = macro.build_crypt()
        # 1b) build surronding muscosa layer
        self.tissue = macro.build_muscosa(self.crypt)

    def get_architecture(self):
        return [
            self.crypt.crypt, self.crypt.crypt_vol_in, 
            self.crypt.crypt_vol_out, self.crypt.crypt_goblet, self.tissue.muscosa]
    
    def random_crop(self, crop):
        objects = self.get_architecture()
        #objects = [crypt, muscosa]
        self.random_translate(objects)
        self.random_rotate(objects)
        for obj in objects:
            bpy.context.view_layer.objects.active = obj
            obj.select_set(True)
            bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)
            
    def random_rotate(self, objects):
        rotation = np.random.uniform([0, 0, 0], [np.pi/2, 0, 2*np.pi])
        for obj in objects:
            # Apply the rotation around the Z axis
            obj.rotation_mode = 'XYZ'
            for i in range(3):
                obj.rotation_euler[i] += rotation[i]
    
    def random_translate(self, objects):
        translation = np.random.uniform(0, 3.5)
        for obj in objects:
            obj.location[2] -= translation
            bpy.context.scene.cursor.location = (0, 0, 0)
            bpy.context.view_layer.objects.active = obj
            obj.select_set(True)
            bpy.ops.object.origin_set(type='ORIGIN_CURSOR', center='MEDIAN')
    

