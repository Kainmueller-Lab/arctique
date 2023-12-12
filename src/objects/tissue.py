import bpy
from src.shading.shading import Material


class Tissue():
    def __init__(self, material: Material, thickness=0.1, size=2, location=(0, 0, 0.5)):
        self.material = material

        # create tissue mesh
        bpy.ops.mesh.primitive_cube_add(size=2, location=location)
        self.tissue = bpy.context.active_object

        # resize to desired scale
        self.tissue.scale = (1, 1, thickness/size)
        bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)

        # add shading
        bpy.context.active_object.data.materials.append(material)
        bpy.context.active_object.active_material = material