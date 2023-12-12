import bpy

class Material():
    def __init__(self):
        # delete all materials
        for material in bpy.data.materials:
            bpy.data.materials.remove(material)
        self.light_source = self.add_light_source()

    def add_light_source(self, brightness=60):
        # add new material and node tree
        material = bpy.data.materials.new(name="light_source")
        material.use_nodes = True
        nodes = material.node_tree.nodes

        # add emission node
        principled_bsdf = nodes.get("Principled BSDF")
        principled_bsdf.inputs['Emission Strength'].default_value = 60.0
        principled_bsdf.inputs['Emission'].default_value = (1.0, 1.0, 1.0, 1.0)  # RGB color

        # link nodes
        material_output = material.node_tree.nodes.get('Material Output')
        material.node_tree.links.new(principled_bsdf.outputs['BSDF'], material_output.inputs['Surface'])

        return material

    def add_nuclei_staining(self, name):
        pass

    def add_tissue_staining():
        pass