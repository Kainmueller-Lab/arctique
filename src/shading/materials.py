import bpy
import os
import sys
import importlib as imp
import numpy as np

dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir)

import src.shading.utils as shading_utils
import src.shading.shaders as shaders
imp.reload(shading_utils)
imp.reload(shaders)



class Material():
    def __init__(self, seed=0, cell_type_params=None, tissue_rips=-0.5, tissue_rips_std=0.1, stroma_intensity=1, stroma_color=(0.62, 0.25, 0.65, 1.0)):
        # delete all materials
        for material in bpy.data.materials:
            bpy.data.materials.remove(material)
        
        # add custom nodes
        np.random.seed(seed)
        self.shift = tuple(np.random.randint([10**3]*3))
        self.custom_nodes = shaders.CustomShaderNodes(shift=self.shift, stroma_intensity=stroma_intensity)

        # add materials
        self.light_source = self.add_light_source()
        rips = np.random.normal(tissue_rips, tissue_rips_std)
        self.muscosa = self.add_mucosa_staining(threshold_rips=-rips, stroma_intensity=stroma_intensity, base_color=stroma_color)
        self.nuclei_mask = self.add_nuclei_mask()
        self.crypt_staining = self.add_crypt_staining(staining_intensity=270*stroma_intensity)
        self.cell_staining = []
        if cell_type_params is None:
            self.nuclei_staining = self.add_nuclei_staining(name="Nucleus")
            self.cytoplasm_staining = self.add_nuclei_staining(name="Cytoplasm", color=(0.605, 0.017, 0.043, 1),
                staining_intensity=50)
        else:
            for cell_type, parts in cell_type_params.items():
                for cell_part, params in parts.items():
                    print(cell_type, cell_part, params)
                    self.cell_staining.append(self.add_nuclei_staining(**params))
                    

    def add_light_source(self, brightness=60, name='light_source'):
        # add new material and node tree
        material = bpy.data.materials.new(name=name)
        material.use_nodes = True
        nodes = material.node_tree.nodes

        # add emission node
        principled_bsdf = nodes.get("Principled BSDF")
        principled_bsdf.inputs['Emission Strength'].default_value = brightness
        principled_bsdf.inputs['Emission'].default_value = (1.0, 1.0, 1.0, 1.0)  # RGB color

        # link nodes
        material_output = material.node_tree.nodes.get('Material Output')
        material.node_tree.links.new(principled_bsdf.outputs['BSDF'], material_output.inputs['Surface'])

        return material

    def add_nuclei_mask(self, transmission=0.5, name="nuclei_mask", base_color=(0.05, 0.0, 0.3, 1.0)):
        # add new material and node tree
        material = bpy.data.materials.new(name=name)
        material.use_nodes = True
        nodes = material.node_tree.nodes

        # add emission node
        principled_bsdf = nodes.get("Principled BSDF")
        principled_bsdf.inputs['Transmission'].default_value = transmission
        principled_bsdf.inputs['Transmission Roughness'].default_value = 1.0
        principled_bsdf.inputs['Base Color'].default_value = base_color  # RGB color

        # link nodes
        material_output = material.node_tree.nodes.get('Material Output')
        material.node_tree.links.new(principled_bsdf.outputs['BSDF'], material_output.inputs['Surface'])

        return material

    def add_crypt_staining(
            self, name="crypt", color=(0.456, 0.011, 0.356, 1),
            staining_intensity=220, start_pos=(0, 0), sep=200):
        
        material, nodes, links = shading_utils.initialize_material(name)
        
        ###### VOLUME
        # add object centered coordinate system
        node = coord = shading_utils.add_node_group(
            nodes, self.custom_nodes.object_coord, pos=(0, 0))
        loc = node.location

        # add noisy intensity
        node = noise = nodes.new('ShaderNodeTexNoise')
        node.inputs['Scale'].default_value = 100
        node.inputs['Detail'].default_value = 2
        node.inputs['Roughness'].default_value = 0.5
        node.inputs['Distortion'].default_value = 0
        loc = node.location = (loc[0]+sep, loc[1])
        node = intensity = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc[0]+sep, loc[1])
        node.color_ramp.elements[0].color = (0.8, 0.8, 0.8, 1)
        node.color_ramp.elements[1].color = (1, 1, 1, 1)
        node.color_ramp.elements[0].position = 0.349
        node.color_ramp.elements[1].position = 1
        links.new(coord.outputs[0], noise.inputs['Vector'])
        links.new(noise.outputs[0], intensity.inputs[0])

        # multiply density
        node = density = nodes.new('ShaderNodeMath')
        node.operation = 'MULTIPLY'
        node.inputs[1].default_value = staining_intensity
        loc = node.location = (loc[0]+sep, loc[1])
        links.new(intensity.outputs[0], density.inputs[0])

        # add volume shader
        node = volume = shading_utils.add_node_group(
            nodes, self.custom_nodes.volume, pos=(loc[0]+sep, loc[1]))
        loc = node.location
        links.new(density.outputs[0], volume.inputs['AbsorptionDensity'])
        volume.inputs['AbsorptionColor'].default_value = color
        volume.inputs['AbsorptionDensity'].default_value = 300
        volume.inputs['ScatterDensity'].default_value = 0

        # ADD TISSUE NOISE
        node = tissue_noise = shading_utils.add_node_group(
            nodes, self.custom_nodes.stacked_noise, pos=(loc[0]+sep, loc[1]))
        loc = node.location
        tissue_noise.inputs['StrengthRips'].default_value = 0
        links.new(volume.outputs[0], tissue_noise.inputs['Shader'])

        ###### SURFACE -> TODO add later again
        # node = surface = nodes.new('ShaderNodeBsdfGlass')
        # loc = node.location = (loc[0], loc[1]+sep)

        # link nodes
        node = material_output = nodes.new('ShaderNodeOutputMaterial')
        links.new(tissue_noise.outputs[0], material_output.inputs['Volume'])
        #links.new(surface.outputs[0], material_output.inputs['Surface'])
        loc = node.location = (loc[0]+sep, loc[1])

        return material
    
    def add_nuclei_staining(
            self, name="Nucleus", color=(0.315, 0.003, 0.631, 1),
            staining_intensity=200, start_pos=(0, 0), sep=200):
        material, nodes, links = shading_utils.initialize_material(name)

        # add object centered coordinate system
        node = coord = shading_utils.add_node_group(
            nodes, self.custom_nodes.object_coord, pos=start_pos)
        loc = node.location

        # add noisy intensity
        node = noise = nodes.new('ShaderNodeTexNoise')
        node.inputs['Scale'].default_value = 12.2
        node.inputs['Detail'].default_value = 2
        node.inputs['Roughness'].default_value = 0.5
        node.inputs['Distortion'].default_value = 0
        loc = node.location = (loc[0]+sep, loc[1])
        node = intensity = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc[0]+sep, loc[1])
        node.color_ramp.elements[0].color = (0.468, 0.468, 0.468, 1)
        node.color_ramp.elements[1].color = (1, 1, 1, 1)
        node.color_ramp.elements[0].position = 0.349
        node.color_ramp.elements[1].position = 0.636
        links.new(coord.outputs[0], noise.inputs['Vector'])
        links.new(noise.outputs[0], intensity.inputs[0])

        # multiply density
        node = density = nodes.new('ShaderNodeMath')
        node.operation = 'MULTIPLY'
        node.inputs[1].default_value = staining_intensity
        loc = node.location = (loc[0]+sep, loc[1])
        links.new(intensity.outputs[0], density.inputs[0])

        # add volume shader
        node = volume = shading_utils.add_node_group(
            nodes, self.custom_nodes.volume, pos=(loc[0]+sep, loc[1]))
        loc = node.location
        links.new(density.outputs[0], volume.inputs['AbsorptionDensity'])
        volume.inputs['AbsorptionColor'].default_value = color
        volume.inputs['ScatterDensity'].default_value = 0

        # add principal noise
        node = noise = shading_utils.add_node_group(
            nodes, self.custom_nodes.principle_noise, pos=(loc[0]+sep, loc[1]))
        loc = node.location
        node.inputs['Size'].default_value = 60
        links.new(volume.outputs[0], noise.inputs[0])

        # link nodes
        node = material_output = nodes.new('ShaderNodeOutputMaterial')
        links.new(noise.outputs[0], material_output.inputs['Volume'])
        loc = node.location = (loc[0]+sep, loc[1])

        return material
    

    def add_mucosa_staining(
            self, name="muscosa", base_color=(0.62, 0.25, 0.65, 1.0),
            start_pos=(0, 0), sep=200, threshold_rips=0.5, stroma_intensity=1):
        
        material, nodes, links = shading_utils.initialize_material(name)
        
        ###### VOLUME
        # object centric coordinate system
        node = coord = shading_utils.add_node_group(
            nodes, self.custom_nodes.object_coord, pos=start_pos)
        loc = node.location
        
        # statining intensity
        node = staining_noise = nodes.new('ShaderNodeTexNoise')
        node.inputs['Scale'].default_value = 2
        node.inputs['Detail'].default_value = 0
        node.inputs['Roughness'].default_value = 0
        node.inputs['Distortion'].default_value = 1.1
        loc = node.location = (loc[0]+sep, loc[1])
        node = staining_intensity = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc[0]+sep, loc[1])
        node.color_ramp.elements[0].color = (1, 1, 1, 1)
        node.color_ramp.elements[1].color = (0.475, 0.475, 0.475, 1)
        node.color_ramp.elements[0].position = 0.1
        node.color_ramp.elements[1].position = 1
        links.new(coord.outputs[0], staining_noise.inputs['Vector'])
        links.new(staining_noise.outputs[0], staining_intensity.inputs[0])
        
        # BACKGROUND
        node = lamina_propria_tissue_base = shading_utils.add_node_group(
            nodes, self.custom_nodes.lamina_propria_tissue_base,
            pos=(loc[0], loc[1]-1.5*sep))
        loc = node.location
        node = staining_base = nodes.new('ShaderNodeMixShader')
        loc = node.location = (loc[0]+1.5*sep, loc[1]+1.5*sep)
        links.new(lamina_propria_tissue_base.outputs[0], staining_base.inputs[2])
        links.new(staining_intensity.outputs[0], staining_base.inputs['Fac'])

        # ADD TISSUE NOISE
        node = tissue_noise = shading_utils.add_node_group(
            nodes, self.custom_nodes.stacked_noise, pos=(loc[0]+sep, loc[1]))
        loc = node.location
        tissue_noise.inputs['ThresholdRips'].default_value = threshold_rips
        links.new(staining_base.outputs[0], tissue_noise.inputs['Shader'])
        
        # RED POINTS
        node = red_points = shading_utils.add_node_group(
            nodes, self.custom_nodes.mixing_red, pos=(loc[0]+sep, loc[1]))
        loc = node.location
        links.new(tissue_noise.outputs[0], red_points.inputs['Shader'])
        
        ###### SURFACE -> TODO add later again
        # node = surface = nodes.new('ShaderNodeBsdfGlass')
        # loc = node.location = (loc[0], loc[1]+sep)
        
        ###### OUTPUT
        node = material_output = nodes.new('ShaderNodeOutputMaterial')
        loc = node.location = (loc[0]+sep, loc[1])
        links.new(red_points.outputs[0], material_output.inputs['Volume'])
        # links.new(surface.outputs[0], material_output.inputs['Surface'])
        
        return material
