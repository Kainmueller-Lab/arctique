import bpy
import os
import sys
import importlib as imp

dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir)

import src.shading.utils as shading_utils
imp.reload(shading_utils)



class CustomShaderNodes():
    def __init__(self, start=(0, 0), sep=200, shift=(0, 0, 0), stroma_intensity=1):
        self.shift = shift
        self.start = start
        self.sep = sep
        
        # general
        self.object_coord = self.add_object_coord(shift=self.shift)
        self.volume = self.add_volume()
        self.mixing_red = self.add_mixing_red_points()
        self.principle_noise = self.add_principle_noise()
        self.stacked_noise = self.add_stacked_noise()
        
        # tissue
        self.lamina_propria_tissue_base = self.add_lamina_propria_tissue_base(stroma_intensity=stroma_intensity)
    
    def add_object_coord(self, node_name='Coord', shift=(0, 0, 0)):
        '''
        creates a node for fetching a object centric coordinate system
        '''
        # create nodes
        node_group, inputs, outputs = shading_utils.create_node_group(node_name, start=self.start)
        nodes = node_group.nodes
        links = node_group.links
        node_group.outputs.new('NodeSocketVector', 'Vector')
        node_group.inputs.new('NodeSocketVector', 'Shift')
        node_group.inputs['Shift'].default_value = shift
        
        # adds object origin
        node = coord = nodes.new('ShaderNodeTexCoord')
        loc = node.location = (self.start[0]+self.sep, self.start[1])
        node = coord_scale = nodes.new('ShaderNodeVectorMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        coord_scale.operation = 'SCALE'
        node = coord_shift = nodes.new('ShaderNodeVectorMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        coord_shift.operation = 'ADD'
        links.new(coord.outputs['Object'], coord_scale.inputs['Vector'])
        links.new(inputs.outputs[0], coord_shift.inputs[1])
        links.new(coord_scale.outputs[0], coord_shift.inputs[0])
        
        # connect to output
        links.new(coord_shift.outputs[0], outputs.inputs[0])
        
        return node_group
    
    def add_volume(self, node_name='Volume'):
        # create nodes
        node_group, inputs, outputs = shading_utils.create_node_group(node_name, start=self.start)
        nodes = node_group.nodes
        links = node_group.links
        node_group.outputs.new('NodeSocketShader', 'Shader')
        node_group.inputs.new('NodeSocketColor', 'AbsorptionColor')
        node_group.inputs.new('NodeSocketFloat', 'AbsorptionDensity')
        node_group.inputs.new('NodeSocketColor', 'ScatterColor')
        node_group.inputs['ScatterColor'].default_value = (1, 1, 1, 1)
        node_group.inputs.new('NodeSocketFloat', 'ScatterDensity')
        node_group.inputs['ScatterDensity'].default_value = 1
        
        # add two volume shaders
        node = abs = nodes.new('ShaderNodeVolumeAbsorption')
        loc = node.location = (self.start[0]+self.sep, self.start[1])
        node = scatter = nodes.new('ShaderNodeVolumeScatter')
        loc = node.location = (loc[0], loc[1]-self.sep)
        node = add = nodes.new('ShaderNodeAddShader')
        loc = node.location = (loc[0]+self.sep, loc[1]+self.sep)
        links.new(abs.outputs['Volume'], add.inputs[0])
        links.new(scatter.outputs['Volume'], add.inputs[1])
        
        # connect to inputs and outputs
        links.new(add.outputs[0], outputs.inputs[0])
        links.new(inputs.outputs[0], abs.inputs[0])
        links.new(inputs.outputs[1], abs.inputs[1])
        links.new(inputs.outputs[2], scatter.inputs[0])
        links.new(inputs.outputs[3], scatter.inputs[1])
        
        return node_group
    
    def add_lamina_propria_tissue_base(self, node_name='LaminaPropBase', stroma_intensity=1, stroma_color=(0.5, 0.5, 0.5)):
        node_group, inputs, outputs = shading_utils.create_node_group(node_name, start=self.start)
        nodes = node_group.nodes
        links = node_group.links
        node_group.outputs.new('NodeSocketShader', 'Shader')
        
        # object centric coordinate system
        node = coord = shading_utils.add_node_group(nodes, self.object_coord, pos=self.start)
        loc = node.location
            
        # texture
        node = musgrave_noise = nodes.new('ShaderNodeTexMusgrave')
        loc_tex = node.location = (loc[0]+self.sep, loc[1])
        links.new(coord.outputs['Vector'], node.inputs['Vector'])
        node.inputs['Scale'].default_value = 16
        
        # color
        node = color_ramp_1 = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc_tex[0]+self.sep, loc_tex[1])
        node.color_ramp.elements[0].color = (1, 1, 1, 1)
        node.color_ramp.elements[1].color = (0, 0, 0, 1)
        node.color_ramp.elements[0].position = 0.145
        node.color_ramp.elements[1].position = 0.190
        node = color_ramp_2 = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc[0]+1.5*self.sep, loc[1])
        node.color_ramp.elements[0].color = (0.07, 0.04, 0.08, 1)
        node.color_ramp.elements[1].color = (0.456, 0.011, 0.356, 1)
        node.color_ramp.elements[1].position = 0.1
        links.new(musgrave_noise.outputs[0], color_ramp_1.inputs[0])
        links.new(color_ramp_1.outputs[0], color_ramp_2.inputs[0])
        
        # density 
        node = density_ramp = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc_tex[0]+self.sep, loc_tex[1]-2*self.sep)
        node.color_ramp.elements[0].color = (1, 1, 1, 1)
        node.color_ramp.elements[1].color = (0, 0, 0, 1)
        node.color_ramp.elements[0].position = 0.4
        node.color_ramp.elements[1].position = 0.48
        node = density_strength = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+1.5*self.sep, loc[1])
        node.operation = 'MULTIPLY'
        node.inputs[1].default_value = 180 * stroma_intensity
        links.new(musgrave_noise.outputs[0], density_ramp.inputs[0])
        links.new(density_ramp.outputs[0], density_strength.inputs[0])
        
        # volume shader
        node = vol = shading_utils.add_node_group(
            nodes, self.volume,
            pos=(color_ramp_2.location[0]+1.5*self.sep, color_ramp_2.location[1]))
        loc = node.location
        links.new(density_strength.outputs[0], vol.inputs['AbsorptionDensity'])
        links.new(color_ramp_2.outputs[0], vol.inputs['AbsorptionColor'])
        
        # connect to inputs and outputs
        links.new(vol.outputs[0], outputs.inputs[0])
        
        return node_group
    
    def add_mixing_red_points(self, node_name='MixingRedPoints'):
        node_group, inputs, outputs = shading_utils.create_node_group(node_name, start=self.start)
        nodes = node_group.nodes
        links = node_group.links
        node_group.outputs.new('NodeSocketShader', 'Shader')
        node_group.inputs.new('NodeSocketShader', 'Shader')       
        
        # volum shader
        node = vol = shading_utils.add_node_group(
            nodes, self.volume,
            pos=(self.start[0]+self.sep, self.start[1]-1.5*self.sep))
        loc = node.location
        node.inputs['AbsorptionColor'].default_value = (0.605, 0.017, 0.043, 1)
        node.inputs['ScatterColor'].default_value = (0.605, 0.019, 0.088, 1)
        node.inputs['AbsorptionDensity'].default_value = 125
        node.inputs['ScatterDensity'].default_value = 0.6
        
        # object centric coordinate system
        node = coord = shading_utils.add_node_group(
            nodes, self.object_coord,
            pos=(self.start[0]+self.sep, self.start[1]))
        loc = node.location
        
        # statining intensity
        node = staining_noise = nodes.new('ShaderNodeTexMusgrave')
        node.inputs['Scale'].default_value = 5
        node.inputs['Detail'].default_value = 2
        loc = node.location = (loc[0]+self.sep, loc[1])
        node = staining_intensity = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.color_ramp.elements[0].position = 0.559
        node.color_ramp.elements[1].position = 0.75
        links.new(coord.outputs[0], staining_noise.inputs['Vector'])
        links.new(staining_noise.outputs[0], staining_intensity.inputs[0])
        
        # Mix Shader
        node = mix_shader = nodes.new('ShaderNodeMixShader')
        loc = node.location = (loc[0]+1.5*self.sep, loc[1])
        links.new(vol.outputs[0], mix_shader.inputs[2])
        links.new(staining_intensity.outputs[0], mix_shader.inputs['Fac'])
        
        # connect to inputs and outputs
        links.new(inputs.outputs[0], mix_shader.inputs[1])
        outputs.location = (loc[0]+self.sep, loc[1])
        links.new(mix_shader.outputs[0], outputs.inputs[0])
        
        return node_group


    def add_principle_noise(self, node_name='PrincipleNoise'):
        node_group, inputs, outputs = shading_utils.create_node_group(node_name, start=self.start)
        nodes = node_group.nodes
        links = node_group.links
        node_group.inputs.new('NodeSocketShader', 'Shader')
        node_group.inputs.new('NodeSocketFloat', 'Size')
        node_group.inputs['Size'].default_value = 40
        node_group.inputs.new('NodeSocketFloat', 'Threshold')
        node_group.inputs['Threshold'].default_value = 0.36
        node_group.inputs.new('NodeSocketFloat', 'Softness')
        node_group.inputs['Softness'].default_value = 0.24
        node_group.inputs.new('NodeSocketFloat', 'Dimension')
        node_group.inputs['Dimension'].default_value = 0
        node_group.inputs.new('NodeSocketFloat', 'Lacunarity')
        node_group.inputs['Lacunarity'].default_value = 1.7
        node_group.inputs.new('NodeSocketFloat', 'Strength')
        node_group.inputs['Strength'].default_value = 0.6
        node_group.inputs.new('NodeSocketFloat', 'CorrStrength')
        node_group.inputs['CorrStrength'].default_value = 1
        node_group.outputs.new('NodeSocketShader', 'Shader')
        
        
        # object centric coordinate system
        node = coord = shading_utils.add_node_group(nodes, self.object_coord, pos=self.start)
        loc = node.location

        # statining intensity
        node = staining_noise = nodes.new('ShaderNodeTexMusgrave')
        node.inputs['Scale'].default_value = 5 # TODO
        node.inputs['Detail'].default_value = 15
        node.inputs['Dimension'].default_value = 0
        node.inputs['Lacunarity'].default_value = 1.7
        links.new(inputs.outputs['Size'], staining_noise.inputs['Scale'])
        links.new(inputs.outputs['Dimension'], staining_noise.inputs['Dimension'])
        links.new(inputs.outputs['Lacunarity'], staining_noise.inputs['Lacunarity'])

        loc = node.location = (loc[0]+self.sep, loc[1])
        links.new(coord.outputs[0], staining_noise.inputs['Vector'])
        
        node = threshold = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.operation = 'SUBTRACT'
        links.new(staining_noise.outputs[0], threshold.inputs[0])
        links.new(inputs.outputs['Threshold'], threshold.inputs[1])
        node = softness = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.operation = 'DIVIDE'
        node.inputs[0].default_value = 1
        links.new(threshold.outputs[0], softness.inputs[0])
        links.new(inputs.outputs['Softness'], softness.inputs[1])
        
        node = staining_intensity = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc[0]+self.sep, loc[1])
        links.new(softness.outputs[0], staining_intensity.inputs[0])

        # Strength
        node = strength = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.operation = 'MULTIPLY'
        links.new(inputs.outputs['Strength'], strength.inputs[1])
        links.new(staining_intensity.outputs[0], strength.inputs[0])

        # Mix Shader
        node = mix_shader = nodes.new('ShaderNodeMixShader')
        loc = node.location = (loc[0]+1.5*self.sep, loc[1])
        links.new(strength.outputs[0], mix_shader.inputs['Fac'])

        # # Mix Shader correction
        # node = mix_shader_correction = nodes.new('ShaderNodeMixShader')
        # loc = node.location = (loc[0]+1.5*self.sep, loc[1])
        # links.new(inputs.outputs['Threshold'], mix_shader_correction.inputs['Fac'])
        # links.new(staining_intensity.outputs[0], mix_shader_correction.inputs[1])

        # # Add correction
        # node = add_correction = nodes.new('ShaderNodeAddShader')
        # loc = node.location = (loc[0]+1.5*self.sep, loc[1])
        # links.new(mix_shader.outputs[0], add_correction.inputs[0])
        # links.new(mix_shader_correction.outputs[0], add_correction.inputs[1])
        # node = strength_correction = nodes.new('ShaderNodeMixShader')
        # loc = node.location = (loc[0]+1.5*self.sep, loc[1])
        # links.new(inputs.outputs['CorrStrength'], strength_correction.inputs['Fac'])
        # links.new(add_correction.outputs[0], strength_correction.inputs[1])

        # connect to inputs and outputs
        links.new(inputs.outputs[0], mix_shader.inputs[1])
        outputs.location = (loc[0]+self.sep, loc[1])
        links.new(mix_shader.outputs[0], outputs.inputs[0])

        return node_group


    
    def add_stacked_noise(self, node_name='PrincipleTissue'):
        node_group, inputs, outputs = shading_utils.create_node_group(node_name, start=self.start)
        nodes = node_group.nodes
        links = node_group.links
        node_group.inputs.new('NodeSocketShader', 'Shader')
        node_group.inputs.new('NodeSocketFloat', 'Strength')
        node_group.inputs['Strength'].default_value = 0.6
        node_group.inputs.new('NodeSocketFloat', 'ThresholdRips')
        node_group.inputs['ThresholdRips'].default_value = 0.3
        node_group.inputs.new('NodeSocketFloat', 'ThresholdGrain1')
        node_group.inputs['ThresholdGrain1'].default_value = 0.3
        node_group.inputs.new('NodeSocketFloat', 'ThresholdGrain2')
        node_group.inputs['ThresholdGrain2'].default_value = 0.3
        node_group.inputs.new('NodeSocketFloat', 'StrengthRips')
        node_group.inputs['StrengthRips'].default_value = 1
        node_group.inputs.new('NodeSocketFloat', 'StrengthGrain1')
        node_group.inputs['StrengthGrain1'].default_value = 0.5
        node_group.inputs.new('NodeSocketFloat', 'StrengthGrain2')
        node_group.inputs['StrengthGrain2'].default_value = 0.1
        node_group.outputs.new('NodeSocketShader', 'Shader')
        
        # add rips
        node = rips = shading_utils.add_node_group(nodes, self.principle_noise, pos=self.start)
        loc = node.location
        node.inputs['Size'].default_value = 4
        rips.inputs['Softness'].default_value = 0.2
        rips.inputs['Dimension'].default_value = 1.5
        rips.inputs['Lacunarity'].default_value = 2
        links.new(inputs.outputs['ThresholdRips'], rips.inputs['Threshold'])
        links.new(inputs.outputs['StrengthRips'], rips.inputs['Strength'])
        links.new(inputs.outputs['Shader'], rips.inputs['Shader'])

        # add grain 1
        node = grain1 = shading_utils.add_node_group(nodes, self.principle_noise, pos=(loc[0]+self.sep, loc[1]))
        loc = node.location
        grain1.inputs['Size'].default_value = 40
        grain1.inputs['Softness'].default_value = 0.1
        grain1.inputs['Dimension'].default_value = 0
        grain1.inputs['Lacunarity'].default_value = 1.7
        links.new(inputs.outputs['ThresholdGrain1'], grain1.inputs['Threshold'])
        links.new(inputs.outputs['StrengthGrain1'], grain1.inputs['Strength'])
        links.new(rips.outputs['Shader'], grain1.inputs['Shader'])

        # add grain 2
        node = grain2 = shading_utils.add_node_group(nodes, self.principle_noise, pos=(loc[0]+self.sep, loc[1]))
        loc = node.location
        grain2.inputs['Size'].default_value = 100
        grain2.inputs['Softness'].default_value = 0.1
        grain2.inputs['Dimension'].default_value = 0
        grain2.inputs['Lacunarity'].default_value = 1.7
        links.new(inputs.outputs['ThresholdGrain2'], grain2.inputs['Threshold'])
        links.new(inputs.outputs['StrengthGrain2'], grain2.inputs['Strength'])
        links.new(grain1.outputs['Shader'], grain2.inputs['Shader'])
        links.new(grain2.outputs['Shader'], outputs.inputs['Shader'])

        return node_group



