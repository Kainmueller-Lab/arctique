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
    def __init__(self, start=(0, 0), sep=200, shift=(0, 0, 0)):
        self.shift = shift
        self.start = start
        self.sep = sep
        
        # general
        self.object_coord = self.add_object_coord(shift=self.shift)
        self.volume = self.add_volume()
        self.mixing_red = self.add_mixing_red_points()
        
        # tissue
        self.lamina_propria_tissue_base = self.add_lamina_propria_tissue_base()
    
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
    
    def add_lamina_propria_tissue_base(self, node_name='LaminaPropBase'):
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
        node.color_ramp.elements[0].color = (0.409, 0.215, 0.430, 1)
        node.color_ramp.elements[1].color = (0.456, 0.011, 0.356, 1)
        links.new(musgrave_noise.outputs[0], color_ramp_1.inputs[0])
        links.new(color_ramp_1.outputs[0], color_ramp_2.inputs[0])
        
        # density 
        node = density_ramp = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc_tex[0]+self.sep, loc_tex[1]-2*self.sep)
        node.color_ramp.elements[0].color = (1, 1, 1, 1)
        node.color_ramp.elements[1].color = (0, 0, 0, 1)
        node = density_strength = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+1.5*self.sep, loc[1])
        node.operation = 'MULTIPLY'
        node.inputs[1].default_value = 91.2
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


    def add_principle_noise(self):
        pass
    
    def add_stacked_noise(self):
        pass