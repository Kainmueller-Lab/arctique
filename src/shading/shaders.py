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
    def __init__(self, start=(0, 0), sep=200, shift=(0, 0, 0), stroma_intensity=1, red_points_strength=0, border_fiber_length=0.7, over_staining=1):
        self.shift = shift
        self.start = start
        self.sep = sep
        
        # general
        self.object_coord = self.add_object_coord(shift=self.shift)
        self.addition = self.add_addition()
        self.volume = self.add_volume()
        self.blood_cells_base = self.add_blood_cells_base()
        self.mixing_red = self.add_mixing_red_points(strength=red_points_strength)
        self.principle_noise = self.add_principle_noise()
        self.stacked_noise = self.add_stacked_noise()
        self.perlin_noise = self.add_perlin_noise()
        
        # tissue
        self.fibers_edge = self.add_edge_detector(node_name='FibersEdge', min=0.331, mid=0.475, max=1, strength=1)
        self.stroma_edge = self.add_edge_detector(node_name='StromaEdge', min=0.422, mid=0.5, max=0.673, strength=0.174)
        self.fibers_base = self.add_fibers_base()
        self.fibers = self.add_fibers()
        self.stacked_fibers = self.add_stacked_fibers(node_name='StackedFibersLarge')
        self.stacked_fibers_tissue = self.add_stacked_fibers(
            node_name='StackedFibersTissue',
            scales=[0.7, 10.8, 12.8, 14.6, 15.3, 27.2, 32.5],
            scale_distortions=[0.5, 35.9, 39, 22.1, 22.1, 22.1, 22.1],
            scale_density_noises=[7.9, 8.7, 5, 5, 7.9, 6, 5])
        self.fiber_network = self.add_fiber_network(fill=0.673*over_staining)
        self.border_fibers = self.add_border_fibers(length=border_fiber_length)
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
    
    def add_edge_detector(
            self, node_name='EdgeDetector', min=0.25, mid=0.5, max=0.75, strength=1):
        
        node_group, inputs, outputs = shading_utils.create_node_group(node_name, start=self.start)
        nodes = node_group.nodes
        links = node_group.links
        node_group.outputs.new('NodeSocketFloat', 'Value')  # between 0 and strength
        node_group.inputs.new('NodeSocketFloat', 'Value')  # between 0 and 1
        node_group.inputs.new('NodeSocketFloat', 'Strength')
        node_group.inputs['Strength'].default_value = strength

        # add coloramp as edge detector
        node = color_ramp = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (self.start[0]+self.sep, self.start[1])
        node.color_ramp.elements[0].position = min
        node.color_ramp.elements[1].position = mid
        # add color ramp element max
        node.color_ramp.elements.new(max)
        node.color_ramp.elements[2].position = max
        node.color_ramp.elements[0].color = (0, 0, 0, 1)
        node.color_ramp.elements[1].color = (1, 1, 1, 1)
        node.color_ramp.elements[2].color = (0, 0, 0, 1)
        links.new(inputs.outputs[0], color_ramp.inputs[0])

        # scale to strength
        node = scale = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        scale.operation = 'MULTIPLY'
        links.new(color_ramp.outputs[0], scale.inputs[0])
        links.new(inputs.outputs[1], scale.inputs[1])

        # connect to output
        links.new(scale.outputs[0], outputs.inputs[0])

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
        node.inputs['ScatterDensity'].default_value = 30
        
        # connect to inputs and outputs
        links.new(vol.outputs[0], outputs.inputs[0])
        
        return node_group
    
    def add_fibers_base(self, node_name='FibersBase'):

        node_group, inputs, outputs = shading_utils.create_node_group(node_name, start=self.start)
        nodes = node_group.nodes
        links = node_group.links
        node_group.outputs.new('NodeSocketFloat', 'Color')
        node_group.inputs.new('NodeSocketFloat', 'Scale')
        node_group.inputs.new('NodeSocketFloat', 'ScaleDistortion')
        node_group.inputs.new('NodeSocketFloat', 'DensityShift')
        node_group.inputs.new('NodeSocketVector', 'CurlVector')
        node_group.inputs['Scale'].default_value = 6
        node_group.inputs['ScaleDistortion'].default_value = 30
        node_group.inputs['DensityShift'].default_value = 0
        node_group.inputs['CurlVector'].default_value = (140, 0, 0)
        loc = self.start
        
        # object centric coordinate system
        node = coord = shading_utils.add_node_group(
            nodes, self.object_coord, pos=(loc[0]+self.sep, loc[1]))
        loc = node.location
        
        # TURBULENCE
        # turbulence basis
        node = noise_1 = nodes.new('ShaderNodeTexNoise')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.noise_dimensions = '4D'
        node.inputs['W'].default_value = 22
        node.inputs['Scale'].default_value = 6
        node.inputs['Detail'].default_value = 0
        node.inputs['Roughness'].default_value = 0.633
        node.inputs['Distortion'].default_value = 0
        node.noise_dimensions = '4D'
        links.new(coord.outputs[0], noise_1.inputs['Vector'])

        # turbulence math
        node = subtract = nodes.new('ShaderNodeVectorMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        subtract.operation = 'SUBTRACT'
        node.inputs[1].default_value = (0.5, 0.5, 0.5)
        links.new(noise_1.outputs['Color'], subtract.inputs[0])
        node = mapping = nodes.new('ShaderNodeMapping')
        loc = node.location = (loc[0]+self.sep, loc[1])        
        node = cross_product = nodes.new('ShaderNodeVectorMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        cross_product.operation = 'CROSS_PRODUCT'
        links.new(inputs.outputs['CurlVector'], mapping.inputs['Vector'])
        links.new(subtract.outputs[0], cross_product.inputs[0])
        links.new(mapping.outputs[0], cross_product.inputs[1])

        # DISTORTION
        # distortion basis
        node = noise_2 = nodes.new('ShaderNodeTexNoise')
        loc = node.location = (coord.location[0]+self.sep, coord.location[1]+self.sep)
        node.inputs['Scale'].default_value = 30
        node.inputs['Detail'].default_value = 0
        node.inputs['Roughness'].default_value = 0.633
        node.inputs['Distortion'].default_value = 0
        links.new(coord.outputs[0], noise_2.inputs['Vector'])

        # distortion math
        node = multiply = nodes.new('ShaderNodeVectorMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        multiply.operation = 'SCALE'
        node.inputs['Scale'].default_value = 7.7
        links.new(noise_2.outputs['Color'], multiply.inputs[0])

        # ADD
        node = add = nodes.new('ShaderNodeVectorMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        add.operation = 'ADD'
        links.new(cross_product.outputs[0], add.inputs[0])
        links.new(multiply.outputs[0], add.inputs[1])

        # APPLY
        # basis
        node = noise_3 = nodes.new('ShaderNodeTexNoise')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.inputs['Scale'].default_value = 0.1
        node.inputs['Detail'].default_value = 11.2
        node.inputs['Roughness'].default_value = 0.5
        node.inputs['Distortion'].default_value = 0
        links.new(add.outputs[0], noise_3.inputs['Vector'])

        # density
        node = density = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        density.operation = 'ADD'
        links.new(noise_3.outputs[0], density.inputs[0])
        links.new(inputs.outputs['DensityShift'], density.inputs[1])

        # sharpening
        node = sharpening = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc[0]+self.sep, loc[1])
        sharpening.color_ramp.elements[0].color = (1, 1, 1, 1)
        sharpening.color_ramp.elements[1].color = (0, 0, 0, 1)
        sharpening.color_ramp.elements[0].position = 0.289
        sharpening.color_ramp.elements[1].position = 0.523
        links.new(density.outputs[0], sharpening.inputs[0])
        
        # connect to inputs and outputs
        links.new(sharpening.outputs[0], outputs.inputs[0])
        links.new(inputs.outputs[0], noise_1.inputs['Scale'])
        links.new(inputs.outputs[1], noise_2.inputs['Scale'])
        links.new(inputs.outputs[2], density.inputs[1])
        
        return node_group
    

    def add_fibers(self, node_name='Fibers'):
        '''
        Takes the base fibers and adds a perlin noise to to control the density/thickness
        '''
        node_group, inputs, outputs = shading_utils.create_node_group(node_name, start=self.start)
        nodes = node_group.nodes
        links = node_group.links
        node_group.outputs.new('NodeSocketFloat', 'Color')
        node_group.inputs.new('NodeSocketFloat', 'Density')
        node_group.inputs.new('NodeSocketFloat', 'Scale')
        node_group.inputs.new('NodeSocketFloat', 'ScaleDistortion')
        node_group.inputs.new('NodeSocketFloat', 'WeightDensityNoise')
        node_group.inputs.new('NodeSocketFloat', 'ScaleDensityNoise')
        node_group.inputs.new('NodeSocketVector', 'CurlVector')
        node_group.inputs['Density'].default_value = 0.5
        node_group.inputs['Scale'].default_value = 6
        node_group.inputs['ScaleDistortion'].default_value = 30
        node_group.inputs['WeightDensityNoise'].default_value = 0.107
        node_group.inputs['ScaleDensityNoise'].default_value = 5
        node_group.inputs['CurlVector'].default_value = (140, 0, 0) 
        loc = self.start

        # object centric coordinate system
        node = coord = shading_utils.add_node_group(
            nodes, self.object_coord, pos=(loc[0]+self.sep, loc[1]))
        loc = node.location

        # density noise
        node = noise_density = nodes.new('ShaderNodeTexNoise')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.inputs['Detail'].default_value = 2
        node.inputs['Roughness'].default_value = 0.5
        node.inputs['Distortion'].default_value = 0
        links.new(coord.outputs[0], noise_density.inputs['Vector'])
        node = weight_density = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        weight_density.operation = 'MULTIPLY'
        links.new(noise_density.outputs[0], weight_density.inputs[0])
        links.new(inputs.outputs['WeightDensityNoise'], weight_density.inputs[1])

        # density offset
        node = density_offset = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.operation = 'MULTIPLY'
        node.inputs[0].default_value = 0.2
        links.new(inputs.outputs['Density'], density_offset.inputs[1])
        node = density_offset_2 = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        density_offset_2.operation = 'ADD'
        node.inputs[1].default_value = -0.1
        links.new(density_offset.outputs[0], density_offset_2.inputs[0])

        # add values
        node = add = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        add.operation = 'ADD'
        links.new(density_offset_2.outputs[0], add.inputs[0])
        links.new(weight_density.outputs[0], add.inputs[1])

        # add fibers base
        node = fibers_base = shading_utils.add_node_group(
            nodes, self.fibers_base, pos=(loc[0]+self.sep, loc[1]))
        loc = node.location
        links.new(add.outputs[0], fibers_base.inputs['DensityShift'])
        links.new(inputs.outputs['Scale'], fibers_base.inputs['Scale'])
        links.new(inputs.outputs['ScaleDistortion'], fibers_base.inputs['ScaleDistortion'])
        links.new(inputs.outputs['CurlVector'], fibers_base.inputs['CurlVector'])

        # connect to inputs and outputs
        links.new(fibers_base.outputs[0], outputs.inputs[0])
        links.new(inputs.outputs['ScaleDensityNoise'], noise_density.inputs['Scale'])

        return node_group

    def add_stacked_fibers(
            self, node_name='StackedFibers',
            scales=[5, 6, 1, 0.7], scale_distortions=[30, 30, 5, 0.5],
            scale_density_noises=[5, 5, 5, 5]):
        '''
        Adds up multiple fiber nodes with different scales
        '''
        node_group, inputs, outputs = shading_utils.create_node_group(node_name, start=self.start)
        nodes = node_group.nodes
        links = node_group.links
        node_group.outputs.new('NodeSocketFloat', 'Color')
        node_group.inputs.new('NodeSocketFloat', 'Density')
        node_group.inputs['Density'].default_value = 0.5
        loc = self.start

        # add fiber nodes
        fibers = []
        for i in range(len(scales)):
            node = fiber = shading_utils.add_node_group(
            nodes, self.fibers, pos=(loc[0]+self.sep, loc[1]))
            loc = node.location
            fiber.inputs['Scale'].default_value = scales[i]
            fiber.inputs['ScaleDistortion'].default_value = scale_distortions[i]
            fiber.inputs['ScaleDensityNoise'].default_value = scale_density_noises[i]
            fibers.append(fiber)

        # add values
        node = add = shading_utils.add_node_group(
            nodes, self.addition, pos=(loc[0]+self.sep, loc[1]))
        loc = node.location
        for i, fiber in enumerate(fibers):
            links.new(fiber.outputs[0], add.inputs[i])
            links.new(inputs.outputs[0], fiber.inputs['Density'])
        
        # connect to inputs and outputs
        links.new(add.outputs[0], outputs.inputs[0])

        return node_group
    
    def add_fiber_network(self, node_name='FiberNetwork', fill=0.673, density=0.7):
        '''
        fiber network with filling tissue in between (softening the effect)
        '''
        node_group, inputs, outputs = shading_utils.create_node_group(node_name, start=self.start)
        nodes = node_group.nodes
        links = node_group.links
        loc = self.start
        node_group.outputs.new('NodeSocketFloat', 'FilledNetwork')
        node_group.outputs.new('NodeSocketFloat', 'Network')
        node_group.inputs.new('NodeSocketFloat', 'VoronoiScale')
        node_group.inputs.new('NodeSocketFloat', 'VoronoiVariation')
        node_group.inputs.new('NodeSocketFloat', 'OverstainingOffset')
        node_group.inputs['VoronoiScale'].default_value = 100
        node_group.inputs['VoronoiVariation'].default_value = 10

        # fiber base
        node = fibers = shading_utils.add_node_group(
            nodes, self.stacked_fibers_tissue, pos=(loc[0]+self.sep, loc[1]))
        loc = node.location
        fibers.inputs['Density'].default_value = density

        # fiber network (adding slightly expanded fibers)
        node = exp_fibers = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc[0]+self.sep, loc[1])
        exp_fibers.color_ramp.elements[1].position = 0.015
        exp_fibers.color_ramp.elements[1].color = (0.709, 0.709, 0.709, 1)
        links.new(fibers.outputs[0], exp_fibers.inputs[0])
        node = add = shading_utils.add_node_group(
            nodes, self.addition, pos=(loc[0]+self.sep, loc[1]))
        links.new(exp_fibers.outputs[0], add.inputs[0])
        links.new(fibers.outputs[0], add.inputs[1])

        # add voronoi fibers
        node = coord = shading_utils.add_node_group(
            nodes, self.object_coord, pos=(loc[0]+self.sep, loc[1]))
        loc = coord.location
        node = scaling_noise = nodes.new('ShaderNodeTexNoise')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.noise_dimensions = '4D'
        node.inputs['W'].default_value = 5.6
        node.inputs['Scale'].default_value = 3.3
        node.inputs['Detail'].default_value = 2
        node.inputs['Roughness'].default_value = 0.5
        node.inputs['Distortion'].default_value = 0
        links.new(coord.outputs[0], scaling_noise.inputs['Vector'])
        node = var = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.operation = 'MULTIPLY'
        links.new(scaling_noise.outputs[0], var.inputs[0])
        links.new(inputs.outputs[1], var.inputs[1])
        node = strength = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        strength.operation = 'ADD'
        links.new(var.outputs[0], strength.inputs[0])
        links.new(inputs.outputs[0], strength.inputs[1])
        node = voronoi = nodes.new('ShaderNodeTexVoronoi')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.feature = 'DISTANCE_TO_EDGE'
        links.new(strength.outputs[0], voronoi.inputs['Scale'])
        links.new(coord.outputs[0], voronoi.inputs['Vector'])
        node = clipping = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.color_ramp.elements[1].position = 0.225
        node.color_ramp.elements[1].color = (0, 0, 0, 1)
        node.color_ramp.elements[0].color = (1, 1, 1, 1)
        links.new(voronoi.outputs[0], clipping.inputs[0])
        links.new(clipping.outputs[0], add.inputs[2])

        # add inverse to fill in between
        node = inverse = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.color_ramp.elements[1].position = 0.382
        node.color_ramp.elements[0].color = (fill, fill, fill, 1)
        node.color_ramp.elements[1].color = (0, 0, 0, 1)
        links.new(add.outputs[0], inverse.inputs[0])
        node = add_inverse = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        add_inverse.operation = 'ADD'
        links.new(inverse.outputs[0], add_inverse.inputs[0])
        links.new(add.outputs[0], add_inverse.inputs[1])

        # connect to outputs
        links.new(add_inverse.outputs[0], outputs.inputs[0])
        links.new(add.outputs[0], outputs.inputs[1])

        return node_group
    
    def add_border_fibers(self, node_name='BorderFibers', length=0.7):
        '''
        adds fibres at border of the tissue
        '''
        node_group, inputs, outputs = shading_utils.create_node_group(node_name, start=self.start)
        nodes = node_group.nodes
        links = node_group.links
        loc = self.start
        node_group.outputs.new('NodeSocketFloat', 'Value')
        node_group.inputs.new('NodeSocketFloat', 'Value')

        # filter region of fibers
        node = filter = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc[0]+self.sep, loc[1])
        base = 0.518
        node.color_ramp.elements[0].position = 0.447
        node.color_ramp.elements[1].position = 0.518
        node.color_ramp.elements.new(base + length*(1-base))
        node.color_ramp.elements.new(1)
        node.color_ramp.elements[0].color = (0, 0, 0, 1)
        node.color_ramp.elements[1].color = (1, 1, 1, 1)
        node.color_ramp.elements[2].color = (1, 1, 1, 1)
        node.color_ramp.elements[3].color = (0, 0, 0, 1)
        links.new(inputs.outputs[0], filter.inputs[0])

        # add fibers
        node = offset = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.operation = 'ADD'
        node.inputs[1].default_value = -1.4
        links.new(filter.outputs[0], offset.inputs[0])
        node = fibers = shading_utils.add_node_group(
            nodes, self.fibers, pos=(loc[0]+self.sep, loc[1]))
        loc = node.location
        links.new(offset.outputs[0], fibers.inputs['Density'])
        node.inputs['Scale'].default_value = 16.6
        node.inputs['ScaleDistortion'].default_value = 34.1
        node.inputs['ScaleDensityNoise'].default_value = 8.2
        node = invert = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.color_ramp.elements[0].color = (1, 1, 1, 1)
        node.color_ramp.elements[1].color = (0, 0, 0, 1)
        links.new(fibers.outputs[0], invert.inputs[0])

        # subtract fibers
        node = subtract = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        subtract.operation = 'SUBTRACT'
        links.new(inputs.outputs[0], subtract.inputs[0])
        links.new(invert.outputs[0], subtract.inputs[1])

        # connect to outputs
        links.new(subtract.outputs[0], outputs.inputs[0])

        return node_group
    
    def add_perlin_noise(self, node_name='PerlinNoise', cutoff_min=0.25, cutoff_max=0.6):
        '''
        adds a perlin noise with a cutoff with object centric coordinate system
        '''
        node_group, inputs, outputs = shading_utils.create_node_group(node_name, start=self.start)
        nodes = node_group.nodes
        links = node_group.links
        node_group.outputs.new('NodeSocketFloat', 'Value')
        node_group.inputs.new('NodeSocketFloat', 'Scale')
        node_group.inputs.new('NodeSocketFloat', 'Detail')
        node_group.inputs.new('NodeSocketFloat', 'Roughness')
        node_group.inputs.new('NodeSocketFloat', 'Distortion')
        node_group.inputs.new('NodeSocketFloat', 'CutoffMin')
        node_group.inputs.new('NodeSocketFloat', 'CutoffMax')
        node_group.inputs['Scale'].default_value = 5.8
        node_group.inputs['Detail'].default_value = 10
        node_group.inputs['Roughness'].default_value = 0.5
        node_group.inputs['Distortion'].default_value = 0
        node_group.inputs['CutoffMin'].default_value = cutoff_min
        node_group.inputs['CutoffMax'].default_value = cutoff_max

        # object centric coordinate system
        node = coord = shading_utils.add_node_group(
            nodes, self.object_coord, pos=self.start)
        loc = node.location
        
        # perlin noise
        node = noise = nodes.new('ShaderNodeTexNoise')
        loc = node.location = (loc[0]+self.sep, loc[1])
        links.new(coord.outputs[0], noise.inputs['Vector'])
        links.new(inputs.outputs['Scale'], noise.inputs['Scale'])
        links.new(inputs.outputs['Detail'], noise.inputs['Detail'])
        links.new(inputs.outputs['Roughness'], noise.inputs['Roughness'])
        links.new(inputs.outputs['Distortion'], noise.inputs['Distortion'])

        # cutoff with new range [0, 1]
        node = subtract = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        subtract.operation = 'SUBTRACT'
        links.new(noise.outputs['Color'], subtract.inputs[0])
        links.new(inputs.outputs['CutoffMin'], subtract.inputs[1])
        node = distance = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        distance.operation = 'SUBTRACT'
        links.new(inputs.outputs['CutoffMax'], distance.inputs[0])
        links.new(inputs.outputs['CutoffMin'], distance.inputs[1])
        node = divide = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        divide.operation = 'DIVIDE'
        links.new(subtract.outputs[0], divide.inputs[0])
        links.new(distance.outputs[0], divide.inputs[1])
        node.use_clamp = True

        # connect to output
        links.new(divide.outputs[0], outputs.inputs[0])

    def add_blood_cells_base(self, node_name='BloodCellsBase', inner_intensity=0.286):
        '''
        creates a bw base for blood cells
        '''
        node_group, inputs, outputs = shading_utils.create_node_group(node_name, start=self.start)
        nodes = node_group.nodes
        links = node_group.links
        node_group.outputs.new('NodeSocketFloat', 'Value')
        node_group.inputs.new('NodeSocketFloat', 'Density')  # uniform voronoi packed 
        node_group.inputs.new('NodeSocketFloat', 'Scale')
        node_group.inputs['Density'].default_value = 28    
        loc = self.start

        # object centric coordinate system
        node = coord = shading_utils.add_node_group(
            nodes, self.object_coord, pos=(loc[0]+self.sep, loc[1]))
        loc = node.location

        # voronoi compartements
        node = noise = nodes.new('ShaderNodeTexVoronoi')
        loc = node.location = (loc[0]+self.sep, loc[1])
        links.new(coord.outputs[0], noise.inputs['Vector'])
        links.new(inputs.outputs['Density'], noise.inputs['Scale'])
        node = noise_subtract = nodes.new('ShaderNodeTexVoronoi')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.feature = 'SMOOTH_F1'
        links.new(coord.outputs[0], noise_subtract.inputs['Vector'])
        links.new(inputs.outputs['Density'], noise_subtract.inputs['Scale'])
        node.inputs['Smoothness'].default_value = 1
        node = cell = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        cell.operation = 'SUBTRACT'
        links.new(noise.outputs['Distance'], cell.inputs[0])
        links.new(noise_subtract.outputs['Distance'], cell.inputs[1])

        # cutoff offset/scale
        node = cutoff = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.color_ramp.elements[0].position = 0.276
        node.color_ramp.elements[1].position = 0.629
        links.new(inputs.outputs['Scale'], cutoff.inputs[0])
        node = scale = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        scale.operation = 'ADD'
        links.new(cutoff.outputs[0], scale.inputs[0])
        links.new(cell.outputs[0], scale.inputs[1])

        # make visibility/intensity of cells
        node = shading = nodes.new('ShaderNodeValToRGB')
        loc = node.location = (loc[0]+self.sep, loc[1])
        node.color_ramp.elements[0].position = 0
        node.color_ramp.elements[1].position = 0.027
        node.color_ramp.elements.new(0.098)
        node.color_ramp.elements.new(0.156)
        node.color_ramp.elements[0].color = (inner_intensity, inner_intensity, inner_intensity, 1)
        node.color_ramp.elements[1].color = (1, 1, 1, 1)
        node.color_ramp.elements[2].color = (1, 1, 1, 1)
        node.color_ramp.elements[3].color = (0, 0, 0, 1)
        links.new(scale.outputs[0], shading.inputs[0])

        # connect to output
        links.new(shading.outputs[0], outputs.inputs[0])

        return node_group
    
    def add_mixing_red_points(self, node_name='MixingRedPoints', strength=0.5):
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
        node.inputs['AbsorptionDensity'].default_value = 125 * (1+0.3*strength)
        node.inputs['ScatterDensity'].default_value = 0.6
        
        # object centric coordinate system
        node = coord = shading_utils.add_node_group(
            nodes, self.object_coord,
            pos=(self.start[0]+self.sep, self.start[1]))
        loc = node.location
        
        # abundance/strength off blood cells
        node = staining_noise = nodes.new('ShaderNodeTexNoise')
        node.inputs['Scale'].default_value = 5
        node.inputs['Detail'].default_value = 0
        node.inputs['Roughness'].default_value = 0
        node.inputs['Distortion'].default_value = 0
        loc = node.location = (loc[0]+self.sep, loc[1])
        node = scale = nodes.new('ShaderNodeMath')
        loc = node.location = (loc[0]+self.sep, loc[1])
        scale.operation = 'MULTIPLY'
        scale.inputs[1].default_value = 1.5-strength
        links.new(coord.outputs[0], staining_noise.inputs['Vector'])
        links.new(staining_noise.outputs[0], scale.inputs[0])

        # add blood cells
        node = blood_cells = shading_utils.add_node_group(
            nodes, self.blood_cells_base,
            pos=(loc[0]+self.sep, loc[1]))
        loc = node.location
        links.new(scale.outputs[0], blood_cells.inputs['Scale'])
        
        # Mix Shader
        node = mix_shader = nodes.new('ShaderNodeMixShader')
        loc = node.location = (loc[0]+1.5*self.sep, loc[1])
        links.new(vol.outputs[0], mix_shader.inputs[2])
        links.new(blood_cells.outputs[0], mix_shader.inputs['Fac'])
        
        # connect to inputs and outputs
        links.new(inputs.outputs[0], mix_shader.inputs[1])
        outputs.location = (loc[0]+self.sep, loc[1])
        links.new(mix_shader.outputs[0], outputs.inputs[0])
        
        return node_group
        
    def add_mixing_red_points_legacy(self, node_name='MixingRedPoints', strength=0):
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
        node.inputs['AbsorptionDensity'].default_value = 125 * (1+0.3*strength)
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
        node.color_ramp.elements[0].position = 0.559 - strength
        node.color_ramp.elements[1].position = 0.75 - strength
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

    def add_addition(self, node_name='Addition', n=7):
        '''
        stacks multiple math 'Add' nodes to add multiple scalar values
        '''
        node_group, inputs, outputs = shading_utils.create_node_group(node_name, start=self.start)
        nodes = node_group.nodes
        links = node_group.links
        loc = self.start
        node_group.outputs.new('NodeSocketFloat', 'Value')
        last_node = node_group
        for i in range(n):
            node = nodes.new('ShaderNodeMath')
            node.operation = 'ADD'
            node.location = (loc[0]+self.sep, loc[1])
            node.inputs[0].default_value = 0
            node.inputs[1].default_value = 0
            loc = node.location
            node_group.inputs.new('NodeSocketFloat', 'Value'+str(i))
            links.new(inputs.outputs[i], node.inputs[0])
            if i == 0:
                links.new(node.outputs[0], outputs.inputs[0])
            else:
                links.new(node.outputs[0], last_node.inputs[1])
            last_node = node
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



