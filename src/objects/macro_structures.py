import bpy
from src.shading.materials import Material



class build_crypt():
    def __init__(self, name='crypt_structure'):
        self.name = name
        self.crypt = self._add_geometry()
        self.node_tree, self.input, self.output = self._add_geometry_node()
        self.hex_input, self.hex_output = self._add_hexagon_structure(self.input.outputs['Geometry'])
        self.crypt_input, self.crypt_otput = self._add_crypts(self.hex_output, self.output.inputs['Geometry'])
        bpy.ops.object.modifier_apply(modifier=self.name)

        # add crypt volumes
        self.crypt_vol_in = self._make_crypt_vol(thickness=0.01, name='crypt_volume_inner')
        self.crypt_vol_out = self._make_crypt_vol(thickness=0.02, name='crypt')

        objects = [self.crypt, self.crypt_vol_in, self.crypt_vol_out]
        for obj in objects:
            #self._cut_geometry(obj)
            self._scale(obj, (10, 10, 10))
            # bpy.context.view_layer.objects.active = obj
            # obj.select_set(True)
            # bpy.ops.object.modifier_apply(modifier="Solidify")

    def _scale(self, obj, scale):
        obj.scale.x = scale[0]
        obj.scale.y = scale[1]
        obj.scale.z = scale[2]
        bpy.context.view_layer.objects.active = obj
        bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
    
    def _make_crypt_vol(self, thickness=0.1, name='crypt_vol'):
        # copy the crypt
        crypt_vol = self.crypt.copy()
        crypt_vol.data = self.crypt.data.copy()
        bpy.context.collection.objects.link(crypt_vol)
        crypt_vol.name = name

        # add solidification modifier
        bpy.context.view_layer.objects.active = crypt_vol        
        bpy.ops.object.modifier_add(type='SOLIDIFY')
        bpy.context.object.modifiers['Solidify'].thickness = thickness
        bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
        bpy.ops.object.modifier_apply(modifier="Solidify")

        return crypt_vol

    def _add_geometry(self):
        bpy.ops.mesh.primitive_plane_add(size=2, location=(0, 0, 0))
        geometry = bpy.context.active_object
        geometry.name = self.name
        return geometry
        
    def _add_geometry_node(self):
        geo_nodes_modifier = self.crypt.modifiers.new(name=self.name, type='NODES')
        node_tree = bpy.data.node_groups.new(name=self.name, type='GeometryNodeTree')
        geo_nodes_modifier.node_group = node_tree
        nodes = node_tree.nodes
        group_input = nodes.new(type='NodeGroupInput')
        group_input.location = (-400, 0)
        group_output = nodes.new(type='NodeGroupOutput')
        group_output.location = (400, 0)
        node_tree.outputs.new('NodeSocketGeometry', 'Geometry')
        node_tree.inputs.new('NodeSocketGeometry', 'Geometry')
        node_tree.links.new(group_output.inputs['Geometry'], group_input.outputs['Geometry'])
        return node_tree, group_input, group_output
    
    def _add_hexagon_structure(self, in_link=None, out_link=None, start_pos=(0,-400), sep=200):
        nodes = self.node_tree.nodes
        links = self.node_tree.links
        
        # Subdivide mesh for many crypts
        subdivide = nodes.new(type='GeometryNodeSubdivideMesh')
        subdivide.location = start_pos
        subdivide.inputs['Level'].default_value = 4 
        if in_link is not None:
            links.new(in_link, subdivide.inputs['Mesh'])
        
        # scale along diagonal for equilateral triangle after triangulation
        hex_scaling = nodes.new(type='GeometryNodeScaleElements')
        hex_scaling.location = (subdivide.location[0]+sep, subdivide.location[1])
        hex_scaling.scale_mode = 'SINGLE_AXIS'
        hex_scaling.inputs['Scale'].default_value = 1.73205
        hex_scaling.inputs['Axis'].default_value = (1, 1, 0) 
        links.new(subdivide.outputs['Mesh'], hex_scaling.inputs['Geometry'])
        
        # triangulate quads in mesh
        tri = nodes.new(type='GeometryNodeTriangulate')
        tri.location = (hex_scaling.location[0]+sep, hex_scaling.location[1])
        links.new(hex_scaling.outputs['Geometry'], tri.inputs['Mesh'])
        
        # apply dual mesh for hexagons as mesh basis
        dual_mesh = nodes.new(type='GeometryNodeDualMesh')
        dual_mesh.location = (tri.location[0]+sep, tri.location[1])
        links.new(tri.outputs['Mesh'], dual_mesh.inputs['Mesh'])
        
        # make mesh quadratic again
        cube_bool = nodes.new(type='GeometryNodeMeshCube')
        cube_bool.location = (dual_mesh.location[0], dual_mesh.location[1]-sep)
        mesh_bool = nodes.new(type='GeometryNodeMeshBoolean')
        mesh_bool.location = (dual_mesh.location[0]+sep, dual_mesh.location[1])
        mesh_bool.operation = 'INTERSECT'
        links.new(dual_mesh.outputs['Dual Mesh'], mesh_bool.inputs['Mesh 2'])
        
        if out_link is not None:
            links.new(mesh_bool.outputs['Mesh'], out_link)
        
        return subdivide.inputs['Mesh'], mesh_bool.outputs['Mesh']
    
    def _add_crypts(self, in_link=None, out_link=None, start_pos=(1500, -400), sep=200, sep_y=400):
        nodes = self.node_tree.nodes
        links = self.node_tree.links
        
        ### A) Rough Geometry
        # first extrude to get a new face
        extrude = nodes.new(type='GeometryNodeExtrudeMesh')
        extrude.inputs['Offset Scale'].default_value = 0
        extrude.location = start_pos
        if in_link is not None:
            links.new(in_link, extrude.inputs['Mesh'])
        
        # scale extrusion for smaller base
        scale = nodes.new(type='GeometryNodeScaleElements')
        scale.location = (extrude.location[0]+sep, extrude.location[1])
        scale.inputs['Scale'].default_value = 0.3
        links.new(extrude.outputs['Mesh'], scale.inputs['Geometry'])
        links.new(extrude.outputs['Top'], scale.inputs['Selection'])
        
        # face area factor
        face_area = nodes.new(type='GeometryNodeInputMeshFaceArea')
        face_area.location = (start_pos[0], start_pos[1]-sep_y)
        factor = nodes.new(type='ShaderNodeMath')
        factor.operation = 'MULTIPLY'
        factor.location = (face_area.location[0]+sep, face_area.location[1])
        factor.inputs[1].default_value = 150
        links.new(face_area.outputs['Area'], factor.inputs['Value'])
        
        # extrude crypts
        extrude_2 = nodes.new(type='GeometryNodeExtrudeMesh')
        extrude_2.location = (scale.location[0]+sep, scale.location[1])
        links.new(factor.outputs['Value'], extrude_2.inputs['Offset Scale'])
        links.new(scale.outputs['Geometry'], extrude_2.inputs['Mesh'])
        links.new(extrude.outputs['Top'], extrude_2.inputs['Selection'])
        
        # scale extrusion 
        scale_2 = nodes.new(type='GeometryNodeScaleElements')
        scale_2.location = (extrude_2.location[0]+sep, extrude_2.location[1])
        scale_2.inputs['Scale'].default_value = 3.5
        links.new(extrude_2.outputs['Mesh'], scale_2.inputs['Geometry'])
        links.new(extrude_2.outputs['Top'], scale_2.inputs['Selection'])
        
        ### B) Noisy distortion
        # noise texture
        noise = nodes.new(type='ShaderNodeTexNoise')
        noise.location = (extrude_2.location[0], extrude_2.location[1]-sep_y)
        noise.inputs['Scale'].default_value = 1.8
        noise.inputs['Roughness'].default_value = 0.0
        
        # scale (zero out z-axis)
        factor_2 = nodes.new(type='ShaderNodeVectorMath')
        factor_2.location = (noise.location[0]+sep, noise.location[1])
        factor_2.operation = 'MULTIPLY'
        factor_2.inputs[1].default_value = [0.3, 0.3, 0]
        links.new(noise.outputs['Color'], factor_2.inputs[0])
        
        # apply
        pos = nodes.new(type='GeometryNodeSetPosition')
        pos.location = (scale_2.location[0]+sep, scale_2.location[1])
        links.new(scale_2.outputs['Geometry'], pos.inputs['Geometry'])
        links.new(factor_2.outputs['Vector'], pos.inputs['Offset'])
        
        ### C) Increase details
        subdivide = nodes.new(type='GeometryNodeSubdivisionSurface')
        subdivide.location = (pos.location[0]+sep, pos.location[1])
        subdivide.inputs['Level'].default_value = 3
        links.new(pos.outputs['Geometry'], subdivide.inputs['Mesh'])
        
        if out_link is not None:
            links.new(subdivide.outputs['Mesh'], out_link)
            
        return extrude.inputs['Mesh'], subdivide.outputs['Mesh']
    
    def _cut_geometry(self, mesh_object, size=1):
        # Ensure the context is correct
        bpy.ops.object.mode_set(mode='OBJECT')

        # 1. Create a helper cube
        bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0))
        helper_cube = bpy.context.object

        # 2. Apply a Boolean Modifier to 'mesh'
        boolean_modifier = mesh_object.modifiers.new(name='BooleanOp', type='BOOLEAN')
        boolean_modifier.solver = 'FAST'
        boolean_modifier.object = helper_cube
        boolean_modifier.operation = 'INTERSECT' 

        # 3. Apply the modifier
        bpy.context.view_layer.objects.active = mesh_object
        bpy.ops.object.modifier_apply(modifier=boolean_modifier.name)

        # 4. Delete the helper cube
        bpy.data.objects.remove(helper_cube, do_unlink=True)
    



class build_muscosa():
    def __init__(self, crypt, name='muscosa', buffer=(0.2, 0.2, 0.1)):
        self.name = name
        self.buffer = buffer
        self.size = crypt.dimensions
        self.crypt = crypt
        self.muscosa = self._add_geometry()
        self._remove_crypts()

    def _add_geometry(self):
        bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, self.size[2]/2+0.0001))
        geometry = bpy.context.active_object
        geometry.name = self.name
        geometry.scale = (self.size[0]*(1-self.buffer[0]), self.size[1]*(1-self.buffer[1]), self.size[2]*(1+self.buffer[2]))
        bpy.context.scene.cursor.location = (0, 0, 0)
        bpy.ops.object.origin_set(type='ORIGIN_CURSOR', center='MEDIAN')
        bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
        return geometry
    
    def _remove_crypts(self):
        # Apply a Boolean Modifier to 'mesh'
        boolean_modifier = self.muscosa.modifiers.new(name='BooleanOp', type='BOOLEAN')
        boolean_modifier.object = self.crypt
        boolean_modifier.operation = 'DIFFERENCE' 

        # Apply the modifier
        bpy.context.view_layer.objects.active = self.muscosa
        bpy.ops.object.modifier_apply(modifier=boolean_modifier.name)