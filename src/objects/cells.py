import bpy 
import bmesh
import random
import math

from typing import List, Optional

from mathutils import Vector

from src.shading.shading import Material
from src.utils.helper_methods import set_orientation

class CellAttribute():
    def __init__(self):
        self.cell_type = None
        self.size = None
        self.scale = None
        self.deformation_strength = None
        self.attribute_name = None
        self.max_bending_strength = None

class CellAttributeA(CellAttribute):
    def __init__(self, cell_type = "A", size = 0.04, scale = (1,1,1), deformation_strength = 0.02, attribute_name = "Cell Type A", max_bending_strength = 0.2):
        self.cell_type = cell_type
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class CellAttributeB(CellAttribute):
    def __init__(self, cell_type = "B", size = 0.07, scale = (2.3,1.3,1), deformation_strength = 0.05, attribute_name = "Cell Type B", max_bending_strength = 0.6):
        self.cell_type = cell_type
        self.size = size
        self.scale = scale
        self.deformation_strength = deformation_strength
        self.attribute_name = attribute_name
        self.max_bending_strength = max_bending_strength

class Cell:
    cell_count = 0

    def __init__(self, location: Vector, arrangement_id: int, arrangement_type: str, cell_attributes: CellAttribute, orientation: Optional[Vector] = None):
        self.cell_attributes = cell_attributes
        self.cell_id = Cell.cell_count
        Cell.cell_count += 1
        self.location = location
        self.orientation = orientation
        self.cell_name = f"Cell_{self.cell_id}_Arr_{arrangement_type}_{arrangement_id}_Type_{self.cell_attributes.cell_type}"
        self.pass_index = self.cell_id
        self.arrangement_type = None
        self.semantic_id = None
        self.material = None
        self.cell_object = None

    def set_material(self, material: Material):
        self.material = material

    def set_material(self, semantic_id: int):
        self.semantic_id = semantic_id

    def add(self):
        # Create a cube
        bpy.ops.mesh.primitive_cube_add(size=self.cell_attributes.size, location=self.location, scale=self.cell_attributes.scale)
        self.cell_object = bpy.context.active_object
        self.cell_object.name = self.cell_name
        
        pass_index = bpy.props.IntProperty(name="Pass Index", subtype='UNSIGNED')
        bpy.context.object.pass_index = self.pass_index
        
        # Orient
        if self.orientation:
            set_orientation(self.cell_object, self.orientation)

        # Apply two levels of subdivision surface
        modifier = self.cell_object.modifiers.new("Subsurface Modifier", "SUBSURF")
        modifier.levels = 2    

        # Deform the vertices proportionally in a random fashion
        self.deform_mesh_old() # NOTE: This is just the random deformation method, the other needs moe care. - ck
        # Bend mesh along the z axis
        self.bend_mesh()


    '''
    translate random mesh vertices using proportional edit
    i.e.: neighboring vertices are also translated proportionally
    '''
    def deform_mesh(self):
        TRANSLATION_RANGE = 0.4
        TRANSFORM_COUNT = 3
        PROPORTIONAL_SIZE = 1       

        # apply subsurface modifier
        bpy.ops.object.modifier_apply({"object": bpy.context.object}, modifier="Subsurface Modifier")
        #extract mesh
        bpy.ops.object.mode_set(mode='EDIT')

        context = bpy.context
        bm = bmesh.from_edit_mesh(context.edit_object.data)
        #deselect all vertices
        for v in bm.verts:
            v.select = False
        bm.verts.ensure_lookup_table() # NOTE: Necessary for accessing vertex list
        for i in range(TRANSFORM_COUNT):
            transform = Vector([random.uniform(-1, 1),
                        random.uniform(-1, 1),
                        random.uniform(-1, 1)])*TRANSLATION_RANGE
            v = bm.verts[random.randint(0, len(bm.verts) - 1)]
            v.select = True
            bpy.ops.transform.translate(value=transform, 
                                    constraint_axis=(False, False, False),
                                    orient_type='GLOBAL',
                                    mirror=False, 
                                    use_proportional_edit = True,
                                    use_proportional_connected =True,
                                    proportional_edit_falloff='SMOOTH',
                                    proportional_size=PROPORTIONAL_SIZE)
            v.select = False

        # # get a reference to the active object
        # mesh_obj = bpy.context.active_object

        # bpy.context.view_layer.objects.active = mesh_obj

        # ###### INITIALIZE BMESH ######
        # # create a new bmesh
        # bm = bmesh.new()
        # # initialize the bmesh data using the mesh data
        # bm.from_mesh(mesh_obj.data)

        # ###### EDIT BMESH ######
        # bm.verts.ensure_lookup_table() # NOTE: Necessary for accessing vertex list
        # for _ in range(TRANSFORM_COUNT):
        #     transform = Vector([random.uniform(-1, 1),
        #                 random.uniform(-1, 1),
        #                 random.uniform(-1, 1)])*TRANSLATION_RANGE
        #     v = bm.verts[random.randint(0, len(bm.verts) - 1)]
        #     v.select_set(True)
        #     print(bpy.context.selected_objects)
        #     v.co += transform
        #     # bpy.ops.transform.translate(value=transform, 
        #     #                         constraint_axis=(False, False, False),
        #     #                         orient_type='GLOBAL',
        #     #                         mirror=False, 
        #     #                         use_proportional_edit = True,
        #     #                         use_proportional_connected =True,
        #     #                         proportional_edit_falloff='SMOOTH',
        #     #                         proportional_size=PROPORTIONAL_SIZE)
        #     # print("Vertex translated")
        #     #v.select_set(True)

        # # bmesh.ops.bevel(
        # #     bm,
        # #     geom=bm.edges,
        # #     offset=0.2,
        # #     segments=4,
        # #     affect="EDGES",
        # #     profile=0.5,
        # # )

        # ###### UPDATE BMESH ######
        # bm.normal_update()
        # # writes the bmesh data into the mesh data
        # bm.to_mesh(mesh_obj.data)
        # # [Optional] update the mesh data (helps with redrawing the mesh in the viewport)
        # mesh_obj.data.update()
        # # clean up/free memory that was allocated for the bmesh
        # bm.free()

        # vertices = self.cell_object.data.vertices
        # # Deselect all vertices
        # for v in vertices:
        #     pass
        # #v.select = False
        # # Transform random vertices
        # for _ in range(TRANSFORM_COUNT):
        #     for obj in bpy.context.selected_objects:
        #         obj.select_set(False)   
        #     transform = Vector([random.uniform(-1, 1),
        #                 random.uniform(-1, 1),
        #                 random.uniform(-1, 1)])*TRANSLATION_RANGE
        #     random_index = random.randint(0, len(vertices) - 1)
        #     v = vertices[random_index]
        #     print_context_details()
            #v.select = True
            # bpy.ops.transform.translate(value=transform, 
            #                         constraint_axis=(False, False, False),
            #                         orient_type='LOCAL',
            #                         use_proportional_edit = True,
            #                         use_proportional_connected =True,
            #                         proportional_edit_falloff='SMOOTH',
            #                         proportional_size=PROPORTIONAL_SIZE)
            #v.select = False
        #bpy.ops.object.mode_set(mode='OBJECT')

    def bend_mesh(self):
        # Bend mesh along Z axis
        modifier = self.cell_object.modifiers.new("Simple Deform Modifier", "SIMPLE_DEFORM")
        modifier_index = len(self.cell_object.modifiers) - 1  # Index of the last added modifier
        modifier = self.cell_object.modifiers[modifier_index]
        modifier.deform_method = 'BEND'

        # Create an empty
        bpy.ops.object.empty_add(type='ARROWS', align='WORLD', location=self.cell_object.location, scale=(1, 1, 1))
        empty = bpy.context.active_object

        # Set the origin and deform axis
        modifier.origin = empty
        modifier.deform_axis = 'Z'
        bending_strength = random.uniform(-1,1)*self.cell_attributes.max_bending_strength
        modifier.angle = 2*math.pi*bending_strength

        # Remove empty object
        bpy.data.objects.remove(empty, do_unlink=True)

    def deform_mesh_old(self):
        # Iterate through each vertex and deform its position
        for vertex in self.cell_object.data.vertices:
            original_position = vertex.co.copy()
            deformation_vector = Vector([
                random.uniform(-1, 1),
                random.uniform(-1, 1),
                random.uniform(-1, 1)
            ])*Vector(self.cell_attributes.scale)*self.cell_attributes.deformation_strength
            vertex.co = original_position + deformation_vector