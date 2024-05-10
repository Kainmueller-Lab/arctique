import bpy


def clear_scene():
    # Deselect all objects
    bpy.ops.object.select_all(action='DESELECT')

    # Delete all objects in the scene
    for obj in bpy.context.scene.objects:
        bpy.data.objects.remove(obj, do_unlink=True)

    # Clear all meshes (geometry data)
    for mesh in bpy.data.meshes:
        bpy.data.meshes.remove(mesh, do_unlink=True)

    # Clear all materials
    for material in bpy.data.materials:
        bpy.data.materials.remove(material, do_unlink=True)

    # Clear all textures
    for texture in bpy.data.textures:
        bpy.data.textures.remove(texture, do_unlink=True)

    # Clear all node groups (includes Geometry Node trees)
    for node_group in bpy.data.node_groups:
        bpy.data.node_groups.remove(node_group, do_unlink=True)


def delete_objects():
    # Select all objects 
    for obj in bpy.context.scene.objects:
        obj.hide_viewport = False
    bpy.ops.object.select_all(action='SELECT')
    # Delete the selected objects
    bpy.ops.object.delete()
    # Remove orphaned meshes
    for mesh in bpy.data.meshes:
        bpy.data.meshes.remove(mesh, do_unlink=True)

def remove_objects(object_list):
    for obj in object_list:
        bpy.data.objects.remove(obj, do_unlink=True)

def get_objects_with(string):
    # Create a list to store matching objects
    object_list = []
    # Iterate through all objects
    for obj in bpy.data.objects:
        # Check if the object name starts with the defined pattern
        if obj.name.startswith(string):
            object_list.append(obj)
    return object_list

def subdivide_object(obj, level, type='SIMPLE'):
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.modifier_add(type='SUBSURF')
    bpy.context.object.modifiers["Subdivision"].subdivision_type = type
    bpy.context.object.modifiers["Subdivision"].levels = level
    bpy.ops.object.modifier_set_active(modifier="Subdivision")
    bpy.ops.object.modifier_apply(modifier="Subdivision")

def get_info_from_cell_name(name):
    id = name.split('_')[1]
    type = name.split('_')[2]
    return id, type