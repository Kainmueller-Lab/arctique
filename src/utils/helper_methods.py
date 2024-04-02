import bpy

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

def get_info_from_cell_name(name):
    id = name.split('_')[1]
    type = name.split('_')[3]
    return id, type