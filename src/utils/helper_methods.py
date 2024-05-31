import bpy
import numpy as np


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
    cell_part = name.split('_')[0]
    type_idx = name.split('_')[-1]
    cell_type = name.split('_')[2]
    unique_string = f"{cell_type}_{type_idx}"
    return cell_part, type_idx, cell_type, unique_string

def get_universal_cell_ids(cell_part_names):
    '''
    each cell part has a name like 'Nucleus_Type_EPI_0' 
    where the last number is the index for that cell type
    this function returns look-up table for unique indices for all cells
    Args:
        cell_part_names: list of strings
    Returns:
        look_up_indices: dict, like {'EPI_0': 1, 'LYM_0': 2, ...}
    '''
    look_up_indices = {}
    idx = 1
    for name in cell_part_names:
        _, _, _, idx_name = get_info_from_cell_name(name)
        if idx_name not in look_up_indices.keys():
            look_up_indices[idx_name] = idx
            idx += 1
    return look_up_indices

def map_16bit_to_index(values, sep=55, tol=10):
    '''
    maps back the values blender assigns integers when saving a 16 bit image
    back to the original indices
    Args:
        values: iteretable of integers
    Returns:
        indices as list of integers
    '''
    values = np.array(values)
    indices = []
    v, i = 0, 0
    while v < np.max(values)+sep:
        filtered = values[values >= v - tol]
        filtered = filtered[filtered <= v + tol]
        if len(filtered) == 0:
            v = v + sep
        else:
            indices.append(i)
            v = filtered[0] + sep
        i += 1
    return indices


def get_cell_location(cell, camera, resolution):
    '''
    returns the location of a cell in pixel coordinates
    Args:
        camera_obj: bpy object, camera object -> quadratic fov
        cell_obj: bpy object, cell object
        resolution: int -> quadratic fov
    Returns:
        loc_object_space: tuple, (x, y)
        loc_pixel_space: tuple, (x, y)
    '''
    scale = camera.data.ortho_scale
    camera_pos = camera.location
    loc_object_space = cell.location
    loc_camera_space = (np.array(loc_object_space) - np.array(camera_pos)) / scale
    loc_camera_space[1] = -loc_camera_space[1]
    loc_pixel_space = (loc_camera_space + 0.5) * resolution
    return list(loc_object_space)[:2], list(loc_pixel_space)[:2]


def copy_object(obj, name=None):
    '''
    copies an object and links it to the scene
    '''
    obj_copy = obj.copy()
    obj_copy.data = obj.data.copy()
    bpy.context.collection.objects.link(obj_copy)
    if name:
        obj_copy.name = name
    return obj_copy


def convert2mesh(obj):
    '''
    converts visual object to mesh object
    '''
    # unselect all objects
    bpy.ops.object.select_all(action='DESELECT')
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.convert(target='MESH')


def add_boolean_modifier(obj, target, name='Boolean Modifier', operation='DIFFERENCE', apply=True):
    '''
    adds a boolean modifier to an object
    '''
    boolean = obj.modifiers.new(name=name, type='BOOLEAN')
    boolean.operation = operation
    boolean.object = target
    boolean.use_self = True # NOTE: This line seems to fix artifacts when using 'DIFFERENCE' operator to create goblet volume. - ck
    bpy.context.view_layer.objects.active = obj
    if apply:
        bpy.ops.object.modifier_apply(modifier=boolean.name)


def check_obj_empty(obj, threshold_polygons=0):
    '''
    checks if the object is empty
    '''
    return len(obj.data.polygons) <= threshold_polygons


def object_in_hull(obj, hull):
    '''
    checks if the object is fully enclosed by hull
    Args:
        obj: bpy object
        hull: bpy object
    Returns:
        bool
    '''
    obj_copy = copy_object(obj)
    convert2mesh(obj_copy)
    add_boolean_modifier(obj_copy, hull, name='check_inside')
    convert2mesh(obj_copy)
    me = obj_copy.data
    print("Polygons: %d" % len(me.polygons))
    # for poly in me.polygons:
    #     print("Polygon index: %d, length: %d" % (poly.index, poly.loop_total))

    enclosed = check_obj_empty(obj_copy)
    bpy.data.objects.remove(obj_copy)
    return enclosed


def delete_cells_outside_tissue(cells, tissue):
    '''
    deletes all cells that are outside the tissue
    Args:
        cells: list of bpy objects
        tissue: bpy object
    '''
    remaining_cells = []
    for cell in cells:
        if not object_in_hull(cell, tissue):
            print(f"Removed {cell.name}")
            bpy.data.objects.remove(cell)
        else:
            remaining_cells.append(cell)
            print(f"Kept {cell.name}")
    return remaining_cells
    

def shade_switch(obj, flat=True):
    '''
    changes the shading of the object
    Args:
        flat: bool, if True flat shading is used, otherwise smooth shading
        material: bpy material, material to assign
    '''
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.shade_flat() if flat else bpy.ops.object.shade_smooth()
    