import bpy
import sys
import os


# IMPORT SOURCES
dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir )

import src.objects.macros_structures as macro

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

# test code
macro.build_crypt()