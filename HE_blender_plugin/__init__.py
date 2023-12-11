import bpy

bl_info = {
    "name": "runscript",
    "author": "",
    "description": "runs a script",
    "blender": (2, 80, 0),
    "location": "View3D",
    "warning": "",
    "category": "Generic"
}


class RunScript(bpy.types.Operator):
    bl_idname = "object.runscript"
    bl_label = "runscript"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        from pathlib import Path

        # chdir to the parent directory of the current file
        import os
        os.chdir(Path(os.path.dirname(os.path.realpath(__file__))).parent)

        # add the parent directory of the current file to the python path
        import sys
        sys.path.append(str(Path(os.path.dirname(os.path.realpath(__file__))).parent))

        import main  # runs the script
        return {'FINISHED'}


def menu_func(self, context):
    self.layout.operator(RunScript.bl_idname)


def register():
    bpy.utils.register_class(RunScript)
    bpy.types.TOPBAR_MT_file.prepend(menu_func)  # Adds the new operator to an existing menu.


def unregister():
    bpy.utils.unregister_class(RunScript)
    bpy.types.TOPBAR_MT_file.remove(menu_func)


# This allows you to run the script directly from Blender's Text editor
# to test the add-on without having to install it.
if __name__ == "__main__":
    register()
