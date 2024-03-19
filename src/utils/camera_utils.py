import bpy
import random
import mathutils

def get_cam(cam=None):
    """
    cam: Camera or str, if is str bpy.data.get(str, camera_add())
    """
    if cam is None:
        cam = bpy.context.scene.camera
    if isinstance(cam, str):
        name = cam
        if name not in bpy.data.objects:
            current_obj = bpy.context.view_layer.objects.active
            bpy.ops.object.camera_add()
            assert len(bpy.context.selected_objects) >= 1, f"Camera_add() failed!"
            cam = bpy.context.selected_objects[-1]
            cam.name = name
            bpy.context.view_layer.objects.active = current_obj
        cam = bpy.data.objects[name]
    return cam

def get_cam_intrinsic(cam=None):
    """
    Refrence:
        https://blender.stackexchange.com/a/120063/86396
    """
    # BKE_camera_sensor_size
    def get_sensor_size(sensor_fit, sensor_x, sensor_y):
        if sensor_fit == "VERTICAL":
            return sensor_y
        return sensor_x

    # BKE_camera_sensor_fit
    def get_sensor_fit(sensor_fit, size_x, size_y):
        if sensor_fit == "AUTO":
            if size_x >= size_y:
                return "HORIZONTAL"
            else:
                return "VERTICAL"
        return sensor_fit

    cam = get_cam(cam)
    camd = cam.data
    if camd.type != "PERSP":
        raise ValueError("Non-perspective cameras not supported")
    scene = bpy.context.scene
    f_in_mm = camd.lens
    scale = scene.render.resolution_percentage / 100
    resolution_x_in_px = scale * scene.render.resolution_x
    resolution_y_in_px = scale * scene.render.resolution_y
    sensor_size_in_mm = get_sensor_size(
        camd.sensor_fit, camd.sensor_width, camd.sensor_height
    )
    sensor_fit = get_sensor_fit(
        camd.sensor_fit,
        scene.render.pixel_aspect_x * resolution_x_in_px,
        scene.render.pixel_aspect_y * resolution_y_in_px,
    )
    pixel_aspect_ratio = scene.render.pixel_aspect_y / scene.render.pixel_aspect_x
    if sensor_fit == "HORIZONTAL":
        view_fac_in_px = resolution_x_in_px
    else:
        view_fac_in_px = pixel_aspect_ratio * resolution_y_in_px
    pixel_size_mm_per_px = sensor_size_in_mm / f_in_mm / view_fac_in_px
    s_u = 1 / pixel_size_mm_per_px
    s_v = 1 / pixel_size_mm_per_px / pixel_aspect_ratio

    # Parameters of intrinsic calibration matrix K
    u_0 = resolution_x_in_px / 2 - camd.shift_x * view_fac_in_px
    v_0 = resolution_y_in_px / 2 + camd.shift_y * view_fac_in_px / pixel_aspect_ratio
    skew = 0  # only use rectangular pixels

    K = mathutils.Matrix(((s_u, skew, u_0), (0, s_v, v_0), (0, 0, 1)))
    return K