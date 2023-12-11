class Camera:
    pass


class LightSource:
    pass


class BioMedicalScene:

    def __init__(self, light_source: LightSource, camera:Camera):
        self.light_source = light_source
        self.camera = camera

    def render(self):
        pass
