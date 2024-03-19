import bpy
import src.objects.macro_structures as macro


class TissueArch():
    def __init__(self):
        # 1) build epithelial layer
        # 1a) build crypt
        self.crypt = macro.build_crypt()
        # 1b) build surronding muscosa layer
        self.tissue = macro.build_muscosa(self.crypt.crypt)

    def get_architecture(self):
        return self.crypt.crypt, self.tissue.muscosa
