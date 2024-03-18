import bpy
import src.objects.macro_structures as macro


class TissueArch():
    def __init__(self):

        # TODO inccoperate scene parameters for cutting

        # 1) build epithelial layer
        # 1a) build crypt
        self.crypt = macro.build_crypt()
        # 1b) build surronding muscosa layer
        self.tissue = macro.build_tissue(self.crypt.crypt)
