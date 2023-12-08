# rendered_HE

## Setup Blender and VS Code

1) Download [https://www.blender.org/download/](https://www.blender.org/download/)
2) Follow the steps in: [5 Steps to setup VSCode for Blender Python (on Windows)](https://www.youtube.com/watch?v=YUytEtaVrrc).
In particular best create a virtual environment where you need to install
```
pip install fake-bpy-module-latest
```
4) Open the .py script you want to run, then CTRL + SHIFT + P to open the VS Code options and choose "Blender: Run Script". This should see the script results in your open Blender application.


## Setup Blender and Pycharm
1) Download [https://www.blender.org/download/](https://www.blender.org/download/)
2) Install blend-charm plugin via: [https://github.com/BlackStartx/PyCharm-Blender-Plugin](https://github.com/BlackStartx/PyCharm-Blender-Plugin) (see releases, make sure you have the correct pycharm version for the plugin!)
3) Follow isntructions in the repository.

## ToDos
- Create a simple end-to-end pipeline script (Generate blob, render it, export png) and check if this works without using Blender extension. This is important in case other users unfamiliar with Blender want to generate images.

## Pipeline
1) Generation of cells
- different classes
- slider between classes (similarity in shape and size)
- deletion of cells if the overlap in 3d-space
2) Add realism with 3D ray-tracing
- structured tissue
- variable staining
- cut-off cells
3) Lightning
- generate image modality (e.g. lightsheed, fluorescence)
4) Rendering: render same image with
- full shading
- classical instance masks
- depth masks (or 3D labels), with more information useful to uncertainty quantification
5) Post processing / add. rendering

