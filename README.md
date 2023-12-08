# rendered_HE

## Setup Blender and VS Code

1) Download [https://www.blender.org/download/](https://www.blender.org/download/)
2) Follow the steps in: [5 Steps to setup VSCode for Blender Python (on Windows)](https://www.youtube.com/watch?v=YUytEtaVrrc)
3) Open the .py script you want to run, then CTRL + SHIFT + P to open the VS Code options and choose "Blender: Run Script". This should see the script results in your open Blender application.

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

