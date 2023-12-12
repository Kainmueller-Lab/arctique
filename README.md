# rHEnder

rHEnder is a package that lets you create synthetic biomedical data, together with their corresponding masks.
The package is specifically designed for H&E stained biomedical images, but could be easily 
used outside its original scope. Thereby, rHEnder allows you to maximally control this creation process and 
offers an easy interface for you to explore dataset creation.

This flexibility allows you to render various biomedical realistically looking images to study various research question
in a downstream analysis. Here, we demonstrate how to use the dataset for two major challenges from the machine learning
world: Domain Adaption and Uncertainty Prediction.

For now, let's dive into the installation procedure.

## Automatic installation with album solution

Install album [https://album.solutions/](https://album.solutions/)
Install the solution form this repository:

`album install rendered_HE/album_solution`

Use with:

`album run rHEnder`

## Manual Installation
1) Download [https://www.blender.org/download/](https://www.blender.org/download/)
2) install rHEnder in the blender python via commandline in the following way:
    Linux:
       `blender_install_path/4.0/python/bin/python3.10 rendered_HE/album_solution/install_rhender.py`
    Windows:
        `blender_install_path\\4.0\\python\\bin\\python.exe rendered_HE/album_solution/install_rhender.py`
    Mac:
        `blender_install_path/Contents/Resources/4.0/python/bin/python3.10 rendered_HE/album_solution/install_rhender.py`
3) Open blender and run the main.py file in blenders scripting block.

## Setup Blender and VS Code for development

1) Download [https://www.blender.org/download/](https://www.blender.org/download/)
2) Follow the steps in: [5 Steps to setup VSCode for Blender Python (on Windows)](https://www.youtube.com/watch?v=YUytEtaVrrc).
3) Install the package with pip in the development option (e.g. run `pip install -e .` in the repository base environment)
In particular best create a virtual environment where you need to install
```
pip install fake-bpy-module-latest
```
4) Open the .py script you want to run, then CTRL + SHIFT + P to open the VS Code options and choose "Blender: Run Script". This should see the script results in your open Blender application.


## Setup Blender and Pycharm
1) Download [https://www.blender.org/download/](https://www.blender.org/download/)
2) Install blend-charm plugin via: [https://github.com/BlackStartx/PyCharm-Blender-Plugin](https://github.com/BlackStartx/PyCharm-Blender-Plugin) (see releases, make sure you have the correct pycharm version for the plugin!)
3) Follow instructions in the repository.
4) Install rHEnder in the python you use to develop (probably a conda environment) with the `-e` option.
5) Mark the folder named "HE_blender_plugin" as plugin directory via rightclick -> new -> Blender-Charm -> Mark as Addon Folder
6) now every time you start blender via Blend-Charm it will automatically install the plugin in the "HE_blender_plugin" folder and you can ecute the plugin.

## Testing:
We use pytest for our testing: [https://docs.pytest.org/en/7.4.x/](https://docs.pytest.org/en/7.4.x/).
To test this suite install pytest in the environment where you additionally installed this repository.

Then run on the cmdline:
```    
pytest
```

## ToDos

### General
- Create a simple end-to-end pipeline script (Generate blob, render it, export png) and check if this works without using Blender extension. This is important in case other users unfamiliar with Blender want to generate images.
- Create config file that contains all input parameters, so that parameters and run script are decoupled

### Cell geometry and distribution
- Generate custom cell field shapes (e.g. rectangle, disc field) in which cells are randomly distributed
- Define a field density parameter (for equidistribution on a cell field)
- Implement different distributions (so far only: equidistribution on field shape; add: Gaussian distribution)
- Pass a list of cell positions in case one needs non-random hand-picked cell positions.
- Create more cell type classes if necessary (use specific cell names?)
- Implement interpolation between cellTypes (to get more and less similar cellTypes)
- Toggle cell occlusion (Yes/No)
- Implement minimum (or average) distance between cells

### Mask creation
- Add label members to cell type class (which will be needed to create (in)correctly labeled GT masks)
- Add parameter between 0 and 1 that steers the ratio of mislabeled or unlabeled cell masks

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

