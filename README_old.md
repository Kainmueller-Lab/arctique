# Arctique

Arctique: ARtificial Colon Tissue for Quantitativ Uncertainty Evaluation

rHEnder is a package that lets you create synthetic biomedical data, together with their corresponding masks.
The package is specifically designed for H&E stained biomedical images, but could be easily 
used outside its original scope. Thereby, rHEnder allows you to maximally control this creation process and 
offers an easy interface for you to explore dataset creation.

This flexibility allows you to render various biomedical realistically looking images to study various research question
in a downstream analysis. Here, we demonstrate how to use the dataset for two major challenges from the machine learning
world: Domain Adaption and Uncertainty Prediction.

For now, let's dive into the installation procedure.

## Setup Blender and VS Code

1) Download [https://www.blender.org/download/](https://www.blender.org/download/)
2) Follow the steps in: [5 Steps to setup VSCode for Blender Python (on Windows)](https://www.youtube.com/watch?v=YUytEtaVrrc).
In particular best create a virtual environment where you need to install
```
pip install fake-bpy-module-latest
```
4) Open the .py script you want to run, then CTRL + SHIFT + P to open the VS Code options and choose "Blender: Run Script". This should see the script results in your open Blender application.
5) Some versions of the script might use the package scipy. In case Blender trows an unkwown module error about scipy although you have scipy install in your local python the problem might be that scipy is not installed in Blender's python. To solve that issue go to the bin folder of Blender's python, e.g. "C:\Program Files\Blender Foundation\Blender 3.6\3.6\python\bin" and then run
```
python.exe -m pip install scipy
```

NOTE: Order to install python libraries inside Blender follow the description here: [pip install in Blender](https://blender.stackexchange.com/questions/56011/how-to-install-pip-for-blenders-bundled-python)


## Setup Blender and Pycharm
1) Download [https://www.blender.org/download/](https://www.blender.org/download/)
2) Install blend-charm plugin via: [https://github.com/BlackStartx/PyCharm-Blender-Plugin](https://github.com/BlackStartx/PyCharm-Blender-Plugin) (see releases, make sure you have the correct pycharm version for the plugin!)
3) Follow instructions in the repository.
4) Mark the folder named "HE_blender_plugin" as plugin directory via rightclick -> new -> Blender-Charm -> Mark as Addon Folder
5) Now every time you start blender via Blend-Charm it will automatically install the plugin in the "HE_blender_plugin" folder and you can ecute the plugin.

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

## Dataset creation from the command line - In testing phase
1) Open the local Git repo in an IDE (commonly VS code)  
2) Via terminal install any missing packages in the Phyton bin folder of the Blender App. Example for Macos:
  ```    
  $ /Applications/Blender.app/Contents/Resources/3.6/python/bin/python3.10 -m pip install <package name>
  ```
4) For a proper functioning, Blender 3.6 and bpy 3.4.0 are required
5) Execute the Python script ```main.py``` with the argument -h to display the options to insert. Recommended line for testing purposes: 
  ```    
  $ /Applications/Blender.app/Contents/Resources/3.6/python/bin/python3.10 main.py --output_dir <user directory> --n_sample 10
  ```

