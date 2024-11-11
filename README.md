# ARCTIQUE

ARCTIQUE: ARtificial Colon Tissue Images for Quantitativ Uncertainty Evaluation

With this tool you can generate synthetic images of microscopic colon tissue images, which A) exhibit a complexity comparable to real images, B) can be manipulated and re-rendered in a controlled fashion, and C) are rendered together with the exact pixelwise nuclei masks; instance and semantic. This makes Arctique a powerful tool for the evaluation of concepts such as uncertainty quantification.

As shown below, to mimic the complex structure of real colon tissue images we a) generate a 3D model of the tissue, b) slice through the tissue and populate it with nuclei of different cell types, c) render the 2D image and its instance and semantic masks (2D and 3D).

![Graphical Abstract](examples/generation_overview.png#gh-light-mode-only)
![Graphical Abstract](examples/generation_overview_dark.png#gh-dark-mode-only)




## Installation

1) Clone the repository
2) Make a new conda/micromamba environment with python 3.10 and activate it
```bash
micromamba create -n arctique python=3.10.0 -c conda-forge
```
```bash	
micromamba activate arctique
```

2) Install the requirements via pip
```bash
pip install -r requirements.txt
```

done :sparkles: Now we can use blender without the need to install the full suite.

We tested ARCTIQUE on Linux and Windows.


## Quickstart

### Generate an example
 
Let's generate an example image and its corresponding masks such as:

![Examples](examples/fig_examples.png)

To do so, you can simply run the following command in the terminal:

```bash
python render.py --start-idx 42 --n-samples 1
```

This will render exactly one example with index 42, which implies the seed in the reandom generation process. To render more subsequent examples, you can increase the `--n-samples` argument. If you want to render examples with different seeds, you can change the `--start-idx` argument. NOTE: The rendering process can take a couple minutes per image. Hence, we recommend to use HPC with GPU resources for generating large-scale datasets.

The rendered image and its corresponding mask will be saved in the `rendered` folder with the following structure:


```bash
render
├── images
│   ├── image_0.png
│   ├── ...
└── masks
    ├── instance
    │   ├── 0.tif
    │   ├── ...
    ├── semantic
    │   ├── 0.tif
    │   ├── ...
    ├── cytoplasm
    │   ├── 0.tif
    │   ├── ...
    ├── instance_3d
    │   ├── 0
    │   │   ├── 0_0.tif
    │   │   ├── 0_1.tif
    │   │   ├── ...
    │   │   ├── 0_stack.npy
```

### Generate a variation of a rendered example

The full power of the ARCTIQUE framework unfolds when varying the scene gradually, which allows to study concepts such as uncertainty in a controlled manner. To this end, we provide parameter sliders to control e.g. the nuclei intensity:

![Graphical Abstract](examples/fig_nuclei_intensity.png)


In order to generate a variation of a rendered examples, you can place the original parameter files in the `rendered/parameters_og` folder and the .json file with the altered parameters (e.g. `changed_params.json`) in the `rendered/parameters_changed` folder. Then, you can run the following command in the terminal

```bash
python render_variation.py --file-changed-parameters changed_params.json
```

and it will re-render all the files in the 'parameters_og' folder with the altered parameters. The rendered images and masks will be saved in a folder with the same name as the changed parameter file, e.g. `rendered/changed_params`.
