# Arctique

Arctique: ARtificial Colon Tissue for Quantitativ Uncertainty Evaluation

## Installation via pip and requirements

1) Clone the repository
2) Make a new conda/micromamba environment with python 3.10
```bash
micromamba create -n arctique python=3.10.0
```
```bash
conda create -n arctique python=3.10.0
```

2) Install the requirements via pip
```bash
pip install -r requirements.txt
```

done :) Now we can use blender without the need to install the full suite.


## Quickstart

### Generate an image with the corresponding mask

In order to generate an image with the corresponding mask, you can use the following command:

```bash
python render.py 
```

This will output a rendered image and its corresponding mask in the `rendered` folder with the following structure:

```bash
render
├── image
│   ├── image_0.png
│   ├── image_1.png
│   ├── ...
└── mask
    ├── mask_0.png
    ├── mask_1.png
    ├── ...
```

### Generate a variation of a rendered example

In order to generate a variation of a rendered example, you can use the following command:

```bash
python render_variation.py 
```