.. _usage:

Usage
=====


Generate an example
-------------------

Let's generate an example image and its corresponding masks as shown in the figure above. To do so, you can simply run the following command in the terminal:

.. code-block:: console

   python render.py --start-idx 42 --n-samples 1


This command will generate a single example with index 42, which serves as the seed for the random generation process.
To generate additional examples, you can increase the ``--n-samples argument``. If you wish to create examples with different seeds, adjust the ``--start-idx`` parameter.

.. note::

   Rendering each image may take a few minutes, so for large-scale dataset generation, we recommend using HPC with GPU resources.


The output image and its corresponding mask will be saved in the rendered folder, following this structure: ::

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

Generate a variation of a rendered example
------------------------------------------

The full power of the ARCTIQUE framework unfolds when varying the scene gradually, which allows to study concepts such as uncertainty in a controlled manner.
To this end, we provide parameter sliders such as:

- example

In order to generate a variation of a rendered example, you can use the following command:

.. code-block:: console

   python render_variation.py