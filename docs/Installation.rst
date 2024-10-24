Installation
============

Make sure you have `micromamba <https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html>`_ installed.

Then, clone the repository `arctique <https://github.com/Kainmueller-Lab/arctique.git>`_ via git or download via browser.

Create a new conda/micromamba environment with python 3.10 and activate it.

.. code-block:: console

    micromamba create -y -n arctique python=3.10 pip -c conda-forge
    micromamba activate arctique

Install the requirements via pip

.. code-block:: console

    pip install -r requirements.txt

done :sparkles: Now we can use blender without the need to install the full suite.
