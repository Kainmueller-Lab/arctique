# ARCTIQUE Documentation
Short description of how to locally look at the documentation after you have done some changes. 
Especially those that refer to the notebooks.


## Build and check the documentation
1. install conda
2. install an environment:
```
micromamba create -n sphinx-arctique python=3.10 pip
micromamba activate sphinx-arctique
pip install -r docs/requirements.txt
```
3. make sure "pandoc" is installed (on ubuntu: `sudo apt install pandoc`)
4. go to docs folder: `cd docs/`
5. build the documentation locally with `make html`. A folder "_build" should appear (e.g. docs/_build). 
This is your local version of the documentation!
6. check for errors in the make process - resolve eventual errors
7. open your browser and drag and drop docs/_build/html/index.html into your window.
This will be the view users have on readthedocs.org once you push changes to the docs branch
8. check documentation locally for errors and style. Especially the notebook files! Repeat step 5-8 until satisfied.
9. DO NOT ADD AND THEN PUSH THE "_build" DIRECTORY! IT SHOULD NOT BE IN THE GITHUB REPOSITORY OF ARCTIQUE!




