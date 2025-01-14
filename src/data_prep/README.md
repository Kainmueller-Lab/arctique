# Preparation of labels for training

Arctique provides exact semantic and instance labels in `tif` format. For using our training pipline these labels need to be converted into png format. For foreground-background segementation binary masks must be created from semantic or instance masks. 

The training pipeline and the uncertainty pipeline both expect the masks to be in `.png` format and located in a folder named `{task}_noise_{noise-percent}` e.g. `semantic_noise_0`, `FG_BG_noise_20` or `instance_noise_100`. 

## Generation of noisy labels: 

While Arctique comes with exact labels for semantic and instance segmentation, we can retrospectively perturb the labels to explicitly study the effect of noisy labels on segmentation performance and uncertainty. 

For generating noisy labels, run `python preprocess_data_and_labels.py`. More detailed parameters for the generation of noisy labels can be specified under `noisy_data_config.yaml`. The task specific noise types and the corresponding aparemeters are described below. The most important task-independent parameters are: 

- `data_source` : location of the Arctique dataset
- `target_location`: where to save the new noisy labels 
- `make_single_masks`: whether to extract masks of individual cells from the instance masks. This is quite time consuming and only needed for the generation of noisy FG-BG masks. Should be `false` by default. 
- `create_new_train_val_set`: whether the train set should be split into a train and a validations set. Further details about the split can be specified in the secondary parameters. 

## Noisy labels by segmentataion task: 
Different segementation tasks can be affected by different types of label noise. 

### Noisy Semantic Masks
The most relevant type of noise in semnatic masks is the missclassification of cell types. In order to introduce errors in the semnatic masks you can set the parameter `global_flip_prob` which corresponds to the probablity that two randomly chosen cells switch their class label (`global_flip_prob=1` means that cells have comeplety random labels). 
If you want to confuse only certain cell types you can specify the dictionary `class_flip_dict`. Here the keys are cell types and the values are tuples containing the probability of that cell type being affected by a label-flip and a list of cell types it could be confused with. eg: the key-value pair `3: (0.3, [2,4])` means with probability 0.3 a cell of type 3 is assigned the label 2 or 4. 

### Noisy Instance Masks
In instance segmentation the most common source of noise is the wrongful delineation of cell-instances. This can happen either by false merges where multiple separate instances are incorrectly characterized as one; or by false splits where a single instance is incrorrectly characterized as multiple instances. 

Currenty, the generation of noisy labels only supports false merges. To generate the noisy mask, we first identify all pairs/clusters of cells that are directy next to each other. The parameter `noise_level` then determines what fraction of all possible merges will in fact be merged. 

**Note:** if there are no directly adjacent cells, even maximum `noise_level` will not generate false merges. 

## Noisy Foreground Masks
For generating noisy masks for foreground-background segmentation we mainly aim to manipulate the shape of the single cell masks. We use `scipy.ndimage` to rescale or shift the single-cell masks or perform elastic transforms. We can also completely remove individual masks to mimic an incomplete segmentation mask. For all transformations we specify the maximum possible value of transformation (e.g the maximum number of pixels to shift a mask) and the base-probability to apply this transformation to any single cell mask. The probabilities for different transformations are mutually independent i.e a cell can experience a mask shift and a elastic transform. The `noise_level` parameters rescales the probabilities of applying any transformation but leaves the maimum values unchanged. 

**Note:** In order to create noisy FG-BG masks we first need to extract masks of individual cells and then iterate through  the one by one applying the transformations. Depending on the number of cells in the original mask, this can take a substantial amount of time. 