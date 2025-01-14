from PIL import Image
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
import numpy as np

from pathlib import Path
import os 


def sanity_check_single_masks(mask_path, target_path, sample_idxs): 

    for IDX in sample_idxs: 
        single_mask_path = mask_path.joinpath(f"single_masks/{IDX}")
        single_masks_idxs = os.listdir(single_mask_path)
        N_masks = len(single_masks_idxs)

        Nx = int(np.floor(np.sqrt(N_masks)))
        Ny = int(np.ceil(N_masks/Nx))

        fig = plt.figure(figsize=(12,12))
        gs = GridSpec(Ny, Nx, figure=fig, hspace=0.05, wspace=0.05)

        for i in range(N_masks): 

            plt.subplot(gs[i//Nx, i%Nx])
            mask = np.array(Image.open(single_mask_path.joinpath(single_masks_idxs[i])))
            plt.imshow(mask, cmap="Greys_r")
            plt.axis("off")

        plt.suptitle(f"sample ID {IDX}", size=16)
        plt.savefig(target_path.joinpath(f"single_masks_sample{IDX}.png"), bbox_inches="tight")



def sanity_check_semantic(mask_path, target_path, sample_idxs): 
    sanity_check_semantic_classes(mask_path, target_path, sample_idxs)
    sanity_check_semantic_noise(mask_path, target_path, sample_idxs)

    data_idx = [int(di.strip(".png")) for di in os.listdir(mask_path.joinpath(f"semantic_noise_0"))]

    with open(target_path.joinpath("semantic_classes.txt"), "w") as file: 
        for idx in data_idx:
            file.write(f"{idx} \t{list(np.unique(np.array(Image.open(mask_path.joinpath(f'semantic_noise_0/{idx}.png')))))} \n")



def sanity_check_semantic_classes(mask_path, target_path, sample_idxs): 
    
    my_cmap = matplotlib.cm.get_cmap('rainbow')
    my_cmap.set_under('black')
    n_samples = len(sample_idxs)

    n_classes = 0
    for IDX in sample_idxs: 
        mask_file = mask_path.joinpath(f"semantic_noise_0/{IDX}.png")
        mask_array = np.array(Image.open(mask_file))
        n_classes = max(n_classes, len(np.unique(mask_array)))

    plt.figure(figsize=(16,10))

    for i, IDX in enumerate(sample_idxs):

        mask_file = mask_path.joinpath(f"semantic_noise_0/{IDX}.png")
        mask_array = np.array(Image.open(mask_file))

        plt.subplot(n_samples,n_classes+1,i*(n_classes+1)+1)
        plt.xticks([])
        plt.yticks([])
        if i==0: plt.title("Full semantic mask")
        plt.imshow(mask_array, interpolation="none", cmap=my_cmap, vmin=1)

        for j, idx in enumerate(np.unique(mask_array)[1:]): 
            plt.subplot(n_samples,n_classes+1,i*(n_classes+1)+j+2)
            plt.imshow(idx*np.ma.masked_array(mask_array, mask_array!=idx), interpolation="none", cmap="Greys_r")
            plt.xticks([])
            plt.yticks([])
            if i==0: plt.title(f"class: {idx}")
        plt.tight_layout()
    plt.savefig(target_path.joinpath("semantic_classes.png"), bbox_inches="tight")


def sanity_check_semantic_noise(mask_path, target_path, sample_idxs): 
    my_cmap = matplotlib.cm.get_cmap('rainbow')
    my_cmap.set_under('black')

    noise_levels = [int(f.split("_")[2]) for f in os.listdir(mask_path) if ((f.startswith("semantic")) and (not f.endswith("indexing"))) ]
    noise_levels.sort()

    for IDX in sample_idxs:

        true_mask = np.array(Image.open(mask_path.joinpath(f"semantic_noise_{0}/{IDX}.png")))
        instance_mask = np.array(Image.open(mask_path.joinpath(f"instance_indexing/{IDX}.png")))

        n_cells = len(np.unique(instance_mask))-1

        plt.figure(figsize=(16, 5))

        for nl_idx, nl in enumerate(noise_levels): 

            mask = np.array(Image.open(mask_path.joinpath(f"semantic_noise_{nl}/{IDX}.png")))

            plt.subplot(2,6,nl_idx+1)
            plt.imshow(mask, cmap=my_cmap, interpolation="none", vmin=1)
            plt.title(f"noise-level {nl}")
            plt.axis("off")
            if nl_idx==0: 
                plt.text(s="noisy label", x=-80, y=500, rotation=90, size=12)

            plt.subplot(2,6,nl_idx+7)
            mislabelled_cells = (mask!=true_mask).astype(int)
            plt.imshow(mislabelled_cells, cmap=my_cmap, interpolation="none", vmin=0.01, vmax=1)
            
            n_misslabelled = len(np.unique(instance_mask * mislabelled_cells))-1
            plt.title(f"frac switched: {np.round(n_misslabelled/n_cells, 2)}")
            plt.axis("off")

            if nl_idx==0: 
                plt.text(s="mislabelled cells ", x=-80, y=500, rotation=90, size=12)

        plt.savefig(target_path.joinpath(f"semantic_noise_{IDX}.png"), bbox_inches="tight")


def sanity_check_FG_BG(mask_path, target_path, sample_idxs):

    noise_levels = [int(f.split("_")[3]) for f in os.listdir(mask_path) if f.startswith("FG") ]
    noise_levels.sort()

    for IDX in sample_idxs:
        plt.figure(figsize=(12,4))

        true_mask_array = np.array(Image.open(mask_path.joinpath(f"FG_BG_noise_0/{IDX}.png")))

        for nl_idx, nl in enumerate(noise_levels): 

            plt.subplot(2, len(noise_levels), nl_idx+1)
            noisy_mask_array = np.array(Image.open(mask_path.joinpath(f"FG_BG_noise_{nl}/{IDX}.png")))
            plt.imshow(noisy_mask_array, cmap="Greys_r", interpolation="none")
            plt.title(f"noise-level {nl}")
            plt.xticks([])
            plt.yticks([])
            if nl_idx==0: plt.ylabel("noisy-mask", size=16)

            plt.subplot(2, len(noise_levels), len(noise_levels)+nl_idx+1)
            plt.imshow(noisy_mask_array!=true_mask_array, cmap="Greys_r", interpolation="none")
            plt.xticks([])
            plt.yticks([])
            if nl_idx==0: plt.ylabel("disagreement", size=16)

        plt.suptitle(f"FG-BG noise sample {IDX}")
        plt.tight_layout()

        plt.savefig(target_path.joinpath(f"FG_BG_noise_{IDX}.png"), bbox_inches="tight")


def sanity_check_instance_3class(mask_path, target_path, sample_idxs):

    noise_levels = [int(f.split("_")[3]) for f in os.listdir(mask_path) if (f.startswith("instance") and (not f.endswith("indexing")))]
    noise_levels.sort()

    for IDX in sample_idxs:

        plt.figure(figsize=(12,4))
        true_mask_array = np.array(Image.open(mask_path.joinpath(f"instance_3class_noise_0/{IDX}.png")))
        noise_levels =  [0, 20, 40, 60, 80, 100]

        for nl_idx, nl in enumerate(noise_levels): 
            plt.subplot(2, len(noise_levels), nl_idx+1)
            noisy_mask_array = np.array(Image.open(mask_path.joinpath(f"instance_3class_noise_{nl}/{IDX}.png")))
            plt.imshow(noisy_mask_array, cmap="Greys_r", interpolation="none")
            plt.xticks([])
            plt.yticks([])
            plt.title(f"noise-level {nl}")
            if nl_idx == 0: plt.ylabel("mask", size=16)

            plt.subplot(2, len(noise_levels), len(noise_levels)+nl_idx+1)
            plt.imshow(true_mask_array!=noisy_mask_array, cmap="Greys_r", interpolation="none")
            plt.xticks([])
            plt.yticks([])
            if nl_idx == 0: plt.ylabel("disagreement", size=16)

        plt.suptitle(f"Border Noise {IDX}")
        plt.tight_layout()

        plt.savefig(target_path.joinpath(f"instance_3class_noise_{IDX}.png"), bbox_inches="tight")