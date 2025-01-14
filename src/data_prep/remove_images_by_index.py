
import shutil
import os 
from pathlib import Path

bad_images = [4, 28, 54, 55, 73, 75, 121, 130, 170]

FOLDER = Path("/home/claudia/synthetic_uncertainty/data/render_v1_4/train")

print(os.listdir(FOLDER.joinpath("masks")))



for img_idx in bad_images: 
    img_name = f"img_{img_idx}.png"
    file = FOLDER.joinpath("images").joinpath(img_name)
    os.remove(file)
    print("removed:", file)

    for subfolder in ['instance', 'semantic_indexing', 'semantic', 'instance_indexing']: 
        mask_name = f"{img_idx}.png"
        file = FOLDER.joinpath("masks").joinpath(subfolder).joinpath(mask_name)
        os.remove(file)
        print("removed:", file)

    for subfolder in ['instance_3d', 'instance_3d_indexing']: 
        folder = FOLDER.joinpath("masks").joinpath(subfolder).joinpath(f"{img_idx}/")
        shutil.rmtree(folder)
        print("removed:", folder)
    

