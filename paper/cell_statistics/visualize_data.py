import numpy as np
import matplotlib.pyplot as plt

n = 2

### Import images and display an image
imgs = np.load('images.npy')
plt.imshow(imgs[n])
plt.show()

### Import labels and display one
labs = np.load('labels.npy')
plt.imshow(labs[n,:,:,0]) #-------------------------- Why 2 channels? - 1. instance seg; 2. semantic channel 
plt.show()