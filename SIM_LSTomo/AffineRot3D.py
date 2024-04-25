import pyclesperanto_prototype as cle
from skimage.io import imread
from skimage.transform import rotate
import numpy as np
from napari_animated_gif_io import napari_write_image as imsave_animated_gif

import NanoImagingPack as nip
import numpy as np

cle.select_device("RTX")

#%%

%gui qt
import napari

# Create an empty viewer
viewer = napari.Viewer()

#%% 

# CREATE CUBE
Npixel = 256
mVolume = np.abs(nip.xx((Npixel,Npixel,Npixel)))<50
mVolume *= np.abs(nip.yy(mVolume.shape))<50
mVolume *= np.abs(nip.zz(mVolume.shape))<50
mVolume[mVolume==0]=(np.random.rand(mVolume[mVolume==0].shape[0])>.99)
viewer.add_image(mVolume)

original = cle.push(1.*mVolume)
original = original * 254.0 / original.max() + 1
original.min(), original.max()

# rotate volume -axis and store it
# in the target 3D stack
temp = cle.create_like(target)
    
angle=45
transform = cle.AffineTransform3D()
transform.center(original.shape)
transform.rotate_around_y_axis(angle)
transform.center(original.shape, undo=True)
#transform.translate(0,translation,0) # xyz
cle.affine_transform(original, temp, transform=transform)


result = cle.pull(temp)
plt.imshow(np.mean(result,-1)), plt.colorbar()

viewer.add_image()

