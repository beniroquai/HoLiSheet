import pyclesperanto_prototype as cle
from skimage.io import imread
from skimage.transform import rotate
import numpy as np
from napari_animated_gif_io import napari_write_image as imsave_animated_gif

import NanoImagingPack as nip
import numpy as np

%gui qt
import napari



from fit3Dsphere import *

# Create an empty viewer
viewer = napari.Viewer()
# select computatonal device
cle.select_device("RTX")



#%% Setup parametres
size3d=(256,256,256)
k0=100; 
kxymax=(70, 70, 70)
img=np.zeros(size3d) #newim(size3d,'scomplex') 
aproj = np.ones((size3d[0],size3d[1]))
NA=.5
lambda0=0.450
pixelsize=lambda0/NA/4
pixelsizes=(pixelsize,pixelsize,pixelsize)
n=1
allowwraparound=True


#  40 subpixel subdivisions
# 2*10+1 kernelsize
# 20 pixel bordersize in all directions
# 1500 iterations
#%% Compute LUT
imatrix=IterateCoefficients(subpix=40,CutoffK=10,TotalSize=size3d[-1],RollOff=20,Iterations=1000)
#%% MAP LUT in 3d space
indexList2D,fullIndex3D,factorList,amask=FillProjSpherePrepare(size3d,lambda0,pixelsizes,NA,imatrix,n=1.0,allowwraparound=1.0)



#%% rotate image in 3D around the rotational axis
Nangles = 10
allAngles = []
for iangle in range(Nangles):
    print(iangle)
    #% Map image to LUT in 3D 
    img = nip.image(size3d)+0j
    img = FillProjSphere(img,aproj,np.int32(indexList2D),fullIndex3D,factorList) # Just a circular wave
    img = np.roll(img, size3d[0]//2, 0)
    #img_Ft = nip.ft(img) 

    
    
    #%    
    original = cle.push(img) # TODO: complex values

    # rotate volume -axis and store it
    # in the target 3D stack
    temp = cle.create_like(target)
        
    angle=iangle*360/Nangles
    transform = cle.AffineTransform3D()
    transform.center(original.shape)
    transform.rotate_around_y_axis(angle)
    transform.center(original.shape, undo=True)
    #transform.translate(0,translation,0) # xyz
    cle.affine_transform(original, temp, transform=transform)
    
    
    result = cle.pull(temp)
    plt.imshow(np.mean(result,-1)), plt.colorbar()
    allAngles.append(result)
    

viewer.add_image(np.mean(np.array(allAngles),0))