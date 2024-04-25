import math
import napari
from skimage.registration import phase_cross_correlation
from skimage import filters
import NanoImagingPack as nip
import matplotlib.pyplot as plt
import numpy as np
import cv2
import sys
import tifffile as tif
#!pip install opencv-contrib-python
#!pip install opencv-python

import pyclesperanto_prototype as cle
from skimage.io import imread
from skimage.transform import rotate
import numpy as np
from napari_animated_gif_io import napari_write_image as imsave_animated_gif

import NanoImagingPack as nip


%gui qt
import napari



from fit3Dsphere import *

# Create an empty viewer
viewer = napari.Viewer()
# select computatonal device
cle.select_device("RTX")


# Read
filename = 'SingleParticleInFocus4_rec_Holography-3.tif'
filename = '10h02m44s_rec_Holography-1.tif'
frames = tif.imread(filename)

# %% Aligning the tube to the center -> anti wiggeling
allShifts = []




for iframe in range(len(frames)-1):
    src = frames[iframe][0:300, :]
    dst = cv2.Canny(src, 0, 50)
    lines = cv2.HoughLines(dst, 1, np.pi / 180, 150, None, 0, 0)

    if lines is not None:
        ypos = np.squeeze(lines[0])[0]
    else:
        ypos = 0
    allShifts.append(ypos)

    frame = cv2.resize(frames[iframe], dsize=None, dst=None, fx=.25, fy=.25)

# fit curve to found shifts using Hough transform 
timetrace = np.array(allShifts)-np.min(np.array(allShifts))
timefct = np.arange(0, timetrace.shape[0])
model = np.poly1d(np.polyfit(timefct, timetrace, 2))

plt.plot(model(timefct))
plt.plot(timetrace)
plt.show()

allBackshifted = []
resizeFactor = 0.25
for iframe in range(len(frames)):
    frame = cv2.resize(frames[iframe],None, None, resizeFactor,resizeFactor)
    allBackshifted.append(np.roll(frame, -int(model(iframe)*resizeFactor), 0))


allBackshifted=np.array(allBackshifted)

viewer.add_image(allBackshifted)
viewer.add_image(frames)

tif.imwrite("allBackshifted.tif", allBackshifted)

#%% Map rotation to 3D fourier space Setup parametres
allBackshiftedMid = nip.extract(allBackshifted, (allBackshifted.shape[0],np.min(allBackshifted.shape[1:]),np.min(allBackshifted.shape[1:])), (154,allBackshifted.shape[2]//2))
size3d=(allBackshiftedMid.shape[1],allBackshiftedMid.shape[1],allBackshiftedMid.shape[2])
lambda0=.405
k0=2*np.pi/lambda0
ftVolume=np.zeros(size3d) #newim(size3d,'scomplex') 
aproj = np.ones((size3d[0],size3d[1]))
NA=0.4
pixelsize=0.5#2.4#lambda0/NA/2
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
Nangles = allBackshifted.shape[0]
allAngles = []
for iangle in range(Nangles):
    print(iangle)
    #% Map image to LUT in 3D 
    ftVolume = nip.image(size3d)+0j
    aproj = nip.ft(allBackshiftedMid[iangle])
    ftVolume = FillProjSphere(ftVolume,aproj,np.int32(indexList2D),fullIndex3D,factorList) # Just a circular wave
    ftVolume = np.roll(ftVolume, size3d[0]//2, 0)
    
    ftVolumeR = np.real(ftVolume)
    ftVolumeI = np.imag(ftVolume)
    
    # convert to CL type
    targetR = cle.create(size3d)
    targetI = cle.create(size3d)
    cle.set(targetR, 0)
    cle.set(targetI, 0)    
    originalR = cle.push(ftVolumeR) 
    originalI = cle.push(ftVolumeI) 

    # rotate volume -axis and store it
    # in the target 3D stack      
    angle=iangle*180/Nangles
    transform = cle.AffineTransform3D()
    transform.center(originalR.shape)
    transform.rotate_around_x_axis(angle)
    transform.center(originalR.shape, undo=True)
    #transform.translate(0,translation,0) # xyz
    cle.affine_transform(originalR, targetR, transform=transform)
    cle.affine_transform(originalI, targetI, transform=transform)
    
    
    resultFt = cle.pull(targetR)+1j*cle.pull(targetI)
    result = nip.ift(resultFt)
    
    if iangle ==0:
        allAngles = np.abs(result)
    else:
        allAngles += np.abs(result)
    
    #viewer.add_image(np.abs(allAngles))
    viewer.add_image(np.abs(cle.pull(originalR)))   
    #viewer.add_image(np.abs(result))   
    

#viewer.add_image(np.mean(np.array(allAngles),0))

viewer.add_image(allAngles)
