# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 16:29:21 2022

@author: bene
"""

import numpy as np 
import NanoImagingPack as nip
import matplotlib.pyplot as plt



#%gui qt
import napari
# Create an empty viewer
viewer = napari.Viewer()

#%% data saving

def savemat(name, array):
    import numpy, scipy.io
    scipy.io.savemat(name+".mat", mdict={name+"_py": array})

def loadmat(name):
    import numpy, scipy.io
    scipy.io.loadmat(name+".mat")


#%%
from fit3Dsphere import *
size3d=(256,256,256)

#  40 subpixel subdivisions
# 2*10+1 kernelsize
# 20 pixel bordersize in all directions
# 1500 iterations
imatrix=IterateCoefficients(subpix=40,CutoffK=10,TotalSize=size3d[-1],RollOff=20,Iterations=1000)
savemat("imatrix", imatrix)
#%%
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

import numpy, scipy.io
aa=scipy.io.loadmat('/Users/bene/Dropbox/Dokumente/Promotion/MATLAB/Toolboxes/ElectricFields/indexList2D.mat')
indexList2D=aa['indexList2D']
fullIndex3D=aa['fullIndex3D']
factorList=aa['factorList']
 
 
indexList2D,fullIndex3D,factorList,amask=FillProjSpherePrepare(size3d,lambda0,pixelsizes,NA,imatrix,n=1.0,allowwraparound=1.0)


#%%
img=nip.image(size3d)+0j
img = FillProjSphere(img,aproj,np.int32(indexList2D),fullIndex3D,factorList) # Just a circular wave
img_Ft = nip.ft(img)


savemat('indexList2D',np.int32(indexList2D))
savemat('fullIndex3D',np.int32(fullIndex3D))
savemat('factorList',np.int32(factorList))
                    
#plt.imshow(np.sum(np.abs(img),0)),plt.colorbar(),plt.show()
plt.imshow(np.log(1+np.abs(img[:,size3d[1]//2,:])))



#%%


import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.ndimage import map_coordinates

# Rotates 3D image around image center
# INPUTS
#   array: 3D numpy array
#   orient: list of Euler angles (phi,psi,the)
# OUTPUT
#   arrayR: rotated 3D numpy array
# by E. Moebel, 2020
def rotate_array(array, orient):
    phi = orient[0]
    psi = orient[1]
    the = orient[2]

    # create meshgrid
    dim = array.shape
    ax = np.arange(dim[0])
    ay = np.arange(dim[1])
    az = np.arange(dim[2])
    coords = np.meshgrid(ax, ay, az)

    # stack the meshgrid to position vectors, center them around 0 by substracting dim/2
    xyz = np.vstack([coords[0].reshape(-1) - float(dim[0]) / 2,  # x coordinate, centered
                     coords[1].reshape(-1) - float(dim[1]) / 2,  # y coordinate, centered
                     coords[2].reshape(-1) - float(dim[2]) / 2])  # z coordinate, centered

    # create transformation matrix
    r = R.from_euler('zxz', [phi, psi, the], degrees=True)
    mat = r.as_matrix()

    # apply transformation
    transformed_xyz = np.dot(mat, xyz)

    # extract coordinates
    x = transformed_xyz[0, :] + float(dim[0]) / 2
    y = transformed_xyz[1, :] + float(dim[1]) / 2
    z = transformed_xyz[2, :] + float(dim[2]) / 2

    x = x.reshape((dim[1],dim[0],dim[2]))
    y = y.reshape((dim[1],dim[0],dim[2]))
    z = z.reshape((dim[1],dim[0],dim[2])) # reason for strange ordering: see next line

    # the coordinate system seems to be strange, it has to be ordered like this
    new_xyz = [y, x, z]

    # sample
    arrayR = map_coordinates(array, new_xyz, order=1)
    return arrayR


img_rot = rotate_array(np.log(1+np.abs(img)), (1,np.pi/4,np.pi/4)) 

from scipy.ndimage import rotate

#%%
all_axes = [(1, 0), (1, 2), (0, 2)]
degree=45
image_rotate=rotate(np.array(np.abs(img.copy())), angle=34, reshape=False, axes=all_axes[0])

plt.imshow(np.sum(np.log(1+image_rotate),2))