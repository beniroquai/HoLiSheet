#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 11:29:21 2022

@author: bene
"""
import NanoImagingPack as nip
import numpy as np
import matplotlib.pyplot as plt
#import h5py
import scipy
import tifffile as tif
nip.setViewer("NIP_VIEW")


#%%

%gui qt
import napari

# Create an empty viewer
viewer = napari.Viewer()
#%%

from stl import mesh
import stl


def propagate_gaussian(size, sampling, lambda0, n0=1., w0=1., zoffset=0):
    """
    propagate_gaussian propagates a gaussian peak through vacuum
    :param size: Input size (Nz, Nx)
    :param sampling: Pixelsizes (dz, dx)
    :param lambda0: wavelength of the gaussian peak
    :param w0: waist of the gaussian peak
    :param: optional: zrange (zmin,zmax) for z-position
    :return: field (ZX)
    """
    Nz,Nx=size[0],size[1]
    dz,dx=sampling[0],sampling[1]
    Dx = Nx*dx                  # Dimesion in X-direction (in mum)
    Dz = Nz*dz                  # Dimesion in Z-direction (in mum)

    zmin,zmax=-Dz/2+zoffset, Dz/2+zoffset

    xp = np.squeeze(Dx*nip.xx((1,Nx), freq='ftfreq'))
    z   = np.linspace(zmin, zmax, Nz)       # sampling of z-cordinates

    # sample grid
    mgz, mgx = np.meshgrid(z, xp)
    
    # compute the intial field 
    lambdaM = lambda0/n0
    k = 2*np.pi/(lambdaM)
    zR = np.pi * w0**2 / lambdaM
    Rz = mgz*(1 + (zR/mgz)**2)
    Rz[np.isnan(Rz)] = np.inf
    wz = w0*np.sqrt(1 + (mgz/zR)**2)
    xiz = np.arctan(mgz/zR)
    
    # propagate
    field = w0 / wz * np.exp(-(mgx/wz)**2) * np.exp(-1j*k*mgx**2/(2*Rz)) * np.exp(-1j*(k*mgz-xiz))
    return field

# settings for scaling/rotating object

Npixels = 100
scalefctor = .6
mWavelength = 405
mNA=0.5
w_LS = .5
pixelsize = mWavelength/mNA/4


# generate widefield PSF
obj_tmp = nip.image((Npixels,Npixels,Npixels))
param = nip.PSF_PARAMS()
param.NA = mNA
param.wavelength = mWavelength
obj_tmp.pixelsize=(pixelsize,pixelsize,pixelsize)
psf_wf = nip.psf(obj_tmp,param)




mGaussian=nip.image(nip.gaussian((Npixels,Npixels),sigma=2))
mGaussian.pixelsize=(pixelsize, pixelsize)
lsfm=nip.abssqr(nip.propagationStack(mGaussian,Npixels,distZ=pixelsize, psf_params = param))

Dz = (Npixels*pixelsize)
field_ds = nip.abssqr(propagate_gaussian(size=(Npixels,Npixels), sampling=(pixelsize,pixelsize), lambda0=mWavelength, n0=1, w0=w_LS, zoffset=0))
psf_illu = np.transpose(np.reshape(np.tile(field_ds,Npixels),[Npixels,Npixels,Npixels]), [2,1,0])


#plt.imshow(psf_wf[:,Npixels//2,]), plt.title("PSF Widefield"), plt.show()
#plt.imshow(psf_illu[:,Npixels//2,]), plt.title("PSF Light-sheet"), plt.show()
#plt.imshow(np.real(nip.ft(psf_wf[:,Npixels//2,]))), plt.title("PTF Widefield"), plt.show()


# export for debug
import tifffile as tif
tif.imsave("psf_illu.tif", np.array(psf_illu))
tif.imsave("psf_wf.tif", np.array(psf_wf))
#tif.imsave("obj.tif", np.array(obj))



#%%

all_OPT = []
all_LSFM = []

Nrot = 180
drot = 1
is_moving = 0


shift_x = Npixels//2
shift_y = Npixels//2 + is_moving * (0-Nrot//2)
shift_z = Npixels//2
rot_angle = np.pi/Nrot*0

# Using an existing stl file: https://sites.google.com/site/drjeffreywkrause/diatom-models
if(0):
    obj_3d = mesh.Mesh.from_file('Ceratium_Model.stl')
    obj_3d.rotate([0,1,0],rot_angle)
    unique_points = np.unique(obj_3d.points.reshape((obj_3d.points.shape[0] * 3, 3)), axis=0)
    
    # scale stl to fit into viewbox
    unique_points-=np.min(unique_points,0)
    unique_points/=np.max(unique_points,0)
    unique_points*=(Npixels*scalefctor-1)
    unique_points = np.int32(np.floor(unique_points))
    obj=np.zeros((Npixels,Npixels,Npixels))
    obj[unique_points[:,0],unique_points[:,1],unique_points[:,2]]=1

else:

    

obj=np.roll(obj,int(shift_x-(Npixels*scalefctor)//2),0)
obj=np.roll(obj,int(shift_y-(Npixels*scalefctor)//2),1)
obj=np.roll(obj,int(shift_z-(Npixels*scalefctor)//2),2)

for iangle in range(0,Nrot,drot):
    print(iangle)
    
    obj_rot= scipy.ndimage.interpolation.rotate(obj, iangle, mode='nearest', axes=(0, 1), reshape=False)
    
    mOPT = np.sum(nip.convolve(obj_rot, psf_wf),0)
    mLSFM = np.sum(nip.convolve(obj_rot*psf_illu, psf_wf),0)

    if(0):
        plt.subplot(131),plt.imshow(np.sum(obj_rot, 0), cmap='gray'), plt.title('Object')
        plt.subplot(132),plt.imshow(mOPT, cmap='gray'), plt.title('Proj. Tomo.')
        plt.subplot(133),plt.imshow(mLSFM, cmap='gray'),plt.title('Light-sheet')
        plt.savefig("plot_"+str(i)+".png")
        plt.show()
    
    # concte data
    all_OPT.append(mOPT)
    all_LSFM.append(mLSFM)


tif.imwrite("test_rot.tif", np.array(all_LSFM))
#%%

# Now reconstructing
# This is probably waaay to much effort, better would be a look up table that maps a single plane into 3D polar coordinates...


myvolume = np.zeros((*all_LSFM[0].shape,Npixels))
myvolume_norm = np.zeros((*all_LSFM[0].shape,Npixels))
myvolume_list = []
for iimage in range(0,len(all_LSFM),1):
    angle = iimage * drot 
    
    
    myvolume_tmp = 0*myvolume.copy()
    myvolume_tmp[:,:,Npixels//2]=all_LSFM[iimage]
    myvolume_tmp = scipy.ndimage.interpolation.rotate(myvolume_tmp, angle, mode='nearest', axes=(0, 2), reshape=False)
    myvolume+=myvolume_tmp
    myvolume_list.append(myvolume)
    if(1):
        # do the same for a slice of ones that can be used 
        myvolume_tmp = 0*myvolume.copy()
        myvolume_tmp[:,:,Npixels//2]=np.ones((Npixels,Npixels))
        myvolume_tmp = scipy.ndimage.interpolation.rotate(myvolume_tmp, angle, mode='nearest', axes=(0, 2), reshape=False)
        myvolume_norm+=myvolume_tmp    
    print(iimage)

#recon_vol_lsfm = np.transpose(myvolume,(2,0,1)) #
recon_vol_lsfm=myvolume/myvolume_norm
recon_vol_lsfm -= np.min(recon_vol_lsfm)
recon_vol_lsfm /= np.max(recon_vol_lsfm)
tif.imwrite("test.tif", np.float32(recon_vol_lsfm), imagej=True,
            metadata={'spacing': 1, 'unit': 'um', 'axes': 'ZYX'})

hf = h5py.File('data.h5', 'w')
hf.create_dataset('recon_vol_lsfm', data=np.float32(recon_vol_lsfm))
hf.create_dataset('all_LSFM', data=np.float32(all_LSFM))
hf.create_dataset('obj', data=np.float32(obj))
hf.close()
plt.imshow(recon_vol_lsfm[:,50,:]), plt.show()
plt.imshow(recon_vol_lsfm[50,]), plt.show()


def plusone(a):
    a=a+1
    return a

a=1
plusone(a)
print(a)