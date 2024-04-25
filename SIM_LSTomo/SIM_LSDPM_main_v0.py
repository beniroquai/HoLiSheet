#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 11:29:21 2022

@author: bene
"""
import NanoImagingPack as nip
import numpy as np
import matplotlib.pyplot as plt
nip.setViewer("NIP_VIEW")

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

def propagate_gaussian_nip(size(128,128), mWavelength=405, w0=5):

    amp_in = nip.gaussian((size), sigma=w0)
    amp_in /= np.max(amp_in)
    ain = nip.image((amp_in+0j)*np.exp(1j*amp_in))
    
    ain.pixelsize=(pixelsize, pixelsize)
    PSFpara = nip.PSF_PARAMS()
    PSFpara.wavelength = mWavelength
    PSFpara.NA=1
    pupil = nip.ft(ain)
    propagated_minus = nip.ift2d(nip.propagatePupil(pupil, sizeZ=Nz/2, distZ=pixelsize/2, psf_params=PSFpara))
    propagated_plus = nip.ift2d(nip.propagatePupil(pupil, sizeZ=Nz/2, distZ=-pixelsize/2, psf_params=PSFpara))
    
    propagated = np.vstack((propagated_minus, propagated_plus))
    
    plt.imshow(np.abs(propagated[:,:,64])), plt.show()
    return propagated


def propagate_cylindrical_nip(size(128,128), pixelsize=1, mWavelength=405, w0=5, f=40):

    #%%
    amp_in = nip.gaussian((size), sigma=w0)
    amp_in /= np.max(amp_in)
    ain = nip.image((amp_in+0j)*np.exp(1j*amp_in))

    ain = nip.image(np.ones(ain.shape)+1j)
    ain *= nip.gaussf(1.*nip.rr(ain)<6,1)
    
    
    k0 = 2*np.pi/(mWavelength)
    
    kx = 1/pixelsize*nip.xx(size)#,freq='ftfreq')
    fi = -k0*(kx**2)/(2*(f*1e-6))
    ang_lens = np.exp(1j * fi)
    
    ain *= ang_lens
    
    
    ain.pixelsize=(pixelsize, pixelsize)
    PSFpara = nip.PSF_PARAMS()
    PSFpara.wavelength = mWavelength
    PSFpara.NA=1
    pupil = nip.ft(ain)
    propagated = nip.ift2d(nip.propagatePupil(pupil, sizeZ=Nz, distZ=pixelsize, psf_params=PSFpara, doDampPupil=True))
    
    plt.subplot(131),plt.imshow(np.abs(propagated[:,64,:]))
    plt.subplot(132),plt.imshow(np.abs(propagated[:,:,64]))
    plt.subplot(133),plt.imshow(np.abs(propagated[64,:,:])), plt.show()    
    #%%
    return propagated



# settings for scaling/rotating object

Npixels = 256
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

Nrot = 360
for i in range(Nrot):
    shift_x = Npixels//2
    shift_y = Npixels//2 +(i-Nrot//2)
    shift_z = Npixels//2
    rot_angle = np.pi/Nrot*i*2
    
    # Using an existing stl file: https://sites.google.com/site/drjeffreywkrause/diatom-models
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
    
    obj=np.roll(obj,int(shift_x-(Npixels*scalefctor)//2),0)
    obj=np.roll(obj,int(shift_y-(Npixels*scalefctor)//2),1)
    obj=np.roll(obj,int(shift_z-(Npixels*scalefctor)//2),2)

    
    mOPT = np.sum(nip.convolve(obj, psf_wf),0)
    mLSFM = np.sum(nip.convolve(obj*psf_illu, psf_wf),0)


    plt.subplot(131),plt.imshow(np.sum(obj, 0), cmap='gray'), plt.title('Object')
    plt.subplot(132),plt.imshow(mOPT, cmap='gray'), plt.title('Proj. Tomo.')
    plt.subplot(133),plt.imshow(mLSFM, cmap='gray'),plt.title('Light-sheet')
    plt.savefig("plot_"+str(i)+".png")
    plt.show()
    
    # concte data
    all_OPT.append(mOPT)
    all_LSFM.append(mLSFM)



