#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import NanoImagingPack as nip
import tifffile as tif
import model as mus
from stl import mesh
import stl

# Optionally, tweak styles.
mpl.rc('figure',  figsize=(9, 6))
mpl.rc('image', cmap='gray')
#plt.switch_backend('agg')
                                                    
''' Define parameters '''
is_padding = False # better don't do it, some normalization is probably incorrect
is_display = True
is_optimization = False 
is_optimization_psf = False
is_flip = False
is_measurement = False

'''Choose between Born (BORN) or BPM (BPM)'''
#psf_model =  'BORN' # MultiSlice
#psf_model =  '3QDPC' # MultiSlice
psf_model =  'BPM' # 1st Born


''' Create the Model'''
Nx=160
Ny=160
Nz=128
lambda0=405e-9
nEmbb=1.33
dn_glass=.2
dn_plankton=.05
dn=0.1
NAo=.5
NAc=.003


dxy=lambda0/NAo/4
dz=dxy*2




# create capillary
r_capillary=50
obj_real =  (1*(nip.rr((Nz,Nx))<r_capillary)-1*(nip.rr((Nz,Nx))<r_capillary-2))*dn_glass+nEmbb
obj_real = np.transpose(np.repeat(np.expand_dims(obj_real,0),Ny,0),[1,2,0])

myallsamples=[]

for iangle in range(1):
    muscat = mus.MuScatModel(size=(Nz,Nx,Ny), dxy=dxy, dz=dz, lambda0=lambda0, 
                         nEmbb=nEmbb, dn=dn, NAo=NAo, NAc=NAc)

    if 1:
        # Using an existing stl file: https://sites.google.com/site/drjeffreywkrause/diatom-models
        rot_angle=iangle/20*np.pi*2
        scalefctor=.4
        obj_3d = mesh.Mesh.from_file('Ceratium_Model.stl')
        obj_3d.rotate([0,0,1],rot_angle)
        unique_points = np.unique(obj_3d.points.reshape((obj_3d.points.shape[0] * 3, 3)), axis=0)
        
        # scale stl to fit into viewbox
        unique_points-=np.min(unique_points,0)
        unique_points/=np.max(unique_points,0)
        unique_points*=(Nx*scalefctor-1)
        unique_points = np.int32(np.floor(unique_points))
        obj=np.zeros((Nz,Nx,Ny))
        obj[unique_points[:,0],unique_points[:,1],unique_points[:,2]]=1
        shift_x = Nx//2
        shift_y = Ny//2 #+(i-Nrot//2)
        shift_z = Nz//2
        obj=np.roll(obj,int(shift_x-(Nx*scalefctor)//2),1)
        obj=np.roll(obj,int(shift_y-(Nx*scalefctor)//2),2)
        obj=np.roll(obj,int(shift_z-(Nx*scalefctor)//2),0)
        obj_real[obj>0]=dn_plankton+nEmbb
        
    else:
        ## add plankton
        obj_tmp = np.transpose(np.repeat(np.expand_dims((nip.rr((Nz,Nx))<23),0),Ny,0),[1,2,0])
        obj_real[obj_tmp*np.random.randn(*muscat.size)>2]=dn_plankton+nEmbb


    obj_real=nip.gaussf(obj_real,1)
    obj_absorption = obj_real*0
    obj = obj_real + 1j*obj_absorption

        
    # This function is a wrapper to compute the Born fwd-model (convolution)
    muscat.computesys(obj=obj, is_tomo = False, is_compute_psf=psf_model, is_mictype='BF')
    # Create Model Instance
    my3dfield=muscat.computemodel(is_allslices=True)

    if 0:
        plt.subplot(231),plt.imshow(np.abs(np.squeeze(np.array(muscat.allslices))[:,:,Ny//2])), plt.colorbar()
        plt.subplot(232),plt.imshow(np.angle(np.squeeze(np.array(muscat.allslices))[:,:,Ny//2])), plt.colorbar()
        plt.subplot(233),plt.imshow(np.squeeze(np.abs(obj[:,:,Ny//2]))), plt.colorbar()
        plt.subplot(236),plt.imshow(np.squeeze(np.mean(np.abs(obj),0))), plt.colorbar()

        plt.subplot(234),plt.imshow(np.squeeze(np.abs(my3dfield[:,:,Ny//2]))), plt.colorbar()
        plt.subplot(235),plt.imshow(np.squeeze(np.abs(my3dfield[Nz//2,:,:]))), plt.colorbar(), plt.show()
        plt.subplot(236),plt.imshow(np.squeeze(np.mean(np.abs(obj),0))), plt.colorbar()


    myallsamples.append(np.abs(my3dfield[Nz//2,:,:]))


tif.imsave("tmp.tif",np.array(myallsamples))