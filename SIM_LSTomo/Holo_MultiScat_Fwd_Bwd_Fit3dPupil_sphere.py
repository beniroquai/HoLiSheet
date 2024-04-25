#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import NanoImagingPack as nip
import tifffile as tif
import model as mus
import scipy
from fit3Dsphere import *

# Optionally, tweak styles.
mpl.rc('figure',  figsize=(9, 6))
mpl.rc('image', cmap='gray')
#plt.switch_backend('agg')


#%%

#%gui qt
import napari
# Create an empty viewer
viewer = napari.Viewer()

#%%                                                    
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
Nx=100
Ny=100
Nz=100
lambda0=405e-9
nEmbb=1.33
dn_glass=.2
dn_plankton=.05
dn=0.1
NAo=.5
NAc=.003


dxy=lambda0/NAo/4
dz=dxy*2


#%% FWD Model => Simulate a rotating sample imaged coherently

# create capillary
r_capillary=50
obj_real =  (1*(nip.rr((Nz,Nx))<r_capillary)-1*(nip.rr((Nz,Nx))<r_capillary-2))*dn_glass+nEmbb
obj_real = np.transpose(np.repeat(np.expand_dims(obj_real,0),Ny,0),[1,2,0])


# select the sample 
IS_STL = False


#%%
# rotate only the computed that has been computed once
allIntensityImages=[]
allSamples=[]

# sample acquisition parameters
dAngle = 360
Nrot = 4
    
allImages = np.zeros((Nx,Ny,Nz))+0j
allImagesRecon = np.zeros((Nx,Ny,Nz))+0j

for iangle in range(Nrot):
    rot_angle=iangle/Nrot*dAngle
    print(rot_angle)


    if 1:
        # 1. Convert the STL-plankton to a 3D volume
    
        scalefctor=.4
        if IS_STL:
            from stl import mesh
    
            # Using an existing stl file: https://sites.google.com/site/drjeffreywkrause/diatom-models
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
            
        elif 0:
            # sideload a volume with proper dimensions
            obj = np.load("obj.npy")
        else:
            # create a sphere
            obj = dn_plankton*(nip.rr((Nx,Ny,Nz))<10)+nEmbb
    
        
        # let the sample "flow" along the fluidics direction?
        if 0:            
            shift_x = Nx//2
            shift_y = Ny//2+(iangle-Nrot//2)
            shift_z = Nz//2
            obj=np.roll(obj,int(shift_x-(Nx*scalefctor)//2),1)
            obj=np.roll(obj,int(shift_y-(Nx*scalefctor)//2),2)
            obj=np.roll(obj,int(shift_z-(Nx*scalefctor)//2),0)
        #obj_real[obj>0]=dn_plankton+nEmbb
        
    else:
        ## add a simple sphere with random background
        obj_tmp = np.transpose(np.repeat(np.expand_dims((nip.rr((Nz,Nx))<23),0),Ny,0),[1,2,0])
        obj_real[obj_tmp*np.random.randn(*muscat.size)>2]=dn_plankton+nEmbb
    
    # rotate sample
    obj_rot = scipy.ndimage.interpolation.rotate(obj, rot_angle, mode='nearest', axes=(0, 1), reshape=False)
    allSamples.append(obj_rot)
    
    # initiliaze the Muscat model
    muscat = mus.MuScatModel(size=(Nz,Nx,Ny), dxy=dxy, dz=dz, lambda0=lambda0, 
                             nEmbb=nEmbb, dn=dn, NAo=NAo, NAc=NAc)
    # simulate the acquisition process (multiple scatterin through sample, then refocus to cameraplane)
    muscat.computesys(obj=obj_rot, is_tomo = False, is_compute_psf=psf_model, is_mictype='AIP')


    # Create Model Instance
    my3dfield_rot = muscat.computemodel(is_allslices=True)
    
    # compute 2d image
    bfImage = nip.abssqr(my3dfield_rot[muscat.size[0]//2,])
    
    # subtract background 
    dGauss = 10
    aipImage = (bfImage-nip.gaussf(bfImage,dGauss))/nip.gaussf(bfImage,dGauss)
    #plt.imshow(np.abs(nip.ft(aipImage))**.15)

    #%% Shift spectrum back to center
    reconResult = np.abs(np.roll(nip.ft(aipImage), -(muscat.size[1]//4-4))*np.squeeze(muscat.Po))
    plt.imshow(np.abs(reconResult)**.15), plt.show()
    
    #plt.imshow(np.real(nip.ift(reconResult)))
               
    #plt.imshow(np.abs(nip.ift(reconResult))), plt.colorbar()   , plt.show()            

        #%%
#    viewer.add_image(np.abs(my3dfield_rot))

    # delete muscat
    del muscat

    if 0:
        plt.subplot(231),plt.imshow(np.abs(np.squeeze(np.array(muscat.allslices))[:,:,Ny//2])), plt.colorbar()
        plt.subplot(232),plt.imshow(np.angle(np.squeeze(np.array(muscat.allslices))[:,:,Ny//2])), plt.colorbar()
        plt.subplot(233),plt.imshow(np.squeeze(np.abs(obj[:,:,Ny//2]))), plt.colorbar()
        plt.subplot(236),plt.imshow(np.squeeze(np.mean(np.abs(obj),0))), plt.colorbar()

        plt.subplot(234),plt.imshow(np.squeeze(np.abs(my3dfield[:,:,Ny//2]))), plt.colorbar()
        plt.subplot(235),plt.imshow(np.squeeze(np.abs(my3dfield[Nz//2,:,:]))), plt.colorbar(), plt.show()
        plt.subplot(236),plt.imshow(np.squeeze(np.mean(np.abs(obj),0))), plt.colorbar()

    
    # collect image from center focus slice (assuming middle section represents focus)
    allIntensityImages.append(bfImage)
    
    # backrotate sample according to sample rotation angle
    allImages+=my3dfield_rot# scipy.ndimage.interpolation.rotate(my3dfield_rot, -rot_angle, mode='nearest', axes=(0, 1), reshape=False)


    # "reconstruct" the data or place them in FT space => backrotate
    allImagesRecon += scipy.ndimage.interpolation.rotate(my3dfield_rot, -rot_angle, mode='nearest', axes=(0, 1), reshape=False)
    

viewer.add_image(np.array(allSamples))
viewer.add_image(np.array(allIntensityImages))

#tif.imsave("tmp.tif",np.array(allIntensityImages))
viewer.add_image(np.log(1+np.abs(nip.ft(allImages))))
viewer.add_image(np.abs(allImagesRecon))
viewer.add_image(np.log(1+np.abs(nip.ft(allImagesRecon))))


#%%PAUSE
asldkfj
#%% Setup parametres
size3d=(Nx,Ny,Nz)
img=np.zeros(size3d) #newim(size3d,'scomplex') 
aproj = np.ones((size3d[0],size3d[1]))
NAo=NAo
lambda0=lambda0*1e6
pixelsize=lambda0/NAo/4 #muscat.dxy
pixelsizes=(pixelsize,pixelsize,pixelsize)
n=1
allowwraparound=True


#  40 subpixel subdivisions
# 2*10+1 kernelsize
# 20 pixel bordersize in all directions
# 1500 iterations
#%% Compute LUT

NA=.5
lambda0=0.450
pixelsize=lambda0/NA/4


imatrix=IterateCoefficients(subpix=40,CutoffK=10,TotalSize=size3d[-1],RollOff=20,Iterations=1000)
#%% MAP LUT in 3d space
indexList2D,fullIndex3D,factorList,amask=FillProjSpherePrepare(size3d,lambda0,pixelsizes,NAo,imatrix,n=n,allowwraparound=allowwraparound)

#%% rotate image in 3D around the rotational axis

allAngles = []
allFT = np.zeros(size3d)+0j

for iangle in range(Nrot):
    rot_angle=iangle/Nrot*dAngle
    print(rot_angle)
    aproj = nip.ft(allIntensityImages[iangle])
    #% Map image to LUT in 3D 
    img = nip.image(size3d)+0j
    img = FillProjSphere(img,aproj,np.int32(indexList2D),fullIndex3D,factorList) # Just a circular wave
    img = np.roll(img, size3d[0]//2-3, 0)
    img = scipy.ndimage.interpolation.rotate(img, rot_angle, mode='nearest', axes=(0, 1), reshape=False)
    allFT += img
    
    #viewer.add_image(np.log(1+np.abs(img)))
    plt.imshow(np.log(1+np.abs(img[:,size3d[1]//2,:]))), plt.show()
    del img
viewer.add_image(np.abs(allFT))
#%%
img_Ft = nip.ft(allFT) 

viewer.add_image(np.abs(img_Ft))