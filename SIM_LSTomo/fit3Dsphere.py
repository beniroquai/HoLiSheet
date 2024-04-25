#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 21:58:24 2022

@author: bene
"""

import numpy as np 
import NanoImagingPack as nip
import matplotlib.pyplot as plt

size3d=(256,256,256)

# ported from Rainer Heintzmanns MATLAB toolbox

def IterateCoefficients(subpix=40,CutoffK=10,TotalSize=(256,256,256),RollOff=20,Iterations=1500):
    '''
    # weight=IterateCoefficients(subpix,CutoffK,TotalSize,RollOff,Iterations) : Calculates a matrix of coefficients optimized for interpolation in Fourier space. To this aim an IFTA algorithm is used.
    # subpix : array of sizes for number of subdivisions between pixels
    # CutoffK : number of pixels to use (in both directions) as a kernel
    # TotalSize : Size of the multidimensional FFT array
    # RollOff : Sizes to use as border pixels
    # Iterations : Number of iterations to use during optimization
    #
    # Example:
    # map=IterateCoefficients(10,10,256,20,1000)  # Generates an interpolation map in 0.1 pixel distances for one fft direction
    # map=map=IterateCoefficients([10 10],[10 10],[256 256],[10 10],100)  # 2D interpolation kernels for a 2D FFT with only 10 pixels rolloff everywhere
    '''
    
    MySize=(subpix, TotalSize)
    if type(TotalSize)==int:
        allramps=nip.xx(MySize,freq='ftfreq') * nip.ramp(MySize,ramp_dim=0,placement='corner')/(subpix-1)
        toprocess=1
        fmask = np.abs(np.abs(nip.xx(MySize))>=CutoffK)
        rmask=(np.abs(nip.xx(MySize))<(np.floor(TotalSize/2)-RollOff))
    elif len(TotalSize.shape)==2:
        # TODO: untested
        print("NotImplemented")
        raise Exception
        '''
        allramps=nip.xx(MySize,freq='ftfreq') * nip.ramp(MySize,ramp_dim=2,placement='corner')/(subpix-1)
        toprocess=1
        fmask = np.abs(np.abs(nip.xx(MySize))>=CutoffK)
        rmask=(np.abs(nip.xx(MySize))<(np.floor(TotalSize/2)-RollOff))
        allramps=xx(MySize,'freq') .* ramp(MySize,3,'corner')/(subpix(1)-1)+yy(MySize,'freq') .* ramp(MySize,4,'corner')/(subpix(2)-1)
        toprocess=[1 1 0 0]
        fmask=~((abs(xx(MySize))<=CutoffK(1)) & (abs(yy(MySize))<=CutoffK(2)))
        rmask=(abs(xx(MySize))<(floor(TotalSize(1)/2)-RollOff(1))) & (abs(yy(MySize))<(floor(TotalSize(2)/2)-RollOff(2)))
        '''
    elif len(TotalSize.shape)==3:
        print("NotImplemented")
        raise Exception
        '''
        allramps=xx(MySize,'freq') .* ramp(MySize,5,'corner')/(subpix(1)-1)+yy(MySize,'freq') .* ramp(MySize,5,'corner')/(subpix(2)-1)+zz(TotalSize,'freq') .* ramp(MySize,6,'corner')/(subpix(3)-1)
        toprocess=[1 1 1 0 0 0]
        fmask=~((abs(xx(MySize))<=CutoffK(1)) & (abs(yy(MySize))<=CutoffK(2)) & (abs(zz(MySize))<=CutoffK(3)))
        rmask=(abs(xx(MySize))<(floor(TotalSize(1)/2)-RollOff(1))) & (abs(yy(MySize))<(floor(TotalSize(2)/2)-RollOff(2))) & (abs(zz(MySize))<(floor(TotalSize(3)/2)-RollOff(3)))
        '''
    else:
        
        print('IterateCoefficients only for a maximum of 3 dimensions')
        raise Exception
    
    perfsignal=np.exp(2*np.pi*1j*allramps)/np.sqrt(allramps.shape[0])  # Dim it down to account for Fourier space 
    
    rspace=nip.image(perfsignal.copy())
    meanPrevErr=0
    for n in range(Iterations):
        fspace = nip.ft(rspace, axes=toprocess)
        fspace[fmask]=0
        fspace=np.real(fspace)  # also force these coefficients to be real valued
        rspace = nip.ift(fspace, axes=toprocess)
        maxErr=np.max(np.abs(rspace[rmask]-perfsignal[rmask]))
        meanErr=np.mean(np.abs(rspace[rmask]-perfsignal[rmask]))
        #fprintf('Error iteration #d/#d is  : max #g mean #g\n',n,Iterations,maxErr,meanErr)
        print(meanErr)
        if np.abs((meanErr-meanPrevErr)/meanErr) < 1e-6:
            print('no more update\n')
            break
        meanPrevErr=meanErr
        rspace[rmask]=perfsignal[rmask]
    
    rspace = nip.ift(fspace, axes=toprocess)
    
    if type(TotalSize)==int:
        weight=fspace[:,TotalSize//2 - CutoffK:TotalSize//2+1 + CutoffK]
    elif len(TotalSize)==2:
        print("Not implemented")
        raise Exception
        #weight=fspace(np.floor(TotalSize(1)/2) - CutoffK(1):floor(TotalSize(1)/2) + CutoffK(1),floor(TotalSize(2)/2) - CutoffK(2):floor(TotalSize(2)/2) + CutoffK(2),:,:)
    elif len(TotalSize)==3:
        print("Not implemented")
        raise Exception
        #weight=fspace(floor(TotalSize(1)/2) - CutoffK(1):floor(TotalSize(1)/2) + CutoffK(1),floor(TotalSize(2)/2) - CutoffK(2):floor(TotalSize(2)/2) + CutoffK(2),floor(TotalSize(3)/2) - CutoffK(3):floor(TotalSize(3)/2) + CutoffK(3),:,:,:)
    
    normfac=np.tile(np.sqrt(np.sum(nip.abssqr(weight),axis=1))[:,np.newaxis],weight.shape[1])
    weight=weight / normfac  # this normalization is needed for energy conservation in projections
    weight=np.real(weight)  # again make sure that the datatype is real.

    return weight


def FillProjSpherePrepare(size3d,lambda0,pixelsizes,NA,imatrix,n=1.0,allowwraparound=1.0):
    '''
    # [indexList2D,fullIndex3D,factorList,amask]=FillProjSpherePrepare(size3d,lambda0,pixelsize,NA,imatrix,n,allowwraparound): Calculates preperatory data for the FillPrjSphere function
    # size3d: size of the 3d array to later fill in the data
    # lambda0: vakuum wavelength
    # pixelsize : size of a pixel in real space (same units as lambda0)
    # NA=sin(alpha) : numerical aperture for refractive index 1.0. Else use NA/n
    # imatrix : the interpolation matrix as a dipimage ranging from 0 to 1 subpixel difference and having the interpolation kernel in each row.
    # n: refractive index of the medium
    # allowwraparound: If True, interpolation can wrap around along Z
    #
    # Use "IterateCoefficients" to generate this matrix.
    # The output arguments (indexList2D,fullIndex3D,factorList) are all one-dimensional index lists into the pixels to access
    # indexList2D : indexes the pixels in the 2d projection (which is supposed to have the same XY size as the 3D array)
    # indexList 3D : indexes the position where the write the interpolation results
    # factorList : contains the interpolation weights for the appropriate pixel
    # amask : the 2D pupil mask used
    '''
    
    mid3d=np.floor(np.array(size3d)/2)
    size2d=size3d[0:2]
    k0=1/(lambda0/n)
    kxymax=k0*(NA/n)
    myr2=(nip.xx(size2d,freq='ftfreq')/pixelsizes[0])**2+(nip.yy(size2d,freq='ftfreq')/pixelsizes[1])**2
    kz=nip.image(myr2.shape)
    tmask=k0*k0 > myr2  # generate the fukll k-circle in 2D out of the existing k-space sphere
    if len(size3d) < 3:
        size3d[-1]=1
    
    
    kz[tmask]=-np.sqrt(k0**2-myr2[tmask])*pixelsizes[2]*size3d[2]  # calculate the hight map for this cap of the sphere
    # changed kz to be negative to correspond to an upward propagating wave.
    
    # amask=myr2 < kxymax*kxymax  # limit to the aperture chosen by the user
    # amask=myr2 < (kxymax*kxymax + k0*k0)/2  # limit to the aperture chosen by the user, but made intentionally larger to account for the edges of the pupil
    amask=np.sqrt(myr2) < (1*kxymax + 3*k0)/4  # limit to the aperture chosen by the user, but made intentionally larger to account for the edges of the pupil
    
    
    if (np.prod(imatrix.shape)==1) and (imatrix==1) and (len(size3d) < 3 or size3d[2]<=1):
        indexarray=nip.image(size2d)
        indexarray=np.reshape(np.arange(0,np.prod(indexarray.shape)), indexarray.shape)
        indexList2D=indexarray[amask]  # To be used for DipImage indexing
        fullIndex3D=indexList2D
        factorList=1
        return
    
    
    if np.sum((mid3d+k0+(imatrix.shape[0]+1)/2) > size3d[2]):
        print('Datavolume is not large enough along Z to accomodate all interpolation values')
        #return
    
    # Write the start z-index into the 2D mask
    
    kzarray=np.floor(kz)  # get the reference pixel position in Z at each XY position
    matrixindex = np.round((kz - kzarray)*(imatrix.shape[1]-1))  # this indicates the subpixel indexing. It denotes the row in the imatrix to use for interpolation
    
    indexarray=nip.image(size2d)
    indexarray=np.transpose(np.reshape(np.arange(0,np.prod(indexarray.shape)), indexarray.shape), (1,0))
    
    indexList2D=np.float32(indexarray.T[amask])  # To be used for DipImage indexing
    ## The code below creates an index list into the 3d array
    kzoffset=mid3d[2]-(imatrix.shape[1]-1)/2
    if np.max(kzarray+kzoffset) >= size3d[2] or np.min(kzarray+kzoffset)<0:
        if (0):
            error('Pupil kz larger than 3D data volume allows. PSF is undersampled')
        else:
            if allowwraparound:
                print('WARNING: Pupil kz larger than 3D data volume allows. Allowing for interpolation to wrap around and cause potential aliasing by Z pixels.')
                kzarray=np.mod(kzarray+kzoffset,size3d[2])
                kzoffset=0 # as it is already accounted for
            else:
                kzoffset=-np.min(kzarray[amask])
                fprintf('WARNING: Pupil kz larger than 3D data volume allows. PSF is undersampled. Shifting applying an offset ')
                if np.max(kzarray[amask])+kzoffset + imatri.shape(0) >= size3d[2]:
                    print('Pupil kz larger than 3D data volume allows. The Problem can be undersampling of the PSF or a too small Z data volume, which also needs to have space for the interpolators to fit completely.')
                    return
    print("something weird is with imatrix!!!!")
    indexlist3d=np.float32((kzarray[amask]+kzoffset) * (size3d[0]*size3d[1]) + indexList2D)   # This is the start index of each line to write in 3D
    
    fullIndex3D=np.transpose(np.tile(indexlist3d,(imatrix.shape[1],1))) + np.expand_dims(np.transpose(indexlist3d*0+1),-1) * (np.arange(0,imatrix.shape[1])-1)*(size3d[0]*size3d[1])
    #fullIndex3D=np.transpose(fullIndex3D)
    fullIndex3D=np.mod(fullIndex3D.flatten(),np.prod(size3d)) # To be used as DipImage indexing
    
    matrixindexlist=np.int32(matrixindex[amask])
    factorList= np.transpose(imatrix[matrixindexlist,:])

    return indexList2D,fullIndex3D,factorList,amask


#%%

def FillProjSphere(img,aproj,indexList2D,fullIndex3D,factorList):
    
    vallist=aproj.flatten()[np.int32(indexList2D)]
    vallist=np.reshape(vallist,(vallist.shape[0], 1));
    vallist=np.transpose(factorList, (1,0)) * vallist # outer product
    
    mshape = img.shape
    img = img.flatten()
    img[np.int32(fullIndex3D)]=img.flatten()[np.int32(fullIndex3D)]+vallist.flatten()
    img = np.reshape(img, mshape)
    
    print("ERROR!!!! Something is wrong with this command!!")
    #img = np.roll(img,3,0)
    return img



