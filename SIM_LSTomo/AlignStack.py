#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 08:22:40 2022

@author: bene
"""


import tifffile as tif
import numpy as np
import NanoImagingPack as nip
import matplotlib.pyplot as plt

filepath = '/Users/bene/Downloads/combined_lightsheet_rec_LightSheet-1_fullrot_composite-1.tif'

myimage_stack = tif.imread(filepath)

#%%

dx = 50
dy = 245
myimage_stack = tif.imread(filepath)
myimage_stack_shifted = []
for iimage in range(myimage_stack.shape[0]):
    myimage_raw = myimage_stack[iimage,]
    myimage_channel1_Shift = np.roll(myimage_raw[0,],dy, axis=0)
    myimage_channel1_Shift = np.roll(myimage_channel1_Shift,dx, axis=1)
    
  t?    myimage = (myimage_channel1_Shift, 
               myimage_raw[1,])
    myimage =  np.transpose(myimage, (1,2,0))
    myimage_stack_shifted.append(myimage)
    
    print(iimage)
    
tif.imwrite("test.tif", np.array(myimage_stack_shifted)[:,:,0])