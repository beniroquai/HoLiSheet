#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 22:30:39 2022

@author: bene
"""

import pyclesperanto_prototype as cle
from skimage.io import imread
from skimage.transform import rotate
import numpy as np
from napari_animated_gif_io import napari_write_image as imsave_animated_gif


cle.select_device("RTX")

#%%

%gui qt
import napari

# Create an empty viewer
viewer = napari.Viewer()

#%%
folder = '/Users/bene/Downloads/BioImageAnalysisNotebooks/docs/18b_reconstruction/'
filename = 'interesting_rot_bead_LS_4_rec_LightSheet-1.tif'
filename = "interesting_rot_algaespiralstable_rec_LightSheet-2.tiff"
filename = "interesting_rot_multialgaeLS_rec_LightSheet180deg2.tiff"

original = cle.push(imread(folder+filename))
raw_angles = np.linspace(0,180,original.shape[0]).astype(np.uint16)

raw_stack = {}
index = 0
for angle in raw_angles:
    raw_stack[angle] = original[index,]
    index+=1
    
# The algorithm shown below will only work if there are no zeros in the raw data.
# Furthermore, the animations only work if the raw data is in range [0, 255]
original = original * 254.0 / original.max() + 1
original.min(), original.max()

def acquire_slice_at_angle(original, angle):
    result = np.expand_dims(original[angle],0)
    return result

example = acquire_slice_at_angle(original, 1)
cle.imshow(example)
example = acquire_slice_at_angle(original, 15)
cle.imshow(example)


# create space for reconstructed image
raw_stack = original
raw_angles = np.linspace(0,180,raw_stack.shape[0])
first_image = raw_stack[0,]
target = cle.create((first_image.shape[1], first_image.shape[0], first_image.shape[1]))
cle.set(target, 0)
print(target.shape)

#demo = np.asarray(list(raw_stack)).astype(np.uint8)
#imsave_animated_gif("demo_1.gif", demo, meta={})

print(raw_angles)
print(raw_stack.shape)
viewer.add_image(cle.pull(raw_stack))


import matplotlib.pyplot as plt

def write_image_at_angle(source, target, angle, translation):
    # rotate a single image around y-axis and store it
    # in the target 3D stack
    temp = cle.create_like(target)
        
    transform = cle.AffineTransform3D()
    transform.center(source.shape)
    transform.rotate_around_y_axis(angle)
    transform.center(target.shape, undo=True)
    transform.translate(0,translation,0) # xyz
    cle.affine_transform(source, temp, transform=transform)
    
    cle.maximum_images(target, temp, target)
    #cle.add_images_weighted(target, temp, target)
    
    return target

def create_stack_from_angular_images(raw_stack, raw_angles, target):
    # Write images and multiple angles in to target stack
    for iimage in range(len(raw_angles)):
        translation = iimage
        angle = int(raw_angles[iimage])
        image = acquire_slice_at_angle(original, iimage) #raw_stack[iimage]
        write_image_at_angle(image, target, angle, translation)
    
    return target 

create_stack_from_angular_images(raw_stack, raw_angles, target)

import tifffile as tif
tif.imwrite("test.tif", np.array(target))
viewer.add_image(cle.pull(target))