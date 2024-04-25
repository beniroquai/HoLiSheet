import cv2
import sys
import tifffile as tif
#!pip install opencv-contrib-python
#!pip install opencv-python
(major_ver, minor_ver, subminor_ver) = (cv2.__version__).split('.')
import numpy as np
import matplotlib.pyplot as plt
import NanoImagingPack as nip
# Set up tracker.
# Instead of CSRT, you can also use



%gui qt
import napari

# Create an empty viewer
viewer = napari.Viewer()

# Read 
filename = 'SingleParticleInFocus4_rec_Holography-3.tif'
filename = '10h02m44s_rec_Holography-1.tif'
frames = tif.imread(filename)
'''
#%%
tracker_types = ['BOOSTING', 'MIL','KCF', 'TLD', 'MEDIANFLOW', 'GOTURN', 'MOSSE', 'CSRT']
tracker_type = tracker_types[1]
 
if int(minor_ver) < 3:
    tracker = cv2.Tracker_create(tracker_type)
else:
    if tracker_type == 'BOOSTING':
        tracker = cv2.TrackerBoosting_create()
    elif tracker_type == 'MIL':
        tracker = cv2.TrackerMIL_create()
    elif tracker_type == 'KCF':
        tracker = cv2.TrackerKCF_create()
    elif tracker_type == 'TLD':
        tracker = cv2.TrackerTLD_create()
    elif tracker_type == 'MEDIANFLOW':
        tracker = cv2.TrackerMedianFlow_create()
    elif tracker_type == 'GOTURN':
         tracker = cv2.TrackerGOTURN_create()
    elif tracker_type == 'MOSSE':
        tracker = cv2.TrackerMOSSE_create()
    elif tracker_type == "CSRT":
        tracker = cv2.TrackerCSRT_create()
 

#plt.imshow((1+np.abs(nip.ft(frame/nip.gaussf(frame,11))))**.2)
# Read first frame.
frame = frames[0]

 
# Define an initial bounding box
bbox = (287, 23, 86, 320)
 
# Uncomment the line below to select a different bounding box
bbox = cv2.selectROI(frame, False)
 
#% Initialize tracker with first frame and bounding box
ok = tracker.init(frame, bbox)
 
# 

allTracks = []
#%
trackedObject = []
for iframe in range(frames.shape[0]):
    print(iframe)
    # Read a new frame
    frame =  frames[iframe]
     
    # Update tracker
    ok, bbox = tracker.update(frame)
 
    trackedObject.append(frame[:,bbox[0]:bbox[0] + bbox[2]])
    
    # Draw bounding box
    if ok:
        # Tracking success
        p1 = (int(bbox[0]), int(bbox[1]))
        p2 = (int(bbox[0] + bbox[2]), int(bbox[1] + bbox[3]))
        frameBoxed = frame.copy()
        cv2.rectangle(frameBoxed, p1, p2, (255,0,0), 2, 1)
    else :
        # Tracking failure
        pass
 
    cv2.imshow("Tracking", frameBoxed)
 
    # Exit if ESC pressed
    if cv2.waitKey(1) & 0xFF == ord('q'): # if press SPACE bar
        break

cv2.destroyAllWindows()

#tif.imwrite("allTracks.tif", np.array(allTracks))

tif.imwrite("trackedObject.tif", np.array(trackedObject))
'''

#%%
from skimage import filters
from skimage.registration import phase_cross_correlation
import math 

# shift back the tube to the center -> anti wiggeling 
allShifts = []
allBackshifted = []
for iframe in range(len(frames)-1):
    src = frames[iframe][0:300,:]
    dst = cv2.Canny(src, 0, 50)
    lines = cv2.HoughLines(dst, 1, np.pi / 180, 150, None, 0, 0)
    
    if lines is not None:
        ypos = np.squeeze(lines[0])[0]
    else:
        ypos = 0
    allShifts.append(ypos)
        
    frame = cv2.resize(frames[iframe], dsize=None, dst=None, fx=.25, fy=.25)
    allBackshifted.append(np.roll(frame,-int(model(iframe)/4) ,0))

    
timetrace = np.array(allShifts)-np.min(np.array(allShifts))
timefct = np.arange(0,timetrace.shape[0])
model = np.poly1d(np.polyfit(timefct, timetrace, 2))

plt.plot(model(timefct))
plt.plot(timetrace)
plt.show()



    
viewer.add_image(np.array(allBackshifted))

#%%
