#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 11:45:28 2017

@author: georgia
"""

#%% Import all the things you need
from __future__ import division
from __future__ import print_function
from builtins import zip
from builtins import str
from builtins import map
from builtins import range
from past.utils import old_div
import cv2
try:
    cv2.setNumThreads(1)
except:
    print('Open CV is naturally single threaded')

try:
    if __IPYTHON__:
        get_ipython().magic('load_ext autoreload')
        get_ipython().magic('autoreload 2')
except NameError:
    print('Not IPYTHON')
    pass

import caiman as cm
import numpy as np
import os
import glob
import time
import matplotlib.pyplot as plt
import psutil
import sys
from ipyparallel import Client
from skimage.external.tifffile import TiffFile
import scipy

from caiman.motion_correction import MotionCorrect, tile_and_correct, motion_correction_piecewise
from caiman.utils.utils import download_demo

#%% First setup some parameters

# dataset dependent parameters
os.chdir('/home/georgia/ToBeMotionCorrected/SC01-1/171027')

fr = 30                             # imaging rate in frames per second
decay_time = 0.4                    # length of a typical transient in seconds

# motion correction parameters
niter_rig = 1               # number of iterations for rigid motion correction
max_shifts = (50, 50)         # maximum allow rigid shift
splits_rig = 50             # for parallelization split the movies in  num_splits chuncks across time
strides = (128, 128)          # start a new patch for pw-rigid motion correction every x pixels
overlaps = (32, 32)         # overlap between pathes (size of patch strides+overlaps)
splits_els = 50             # for parallelization split the movies in  num_splits chuncks across time
upsample_factor_grid = 50    # upsample factor to avoid smearing when merging patches
max_deviation_rigid = 6     # maximum deviation allowed for patch with respect to rigid shifts
dview = Client #need to add the ipyparallel object??

#%% select video 
fname = os.path.join(os.getcwd(),'Substackmini.tif')  # filename to be processed

#%% start a cluster for parallel processing
c, dview, n_processes = cm.cluster.setup_cluster(
    backend='local', n_processes=None, single_thread=False)

#%%% MOTION CORRECTION
# first we create a motion correction object with the parameters specified
min_mov = cm.load(fname[0], subindices=range(200)).min() 
        # this will be subtracted from the movie to make it non-negative 

mc = MotionCorrect(fname[0], min_mov,
                   dview=dview, max_shifts=max_shifts, niter_rig=niter_rig,
                   splits_rig=splits_rig, 
                   strides= strides, overlaps= overlaps, splits_els=splits_els,
                   upsample_factor_grid=upsample_factor_grid,
                   max_deviation_rigid=max_deviation_rigid, 
                   shifts_opencv = True, nonneg_movie = True)
# note that the file is not loaded in memory




