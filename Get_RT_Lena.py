#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 17:12:36 2018

@author: Zou-Williams, L.
"""


import mne
import numpy as np
import os.path as op
import glob
import csv
mne.set_log_level('ERROR')



def rt(set_file):
    raw  = mne.io.read_raw_brainvision(set_file,eog=['SO1', 'IO1', 'LO1', 'LO2'],montage='standard_1020',preload=True)
    events = mne.find_events(raw)
    t = events[:,0].tolist()
    r = events[:,2].tolist()
    rt = zip(t,r)
    with open('rt/' + set_file + '.csv', "w") as f:
        writer = csv.writer(f)
        for row in rt:
            writer.writerow(row)    
