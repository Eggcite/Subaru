#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 17:09:51 2018

@author: Zou-Williams
"""


from multiprocessing import Pool
import glob 
import mne
import numpy 
import pandas as pd
import warnings
import philistine as phil
import csv
import os
import argparse
import os.path as op
import re
warnings.filterwarnings("ignore", category=DeprecationWarning)


def grand_epochs(subjects,event_id):
    raw = mne.io.read_raw_fif(subjects)
    events = mne.find_events(raw)
    d=events[:,2].tolist()
    
    for x in d:
        pick=[230,231,232,233,234,235,236,237,238,239,240,241]
        if x in pick:
            if d[d.index(x)+1]<8:
                d[d.index(x)]=x*100
                continue
            else:
                for n in range(1,8):
                    if d[d.index(x)+n+1]==n:
                        d[d.index(x)+n+1]=x*100+n
                    else:
                        continue
                d[d.index(x)]=x*100
            continue        
            
    events[:,2]=d
    tmin = -0.2
    tmax = 1.2
    baseline = None
   # raw.info['bads'] = self.bads
    eog_events = mne.preprocessing.find_eog_events(raw)
    n_blinks = len(eog_events)
    onset = eog_events[:, 0] / raw.info['sfreq'] - 0.25
    duration = numpy.repeat(0.5, n_blinks)
    raw.annotations = mne.Annotations(onset, duration, ['bad blink'] * n_blinks)
                             # orig_time=raw.info['meas_date'])
    picks = mne.pick_types(raw.info, meg=False, eeg=True, eog=True, stim=False, exclude='bads')
    reject = dict(eeg=150e-6)
    flat = dict(eeg=5e-6)
    epochs = mne.Epochs(raw, events, event_id=event_id,
                        tmin=tmin, tmax=tmax,
                        proj=False, picks=picks,
                        baseline=baseline,
                        detrend=0,  # DC offset
                        reject_by_annotation=True,
                        flat=flat,
                        reject=reject,
                        preload=True)
    bad_epoch_mask = phil.mne.abs_threshold(epochs, 75e-6)
    epochs.drop(bad_epoch_mask, reason="absolute threshold")
    try:
        epochs.save('delay/epoched_grand/CP' + k[0] + '-epo.fif.gz')
    except IndexError:
        print(k[0], 'skipped for no epochs')


subjects = glob.glob('delay/processed/CP*_raw.fif.gz')

pool = Pool()
pool.map(grand_epochs, subjects,event_id)


def extract_plot(filenames, compare, outfolder, group, event_id):
    sleep_group=list()
    wake_group=list()
    for line in filenames:
        k=line[2:4].split()
        if k[0] in sleep_list:
            sleep_group.append(line)
        elif k[0] in wake_list:
            wake_group.append(line)
        if group == 'sleep':
            filenames = sleep_group
        elif group == 'wake':
            filenames = wake_group
        elif group == 'all':
            filenames = filenames
 #       print(filenames)
    
    ss_avg = dict()  # new empty dict
    ss_avg_by_cond = {e: list() for e in event_id}
    ss_avg_all = list()        
    try:
        for filename in filenames:
            epochs = mne.read_epochs(filename)
            evokeds = {str(cond): epochs[str(cond)].average() for cond in epochs.event_id}
            ss_avg[filename] = evokeds
            for cond in event_id:
                ss_avg_by_cond[cond].append(ss_avg[filename][cond])
                ss_avg_all.append(ss_avg[filename][cond])
    except (IOError, ValueError) as e:
        print(filename, 'skipped for no epochs')
        
    grand_avs = {cond:mne.grand_average(ss_avg_by_cond[cond]) for cond in event_id}
    grand_all=mne.grand_average(ss_avg_all)
    com = {k:grand_avs[k] for k in compare}
    
    for name in grand_all.ch_names:
        name_int=grand_all.ch_names.index(name)
        fig= mne.viz.plot_compare_evokeds((com),picks=[name_int],invert_y=True, show= False)
        fig.savefig('figures/' + zone + '-mon' + '_' + 'g-avg.pdf')
        
        A=fig.axes[0].lines[0]
        U=fig.axes[0].lines[1]
        Ax=A.get_xdata()
        Ay=A.get_ydata()
        Uy=U.get_ydata()
        x=Ax.tolist()
        Ay=Ay.tolist()
        Uy=Uy.tolist()
        x.insert(0,'x')
        Ay.insert(0,'Ay')
        Uy.insert(0,'Uy')
        d = zip(x, Ay, Uy)
         
        with open('Baseline/plotdata/' + outfolder + '/' + group + '/' + name + '.csv', "w") as f:
            writer = csv.writer(f)
            for row in d:
                writer.writerow(row)    

filenames = glob.glob('delay/epoched_grand/*' + '-epo.fif.gz')
filenames = glob.glob('*' + '-epo.fif.gz')

# for plotting 
# creating group list for high/low English rating groups
sleep_list=['02','05','08','09','11','12','13','16','18','20','22','25','26','29','35','37','38']
wake_list=['01','03','04','06','07','10','14','15','17','19','21','23','24','30','31','33','36']


#define comparisons
AU_V = ("23003","23103")

G_C = ('AbaUV/C','UbeAV/C')
G_N = ('AbaUV/U','UbeAV/A')
W_N = ('AbeUV/U', 'UbaAV/A')
W_C = ('AbeUV/C', 'UbaAV/C')
CV_C = ('AbaUV/C', 'UbeAV/C','AbeUV/C','UbaAV/C')
CV_N = ('AbaUV/U','UbeAV/A','AbeUV/U','UbaAV/A')
Cv = ("AbaUV/ba", "UbeAV/be")
N = ("AbaUV/A", "UbeAV/U")
V = ("AbaVU/V", "UbeVA/V")

An=('23203','23403')

order=('23204','23604')


event_id = {}
f = open('event_id_grand.txt', 'r')
for line in f:
    data=line.split()
    key, value = data[0], data[0]
    event_id[key] = int(value)
print(event_id)
   
#group can be 'high', 'low', 'all'
extract_plot(filenames=filenames,compare=An, outfolder='V', group='all', event_id=event_id)




extract_plot(filenames=filenames,compare=, outfolder='V', group='all', event_id=event_id)
