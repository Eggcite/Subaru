#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 09:25:26 2018

@author: Zou-Williams, L.
"""

import glob 
import mne
import numpy 
import pandas as pd
import warnings
import philistine as phil
import csv
import os

warnings.filterwarnings("ignore", category=DeprecationWarning)
      

    
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
        R=fig.axes[0].lines[0]
        W=fig.axes[0].lines[1]
        rx=R.get_xdata()
        ry=R.get_ydata()
        wy=W.get_ydata()
        x=rx.tolist()
        ry=ry.tolist()
        wy=wy.tolist()
        x.insert(0,'x')
        ry.insert(0,'ry')
        wy.insert(0,'wy')
        d = zip(x, ry, wy)
        
        with open('plotdata/' + group + '/' + name + '.csv', "w") as f:
            writer = csv.writer(f)
            for row in d:
                writer.writerow(row)    


filenames = glob.glob('epochs/*' + '-epo.fif.gz')

# getting data
get_data(filenames,event_id_nitem)

# for plotting 
# creating group list for high/low English rating groups
high_list=['11','01','02','04','12','16','18','21','22','23','24','28']
low_list=['03','07','08','09','15','17','19','20','26','27','29','30','31','32']

''' choose a possible comparison
plausibility = ("SAR","SAW")
AV_V = ("AVR/V","AVW/V")
phrase=("PSVR","PSVW")
IV_V = ("IVR/V", "IVW/V")
AV_N1=("AVR/N1","AVW/N1")
AV_N2=("AVR/N2","AVW/N2")
'''
#choose a comparison
compare=AV_V

#group can be 'high', 'low', 'all'
extract_plot(filenames=filenames,compare=compare,group='low')




    
    



