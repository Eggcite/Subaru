"""

Created on Circa August 2017
@Author: Zou-Williams, L. & Chan, R.W.

"""
#Packages used:
import mne
import glob
import numpy as np
import os.path as op
import pandas as pd
import statsmodels
import csv


mne.set_log_level('ERROR')

# Rolling definitions for functions for Russell's sanity and anyone learning -
# def abs_threshold is to drop epochs based on absolute voltages

# Here you can rely on baselining - but no?
def ss_epo(set_file, event_id, windows):
    # reading files
    raw = mne.io.read_raw_fif(set_file)

    # finding triggered events
    events = mne.find_events(raw)
    tmin = -0.2
    tmax = 0.8
    baseline = None

    # if set_file == 'MenTrainP33S1.vhdr-raw.fif':
    if set_file == 'MenTrainP33S1.vhdr-raw.fif':
       raw.info['bads'] = ['SO1']

    picks = mne.pick_types(raw.info, meg=False, eeg=True, eog=True,
                           stim=False, exclude='bads')

    # reject blinks and movement
    reject = dict(eeg=150e-6)
    flat = dict(eeg=5e-6)


    d = events[:, 2].tolist()
    d[0:200]

    for x in d:
        # for Visual-stimuli-onset 
        if x == 241:
            y = d[d.index(x) + 2]
            n = x * 1000 + y
            d[d.index(x)] = n

        # for Response-onset
        if x == 242:
            y = d[d.index(x) + 1]
            n = x * 1000 + y
            d[d.index(x)] = n

    pos_list = []

    for x in range(len(d)):
        if d[x] in (216, 217):
            pos_list.append(x)
    print(pos_list)

    for i in pos_list:
        if d[0] == 255:
            if i % 2 == 0:
                lst = d[pos_list[pos_list.index(i) - 1]: i]
                for x in lst:
                    if x in (241250, 242243):
                        n = d[i] * 1000000 + x
                        lst[lst.index(x)]=n
                        d[pos_list[pos_list.index(i)-1]:i]=lst
        else:
            if not i % 2 == 0:
                lst= d[pos_list[pos_list.index(i)-1]:i]
                for x in lst:
                    if x in (241250, 242243):
                        n = d[i] * 1000000 + x
                        lst[lst.index(x)] = n
                        d[pos_list[pos_list.index(i) - 1]:i] = lst

    events[:, 2] = d

    epochs = mne.Epochs(raw, events, event_id=event_id,
                        tmin=tmin, tmax=tmax,
                        proj=False, picks=picks,
                        baseline=baseline,
                        detrend=0, # DC offset
                        reject_by_annotation=True,
                        flat=flat,
                        reject=reject,
                        preload=True)


    bad_epoch_mask = abs_threshold(epochs, 75e-6)
    epochs.drop(bad_epoch_mask,reason = "absolute threshold")
    epochs.save(filename + '-epo.fif')
    return events

def abs_threshold(epochs, threshold):
    '''Compute boolean mask for dropping epochs based on absolute voltage threshold'''

    data = epochs.pick_types(eeg=True,misc=False,stim=False).get_data()
    # channels and times are last two dimension in MNE ndarrays,
    # and we collapse across them to get a (n_epochs,) shaped array
    rej = np.any( np.abs(data) > threshold, axis=(-1,-2))
    return rej

# Rolling definitions for functions for Russell's sanity and anyone learning -  def ss_preproc is a definition 


def retrieve(epochs, windows, subj=None, items=None):
    df = epochs.to_data_frame(picks=None, index=['epoch','time'],
                              scaling_time=1e3)
    eeg_chs = [c for c in df.columns if c not in ('condition')]
    # the order is important here! otherwise the shortcut
    # with items later won't work
    factors = ['epoch', 'condition']
    sel = factors + eeg_chs
    df = df.reset_index()

    retrieve = []
    for w in windows:
        temp = df[ df.time >= windows[w][0] ]
        dfw = temp[ temp.time <= windows[w][1] ]
        dfw_mean = dfw[sel].groupby(factors).mean()
        if subj:
            dfw_mean["subj"] = subj
        if items:
            dfw_mean["item"] = items
        dfw_mean["win"] = "{}..{}".format(*windows[w])
        dfw_mean["wname"] = w
        retrieve.append(dfw_mean)

    retrieve = pd.concat(retrieve)
    return retrieve


def ss_preproc(set_file):
    # read the raw file
    raw  = mne.io.read_raw_brainvision(set_file,
                                       eog=['SO1', 'IO1', 'LO1', 'LO2'],
                                       montage='standard_1020', preload=True)

    # making reference-- averaging (mean of A1+A2)
    raw = mne.io.add_reference_channels(raw, ['A1'])
    raw.set_eeg_reference(['A1','A2'])
    # filtering sensitivity either 0.1 or 0.3
    raw = raw.copy().filter(0.1, 30, l_trans_bandwidth='auto',
                            h_trans_bandwidth='auto', method='fir',
                            phase='zero', fir_window='hamming', n_jobs=2)
    # find eog events-- blink related
    eog_events = mne.preprocessing.find_eog_events(raw)
    n_blinks = len(eog_events)
    onset = eog_events[:, 0] / raw.info['sfreq'] - 0.25
    duration = np.repeat(0.5, n_blinks)
    raw.annotations = mne.Annotations(onset, duration,
                                      ['bad blink'] * n_blinks,
                                      orig_time=raw.info['meas_date'])
    raw.save(set_file + '-raw.fif', overwrite=True)


if __name__ == '__main__':
    data_path = '/home/dmalt/Data/Russel/Participant 11/'
    set_files = glob.glob(data_path + '*.vhdr')
    for filename in set_files:
        ss_preproc(filename)
        print (filename[0:-4])

    event_id={
            'S4b_Sti_correct':241250,'S4b_Sti_incorrect':241243,
            #'S4b_ERN_16':216242243,#'S4b_ERN_17':217242243,
            # 'S4b_correct_16':216241250,'S4b_correct_17':217241250,
            'S4b_Res_correct':242250, 'S4b_Res_incorrect':242243
            }

    # event_id = {
    #         'S1_Sti_correct':241250,'S1_Sti_incorrect':241243,
    #         'S1_Res_correct':242250, 'S1_Res_incorrect':242243
    #         }

    windows = {
        "Prestim100": (-100, 0),
        "ERN": (0, 150),
        "N2" : (160, 250),
        "P3" : (250, 400),
        }

    # glob.glob is the function to run through the files
    raw_files = [f for f in glob.glob(data_path + '*.fif')
                 if not f.endswith('-epo.fif')]

    for raw_fname in raw_files:
            evs = ss_epo(raw_fname, event_id=event_id, windows=windows)
            print (raw_fname[0:-4])



    ss_avg = dict() # new empty dict
    ss_wins = dict()
    # dict with empty list for each condition
    ss_avg_by_cond = {e:list() for e in event_id}
    ss_avg_all = list()
    retrieve_test = list()

    epo_files = glob.glob(data_path + '*-epo.fif')
    #set_files=['FLC_22.cnt-epo.fif']


    for epo_file in epo_files:
        epochs=mne.read_epochs(epo_file)
        evokeds = {str(cond):epochs[str(cond)].average() for cond in event_id}
        wins = retrieve(epochs, windows, epo_file)
        ss_avg[epo_file]= evokeds
        ss_wins[epo_file] = wins
        retrieve_test.append(ss_wins[epo_file])

        for cond in event_id:
            ss_avg_by_cond[cond].append(ss_avg[epo_file][cond])
            ss_avg_all.append(ss_avg[epo_file][cond])

        print(epo_file)

    # To retrieve the data into CSV file
    df=pd.concat(retrieve_test, ignore_index=False)
    df.to_csv('EveryoneS4ERP.csv', sep=',')

    # grand_all = mne.grand_average(ss_avg_all) #sanity check using pilot data
    # grand_avs = {cond:mne.grand_average(ss_avg_by_cond[cond]) for cond in event_id}

    # grand_all=mne.grand_average(ss_avg_all)
    # mne.viz.plot_compare_evokeds((grand_all), invert_y=True,picks=4)

    # plausibility = ("S1_Sti_correct","S1_Sti_incorrect")
    # stimuli = {k:grand_avs[k] for k in plausibility}

    # ERN = ("S4b_Res_correct","S4b_Res_incorrect")
    # ERN_all = {k:grand_avs[k] for k in ERN}

    # ERN = ("S1_Res_correct","S1_Res_incorrect")
    # ERN_all = {k:grand_avs[k] for k in ERN}

    # mne.viz.plot_compare_evokeds((stimuli),picks=2,invert_y=True)
    # mne.viz.plot_compare_evokeds((ERN_all),picks=2,invert_y=True)
    # mne.viz.plot_compare_evokeds((grand_all),invert_y=True)
