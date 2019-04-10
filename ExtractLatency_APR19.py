import os.path as op
import numpy as np
from mne.io import Raw
import mne


def get_block(b_end_samples, sample):
    for i_block in b_end_samples:
        if b_end_samples[i_block][1] > sample > b_end_samples[i_block][0]:
            return i_block
        else:
            continue


if __name__ == '__main__':
    data_path = './data/Participant 11'
    raw_fname = op.join(data_path, 'MenTrainP11S4b.vhdr-raw.fif')
    raw = Raw(raw_fname)
    events = mne.find_events(raw)
    times = raw.times

    event_codes = {'start': 241, 'end': 242}
    samples_start = events[events[:, 2] == event_codes['start'], 0]
    times_start = times[samples_start]

    samples_end = events[events[:, 2] == event_codes['end'], 0]
    times_end = times[samples_end]


    latencies = times_end - times_start
    print('Computed latencies for %d epochs' % len(latencies))

    # -------- compute epoch block indices -------- #
    block_codes = [200 + i for i in range(4, 18)]
    print(block_codes)

    b_end_samples = {}
    for b_code in block_codes:
        b_end_samples[b_code] = events[events[:, 2] == b_code, 0]

    epoch_block_ids = []
    for sample in samples_start:
        epoch_block_ids.append(get_block(b_end_samples, sample))
    epoch_block_ids = np.array(epoch_block_ids)
    print('Computed epoch block indices for %d epochs' % len(epoch_block_ids))
    # --------------------------------------------- #

    # sanity check
    # check if number of epochs in each block is valid
    n_epochs = 120
    for i_block in b_end_samples:
        n_epochs_found = np.sum(epoch_block_ids == i_block)
        if  n_epochs_found == n_epochs:
            print('Checking number of epochs for block #',
                  str(i_block), ' --> ', str(n_epochs_found), '(OK)')
        else:
            print('ERROR: Block {} has {} epochs instead of {}'.format(
            i_block, n_epochs_found, n_epochs))
