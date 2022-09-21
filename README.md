# Microstates Analysis of a 64-channels EEG signal recorded from an healthy subject during closed-eyes resting condition.

Electroencephalography (EEG) is a very useful method to study the electrophysiology of the brain in a resting state with high temporal resolution. In particular, it allows characterizing spontaneous brain activity that is not detected when a subject performs cognitive tasks.

Several analytical approaches have been proposed to extract information from the EEG signal, one of which, called Microstate Analysis, considers the multichannel EEG as a series of semi-stable microstates each of which is characterized by a unique topography of electrical potentials with respect to channels. 

Aim of this project is to identify the microstates of a 64-channel EEG signal of a subject at rest and with eyes closed and then study its variation over time. 

## Contents

### Analysis
- `Progetto_4C.mlx` is the Matlab Live Script with data processing and analysis of the parameters of interest, using clustering algorithm such as *k-means* and *Gaussian Mixtures*. This file contains detailed explanations of every steps, the results and discussions (in Italian).
- `script.m` contains the plain code to obtain the results. 

### Data
- `S092R02.edf` Data are publicly available on [PhysioNet](https://physionet.org) in .edf format. The 64 channels correspond to 64 electrodes placed according to the international standard 10-10 system, whose signals were recorded using the BCI2000 system. 

### Utility 
- `edfread.m` script used to import .edf files.
- `gfp_plot.mlx`
- `init_ms_seq.m`
- `plot_topography.m`
