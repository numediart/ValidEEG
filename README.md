# Validation Framework Source Reconstruction
 This repository contains all the codes and data relative to "La Fisca et al.,  A Versatile Validation Framework for ERP and Oscillatory Source Reconstruction Using FieldTrip, 2021"

Configuration description and template availability: TO BE DONE.

Description of configuration parameters:
Parameter | Pattern | Definition
--- | --- | --- 
benchmark | true/false | run benchmark (template parameters) 
evaluation | true/false | run the evaluation (set to true when reconstructed sources are computed) 
PATH_TO_FIELDTRIP | path name | Path to FieldTrip toolbox 
PATH_TO_SEREEGA | path name | Path to SEREEGA toolbox 
n_sessions | int | number of sessions over which to generate pseudo-data (different dipoles for each session)
n_dipoles | int | number of dipoles defined as source (recommended between 2 and 5)
n_trials | int | number of occurences of source activation within one session
n_artifacts | int | number of artifactual segments occuring within one session
fsample | int | sampling frequency of generated EEG
session_duration | float | duration in minutes of each session
pseudo_length | float | duration in seconds of each trial (source activation)
event | n_sessions x n_trials matrix | matrix defining the starting time of each trial
type | ERP/OSCIL | definition of source signal type (event-related or oscillatory)
ERP.peaks | Pxxx/Nxxx series | definition of ERP peaks as series of positive/negative peaks (e.g. P100,N200,P300)
ERP.ampli | float series | maximum amplitude of each peak
ERP.width | int | width of each peak corresponding to 6 standard deviation
OSCIL.freq | "f1": [min_f, max_f], "f2": [min_f2,max_f2],... | definition of the frequency band of each desired oscillation
OSCIL.max_freq | int | maximum default frequency if OSCIL.freq is not defined 
OSCIL.ampli | float series | maximum amplitude of each oscillation
OSCIL.modulation | none/ampmod | enable or not amplitude modulation of the predefined osciillations 
atlas | mesh structure (FieldTrip) | FieldTrip atlas structure defining targeted brain regions (cf. ft_read_atlas)
dipoles_selection | n_sessions x n_dipoles matrix | matrix defining the index of the selected dipoles (relatively to the atlas structure) 
avoid_2nd_neighb | true/false | avoid or not the selected dipoles to be second neighbors (i.e., neighbor of neighbor) of each other
headmodel | headmodel structure (FieldTrip) | head model FieldTrip structure corresponding to the volume conduction model (cf. ft_prepare_headmodel)
elec | elec structure (FieldTrip) | electrode FieldTrip structure describing the EEG sensors (cf. ft_datatype_sens)
channels | 1 x n_channels vector | vector defining which channels to work with (relatively to the elec structure)
reconstr_source | n_sessions x 1 structure | signal and label of the reconstructed sources from the generated pseudo-EEG (cf. template/reconstr_source.mat)
artifacts | struct(artf,time,fsample,label) | structure containing the artifactual segments (artf field) with specific longest time vector (1 x n_sample), sampling frequency (fsample field) and channel labels (n_channel x 1 cell) (cf.template/artf_template.mat)
snr_source | int | signal-to-noise ratio of generated source signal 
snr_eeg | int | signal-to-noise ratio of generated EEG signal 
