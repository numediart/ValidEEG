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
OSCIL.max_freq | TO DO | TO DO 
OSCIL.ampli | TO DO | TO DO 
OSCIL.modulation | TO DO | TO DO 
atlas | TO DO | TO DO 
dipoles_selection | TO DO | TO DO 
dipole_idx | TO DO | TO DO 
avoid_2nd_neighb | TO DO | TO DO 
headmodel | TO DO | TO DO 
elec | TO DO | TO DO 
channels | TO DO | TO DO 
benchmark_artf | TO DO | TO DO 
reconstr_source | TO DO | TO DO 
artifacts | TO DO | TO DO 
snr_source | TO DO | TO DO 
snr_eeg | TO DO | TO DO 
