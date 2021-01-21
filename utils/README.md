# Utils functions description

## Neighboring matrices
eeg_neighbmat and source_neighbmat functions compute the neighboring matrix in sensor and source domain respectively.

eeg_neighbmat takes as input a cfg parameter composed of elec, channel and template_chan fields defined as follow:
* elec: sensors description (cf. ft_read_sens)
* channel: n_chan_elec x 1 cell with the labels of the desired channels in elec
* template_chan: n_chan_artf x 1 cell with the labels of the template artifact signals (in benchmark data, artifact signals are given on a 19-channels set-up (Hamid et al., 2021))

Here is an example:
```
cfg = [];
cfg.elec = elec;
cfg.channel = labels;
cfg.template_chan = artf.label;
neighboring_artf = eeg_neighbmat(cfg);
```

source_neighbmat takes as input the atlas and a boolean defining whether or not the neighboring matrix should be displayed at the end.

Here is an example:
`neighboring_matrix = source_neighbmat(atlas,0);`

## Custom atlas
There are 2 possible ways to customize your atlas: modify an existing one through `select_roi` function or create a new one from an atlas in NIfTI format through `prepare_atlas` function as showed in the following examples
1. ```
	roi_mat = {[3,5];[4,6];[7,9];[8,10];[11,13,15];[12,14,16];61;62};
    [atlas] = select_roi(atlas,roi_mat);
    ```
2. ```
	atals_path ='FieldTrip\template\atlas\aal\ROI_MNI_V4.nii';
	rois = [1:90];
	atlas = prepare_atlas(atlas_path,rois)
	```
Note: you can combine both functions by passing `roi_mat` as third argument of the `prepare_atlas` function.

To realign the atlas to a specific head model, rename or delete `transform_sourcemodel_atlas.mat` from utils.

## Compute ROI-by-ROI source activity
To perform the evaluation through the proposed framework, you need to get the source activity as 1 signal per brain region. This is what the `dipole2roi` function does. Basically, this function takes the output of `ft_sourceanalysis` and the atlas as inputs to compute ROI-by-ROI source signal. Example:
```
...
source_dipole = ft_sourceanalysis(cfg, timelock);

source_roi = dipole2roi(source_dipole,atlas);
reconstr_source.signal = source_roi.pow';
reconstr_source.label = source_roi.label';
```
`reconstr_source` is the required structure to run the evaluation process.