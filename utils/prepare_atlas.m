function [atlas] = prepare_atlas(atlas_path,rois,roi_mat)
%% convert nifti atlas into FieldTrip sourcemodel structure

if nargin < 2
    error('Missing argument')
elseif nargin == 2
    roi_mat = {};
end

%% Realign the source model to the atlas
%load the raw source model to realign
try
    filename = 'utils/L.midthickness.8k_fs_LR.surf.gii';
    sourcemodel_raw = ft_read_headshape({filename, strrep(filename, 'L.', 'R.')});
catch
    error('No source model defined... Please download the default source model and add it to "utils" subfolder')
end



atlas = ft_read_atlas(atlas_path);
atlas = ft_convert_units(atlas,'mm');

atlas.tissue(~ismember(atlas.tissue,rois)) = 0; %Remove cerebellum and vermis (not in the sourcemodel)

if ~exist('transform_sourcemodel_atlas.mat','file')
    global PATH_TO_FIELDTRIP
    load(fullfile(PATH_TO_FIELDTRIP,'template/sourcemodel/standard_sourcemodel3d4mm.mat'))
    template_grid = sourcemodel;
    template_grid = ft_convert_units(template_grid,'mm');
    cfg            = [];
    cfg.atlas      = atlas;
    cfg.roi        = atlas.tissuelabel;  % here you can also specify a single label, i.e. single ROI
    cfg.inputcoord = 'mni';
    mask           = ft_volumelookup(cfg, template_grid);

    template_grid.inside = false(template_grid.dim);
    template_grid.inside(mask==1) = true;

    MRI_DATA = fullfile(PATH_TO_FIELDTRIP,'template/headmodel/standard_mri.mat');
    mri = ft_read_mri(MRI_DATA);
    cfg                = [];
    cfg.warpmni   = 'yes';
    cfg.template  = template_grid;
    cfg.nonlinear = 'yes'; % use non-linear normalization
    cfg.mri            = mri;
    sourcemodel        = ft_prepare_sourcemodel(cfg);
    sourcemodel = ft_convert_units(sourcemodel,'mm');
    template_sourcemodel = sourcemodel;
    template_sourcemodel.pos = sourcemodel.pos(sourcemodel.inside,:);

    cfg = [];
    cfg.individual.headshape= sourcemodel_raw;
    cfg.individual.headmodelstyle = 'surface';
    cfg.template.headshape= template_sourcemodel;
    cfg = ft_interactiverealign(cfg);
    transform_sourcemodel_atlas = cfg.m;
else
    load('utils/transform_sourcemodel_atlas.mat')
end
sourcemodel_realigned = ft_transform_geometry(transform_sourcemodel_atlas,sourcemodel_raw);

clear cfg datapath filename sourcemodel_raw

% save('transform_sourcemodel_atlas.mat', 'transform_sourcemodel_atlas');

%% Split the reconstructed sources in atlas
% load(fullfile(PATH_TO_SAVED_MAT,SUBJECT_NAME,'sourcemodel_realigned.mat'));
% load(fullfile(PATH_TO_SAVED_MAT,SUBJECT_NAME,'transform_sourcemodel_atlas.mat'));

figure()
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'tissue';
cfg.funcolormap  = 'jet';
ft_sourceplot(cfg, atlas)

% and call ft_sourceinterpolate:
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
atlas = ft_sourceinterpolate(cfg, atlas, sourcemodel_realigned);

T = transform_vox2ctf/transform_vox2acpc/transform_sourcemodel_atlas;
atlas = ft_transform_geometry(T,atlas);

if ~isempty(roi_mat)
    %example: roi_mat = {[3,5];[4,6];[7,9];[8,10];[11,13,15];[12,14,16];61;62};
    [atlas] = select_roi(atlas,roi_mat);
end

% save('atlas.mat', 'atlas');