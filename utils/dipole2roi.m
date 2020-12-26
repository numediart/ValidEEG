function [source_roi] = dipole2roi(source_dipole,roi_atlas)

% This function transforms the dipole-by-dipole source activity into
% region-by-region source activity using SVD

%% Convert 3-dimensional moments into 1-dimensional moments using SVD
ndipole = length(source_dipole.pos);
ntime = length(source_dipole.time);
tmp = NaN(ndipole,ntime);
for j = 1:ndipole
    moments = source_dipole.avg.mom{j};
    if ~isempty(moments)
        [u, ~, ~] = svd(moments, 'econ'); %Find the dimension explaining/representing the most variance;
        m = u(:,1)'*moments;
        tmp(j,:) = m;
    end
end
source_dipole.avg.pow = source_dipole.avg.pow;
source_dipole.mom = tmp;


%% Parcellate the source signal using eigenvectors
roi_atlas.pos = source_dipole.pos; % otherwise the parcellation won't work

cfg = [];
cfg.method       = 'eig';
cfg.parcellation = 'tissue';
cfg.parameter    = [{'pow'};{'mom'}];
source_roi = ft_sourceparcellate(cfg, source_dipole, roi_atlas);

% properly reshape the parameters: the 1st dimension is the one on which we
% will concatenate the data (trial)
param = cfg.parameter;
for i = 1:length(param)
    source_roi.(param{i}) = permute(source_roi.(param{i}), [2 1 3]);
    source_roi.([param{i} 'dimord']) = 'time_chan';
end

% add the position of the ROI centers
source_roi.pos = NaN(length(source_roi.label),3);
for i = 1:length(source_roi.label)
    source_roi.pos(i,:) = mean(roi_atlas.pos(roi_atlas.tissue==i,:));
end

end