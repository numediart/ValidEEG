function [roi_atlas] = select_roi(sourcemodel_atlas,roi_mat)

% This function transform the input sourcemodel_atlas into a new one 
% considering the targeted ROIs defined in roi_mat.

% roi_mat is a Nx1 cell-array where each cell contain the ID of the old
% label(s) that has/have to be relabeled to the corresponding row's ID.
% Example:
%     old labels IDs: 1,2,3,4,5
%     new labels: 1,[2+3],5 -> the 2nd and 3rd labels are marged and the 4th is
%     removed.
%     => roi_mat = {1;[2,3];5};

roi_atlas = sourcemodel_atlas;
roi_atlas.tissue = zeros(length(sourcemodel_atlas.tissue),1);
roi_atlas.tissuelabel = cell(length(roi_mat),1);
for i = 1:length(roi_mat)
    roi_atlas.tissuelabel{i} = '';
    for j = roi_mat{i}
        roi_atlas.tissue(sourcemodel_atlas.tissue == j) = i;
        if ~isempty(roi_atlas.tissuelabel{i})
            roi_atlas.tissuelabel{i} = [roi_atlas.tissuelabel{i} ' and '];
        end
        roi_atlas.tissuelabel{i} = [roi_atlas.tissuelabel{i}, sourcemodel_atlas.tissuelabel{j}];
    end
end
        
end