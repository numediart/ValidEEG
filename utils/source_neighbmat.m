function [neighboring_matrix] = source_neighbmat(sourcemodel_atlas,to_disp)

if nargin == 1
    to_disp = 1;
end

dipole_pos = sourcemodel_atlas.pos;
dist = sqrt(sum(diff(dipole_pos).^2,2));
thresh = median(dist); %the threshold is defined as the median distance between 2 following vertices

neighboring_matrix = zeros(length(sourcemodel_atlas.tissuelabel));
for i = 1:length(dipole_pos)
    dist = sqrt(sum((dipole_pos(i,:)-dipole_pos).^2,2));

    neighb_idx = find(dist < thresh);

    for idx = neighb_idx'
        id1 = sourcemodel_atlas.tissue(i);
        id2 = sourcemodel_atlas.tissue(idx);
        if  id1 ~= id2 && all([id1,id2])
    %         are_neighb = [are_neighb; [id1,id2]];
            neighboring_matrix(id1,id2) = 1;
            neighboring_matrix(id2,id1) = 1;
        end
    end
end

% neighborhood correction following regions properties
for i = 1:length(neighboring_matrix)
    % forbid left/right neighbours
    if mod(i,2)
        neighboring_matrix(i,2:2:end)=0;
    else
        neighboring_matrix(i,1:2:end)=0;
    end
    % ensure left/right symmetry
    for j = 1:i
        if neighboring_matrix(i,j)==1
            if mod(i,2)
                neighboring_matrix(i+1,j+1)=1;
                neighboring_matrix(j+1,i+1)=1;
            else
                neighboring_matrix(i-1,j-1)=1;
                neighboring_matrix(j-1,i-1)=1;
            end
        end
    end
end

if to_disp
    figure()
    imagesc(neighboring_matrix)
    title('Neighbouring matrix')
    xlabel('region index')
    ylabel('region index')
end
end