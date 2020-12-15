function neighboring_matrix = eeg_neighbmat(cfg)

elec = cfg.elec;
channel = cfg.channel;
template_chan = cfg.template_chan;
neighbordist = 10; %start with neighbor distance of 1cm

nsensors = length(channel);
chan_idx = find(ismember(lower(elec.label),lower(channel)));

neighboring_matrix = zeros(nsensors,nsensors);

% increase neighbordist until we have enough neighbors to interpolate the
% template channels (to the elec channels)
while(any(sum(neighboring_matrix,2)<(ceil(length(channel)/length(template_chan)))))
    % compute the distance between all sensors
    dist = zeros(nsensors,nsensors);
    idx = 1;
    for i = chan_idx'
        dist(idx,:) = sqrt(sum((elec.chanpos(chan_idx,:) - repmat(elec.chanpos(i,:), nsensors, 1)).^2,2))';
        idx = idx+1;
    end

    % find the neighboring electrodes based on distance
    % later we have to restrict the neighboring electrodes to those actually selected in the dataset
    neighboring_matrix = (dist<neighbordist);

    % electrode istelf is not a neighbour
    neighboring_matrix = (neighboring_matrix .* ~eye(nsensors));
    
    neighbordist = neighbordist + 10; %add 10cm for each iteration
end
neighboring_matrix = logical(neighboring_matrix);

% figure;imagesc(neighboring_matrix)