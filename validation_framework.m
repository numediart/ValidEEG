%% Fashionable validation framework for ERP and oscillatory source reconstruction using Fieldtrip

%% set paths and load config
mfilename = 'validation_framework.m';
root = fileparts(which(mfilename));
addpath(genpath(root))
cd(root)

%read config file
fname = 'config.json';
fid = fopen(fname);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
config = jsondecode(str);

addpath(config.PATH_TO_FIELDTRIP)
addpath(genpath(config.PATH_TO_SEREEGA))
ft_defaults
cd(root)

%% create pseudo-EEG signal
fs = config.fsample;
n_dip = config.n_dipoles;
n_sess = config.n_sessions;
n_trial = config.n_trials;
n_artf = config.n_artifacts;
time = (0:config.pseudo_length*fs-1)/fs;

% define epochs (required by SEREEGA)
epochs        = struct();
epochs.n      = n_trial; % the number of epochs to simulate
epochs.srate  = fs; % their sampling rate in Hz
epochs.length = config.pseudo_length * 1000; % their length in ms

% load atlas and compute neighboring matrix
atlas = load(config.atlas);
atlas = atlas.(cell2mat(fieldnames(atlas)));
neighboring_matrix = source_neighbmat(atlas,0);
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

% prepare dipole simulation config
vol = load(config.headmodel);
vol = vol.(cell2mat(fieldnames(vol)));
elec = load(config.elec);
elec = elec.(cell2mat(fieldnames(elec)));
cfg      = [];
cfg.headmodel = vol;
cfg.elec = elec;
cfg.channel = elec.label(str2num(config.channels));

labels = upper(cfg.channel);

% create event array: samples at which the template starts
duration = config.session_duration*60*fs;
if isempty(config.event)
    min_space = .25*fs; %events are spaced by 0.25s minimum 
    count = 0;
    event = -1;
    while(any(event(:)<0)) %prevent negative samples
        event = randi(duration-length(time),n_sess,n_trial);
        event = sort(event,2); %sort the events by appearance time
        event_space = diff(event,[],2);
        while any(event_space(:) < (length(time)+min_space)) %ensure minimum space
            event(event_space < (length(time)+min_space)) = ...
                event(event_space < (length(time)+min_space)) - (length(time)+min_space);
            event_space = diff(event,[],2);
        end
        count = count + 1;
        if count > 1000
            error('Too many trials! Tips: increase duration or decrease n_trials')
        end
    end
    save('pseudo_data/event.mat','event')
else
    load(config.event)
end

% define the artifacts from the template
artf = load(config.artifacts);
artf = artf.(cell2mat(fieldnames(artf)));

if ~all(ismember(lower(labels),lower(artf.label))) %deal with channel difference between artifact template and pseudo-eeg
    % from MCN to 1010-1020 systems
    if strcmp(elec.type,'eeg1010') || strcmp(elec.type,'eeg1020')
        for i = 3:4
            artf.label(strcmp(artf.label,['T' num2str(i)])) = {['T' num2str(i+4)]};
            artf.label(strcmp(artf.label,['T' num2str(i+2)])) = {['P' num2str(i+4)]};
        end
    end
    artf.label = upper(artf.label(ismember(lower(artf.label),lower(labels))));
    
    % compute neighboring matrix of pseudo-eeg sensors to further
    % interpolate the artifact template to every channels
    cfg_neighb = [];
    cfg_neighb.elec = elec;
    cfg_neighb.channel = labels;
    cfg_neighb.template_chan = artf.label;
    neighboring_artf = eeg_neighbmat(cfg_neighb);
end    

if artf.fsample ~= fs %ensure same sampling rate
    time_artf = 0:1/fs:artf.time(end);
    for field = fieldnames(artf.artf)'
        field = field{1};
        for c = 1:length(artf.artf.(field))
            current_artf = artf.artf.(field){c};
            current_artf = interp1(artf.time(1:size(current_artf,2)),current_artf',time_artf(time_artf<=artf.time(size(current_artf,2))))';
            current_artf(isnan(current_artf)) = 0;
            artf.artf.(field)(c) = {current_artf};
        end
    end
    artf.time = time_artf;
    artf.fsample = fs;
end

all_artf = struct2cell(artf.artf);
artf_type.name = fieldnames(artf.artf);
artf_type.idx = [];
for i = 1:length(all_artf)
    artf_type.idx = [artf_type.idx;i*ones(length(all_artf{i}),1)];
end
all_artf = [all_artf{:}];
for i = 1:length(all_artf)
    all_artf{i} = all_artf{i}./prctile(abs(all_artf{i}(:)),90);
    
    %interpolate missing channels
    missing_chan = ~ismember(lower(labels),lower(artf.label));
    if any(missing_chan) %if some channels are missing
        tmp = zeros(length(labels),size(all_artf{i},2));
        for chan = 1:length(labels)
            if missing_chan(chan)
                neighb_lab = labels(neighboring_artf(chan,:));
                neighb = neighb_lab(ismember(neighb_lab,artf.label));

                tmp(chan,:) = mean(all_artf{i}(ismember(artf.label,neighb),:));
            else
                tmp(chan,:) = all_artf{i}(strcmpi(artf.label,labels(chan)),:);
            end
        end
    end
    all_artf{i} = tmp;
end

pseudo_eeg = cell(n_sess,1);
pseudo_dipole = cell(n_sess,1);
pseudo_source = cell(n_sess,1);
pseudo_source(:) = {zeros(n_dip,duration)};
pseudo_artf = cell(n_sess,1);
for s = 1:n_sess
    % random dipole position
    rng('shuffle');
    neighb_cond = 0;
    non_zero_dip = find(atlas.tissue ~= 0);
    count = 0;
    clc;
    while(~neighb_cond)
        dipole_idx = non_zero_dip(randi(length(non_zero_dip),1,n_dip));
        % avoid dipoles to be in neighbor regions (or neighbor of neighbor)
        for dip = 1:n_dip
            neighb_cond = 1;
            [~,neighbors] = find(neighboring_matrix(atlas.tissue(dipole_idx(dip)),:));
            if any(ismember(atlas.tissue(dipole_idx),neighbors))
                neighb_cond = 0;
                break
            elseif config.avoid_2nd_neighb
                for neighb = neighbors
                    [~,second_neighb] = find(neighboring_matrix(neighb,:));
                    if sum(ismember(atlas.tissue(dipole_idx),second_neighb))>1
                        neighb_cond = 0;
                        break
                    end
                end
            end
            if ~neighb_cond
                break
            end
        end
        dipole_pos = atlas.pos(dipole_idx, :);
    end
    pseudo_dipole{s}.pos = dipole_pos;
    pseudo_dipole{s}.region = atlas.tissue(dipole_idx);

    % random dipole moment
    dipole_mom = zeros(3*n_dip,1);
    for i = 1:n_dip
        r = rand(3,1);
        r = r/norm(r);
        dipole_mom(3*(i-1)+1:3*i) = r';
    end
    
    switch(config.type)
        case {'erp', 'ERP'}
            peaks = split(config.ERP.peaks,',');
            latency = regexp(peaks,'\d*','match');
            latency = cellfun(@(x) str2double(x{1}),latency)';
            amplitude = str2num(config.ERP.ampli) .* cell2mat(cellfun(@(x) ...
                contains(x,'p','IgnoreCase',true)-contains(x,'n','IgnoreCase',true),...
                peaks,'un',0))';
            width = str2num(config.ERP.width);

            erp = struct();
            erp.peakLatency = latency;      % in ms, starting at the start of the epoch
            erp.peakWidth = width;        % in ms
            erp.peakAmplitude = amplitude;      % in microvolt
            erp.peakLatencyDv = repmat(50,1,length(peaks));
            erp.peakAmplitudeDv = .2*amplitude;
            erp.peakWidthDv = .5*width;
            erp.peakAmplitudeSlope = -.75*amplitude; %introduce stimuli habituation

            erp = utl_check_class(erp, 'type', 'erp');

            tmp = pseudo_source{s};
            for trial = 1:n_trial
                idx = event(s,trial):event(s,trial)+length(time)-1;
                for dip = 1:n_dip
                    signal = erp_generate_signal_fromclass(erp, epochs, 'epochNumber', trial);
                    if atlas.pos(dipole_idx(dip), 1) > 0 %frontal region -> polarity inversion
                        tmp(dip,idx) = tmp(dip,idx) - signal;
                    else
                        tmp(dip,idx) = tmp(dip,idx) + signal;
                    end
                end
            end
            tmp = awgn(tmp,config.snr_source,'measured'); %add white noise
            pseudo_source{s} = tmp;

        case {'oscil', 'OSCIL'}
            tmp = pseudo_source{s};
            if isempty(config.OSCIL.freq)
                freq = (1:n_dip).*config.OSCIL.max_freq/n_dip;
                % impose 2 dipoles having the same frequency 
                % (for connectivity analysis purpose)
                rand_idx = randi(n_dip);
                new_idx = rand_idx+1;
                if new_idx > n_dip
                    new_idx = 1;
                end
                freq(new_idx) = freq(rand_idx);
                freq = [freq'-2,freq'+2];
            else
                freq = [];
                for f = fieldnames(config.OSCIL.freq)'
                    freq = [freq; config.OSCIL.freq.(f{1})'];
                end
            end
            amplitude = str2num(config.OSCIL.ampli);
            
            for dip = 1:n_dip
                if dip > size(freq,1)
                    f = freq(mod(dip-1,size(freq,1))+1,:);
                else
                    f = freq(dip,:);
                end
                
                ersp = struct();
                ersp.frequency = [f(1)-.2*diff(f), f, f(2)+.2*diff(f)]; % frequency band
                ersp.amplitude = amplitude(dip); % in microvolt
                ersp.modulation = config.OSCIL.modulation;
                ersp.modFrequency = 3;
                ersp.modPhase = .25;

                ersp = utl_check_class(ersp, 'type', 'ersp');

                for trial = 1:n_trial
                    ersp.phase = rand(1)*2*pi;
                    idx = event(s,trial):event(s,trial)+length(time)-1;
                    signal = ersp_generate_signal_fromclass(ersp, epochs, 'epochNumber', trial);
                    tmp(dip,idx) = tmp(dip,idx) + signal;
                end
            end
            tmp = awgn(tmp,config.snr_source,'measured'); %add white noise
            pseudo_source{s} = tmp;            
    end

    cfg.dip.pos = dipole_pos;
    cfg.dip.mom = dipole_mom;
    cfg.dip.signal = pseudo_source(s);
    cfg.fsample = fs;
    pseudo_eeg{s} = ft_dipolesimulation(cfg);
    tmp = pseudo_eeg{s}.trial{1};
%     pseudo_eeg{s}.trial{1} = tmp/max(tmp(:)); %Normalization
    pseudo_eeg{s}.trial{1} = tmp./max(abs(tmp(:))); %Normalization
    
    % add artifacts randomly selected from the templates
    pseudo_artf{s} = cell(n_artf,3);
    rand_artf = randi(length(all_artf),n_artf,1);
    rand_time = sort(randi(duration-length(artf.time),n_artf,1));

    pseudo_artf{s} = {rand_time, rand_artf, artf_type.name(artf_type.idx(rand_artf))};
    pseudo_artf{s} = cell2struct(pseudo_artf{s},{'sample' 'index' 'type'},2);
    pseudo_artf{s} = struct2table(pseudo_artf{s});
    
    tmp = pseudo_eeg{s}.trial{1};
    for i = 1:n_artf
        tmp(:,rand_time(i):rand_time(i)+size(all_artf{rand_artf(i)},2)-1) = ...
            tmp(:,rand_time(i):rand_time(i)+size(all_artf{rand_artf(i)},2)-1) + all_artf{rand_artf(i)};
    end
    
    % add pink noise
    noise = pinknoise(size(tmp,1),size(tmp,2));
    pRMS = rms(noise,2).^2;
    tmp = tmp + noise./(pRMS.*10^(config.snr_eeg/20));

    pseudo_eeg{s}.trial{1} = tmp;
    
    % save pseudo-eeg of each session (to avoid out of memory)
    eeg = pseudo_eeg{s};
    save(['pseudo_eeg/session_' num2str(s) '.mat'],'eeg')
%     save(['pseudo_eeg/session' num2str(s)],'pseudo_eeg')
end
% eegplot(pseudo_source{s}, 'srate', fs,'position',[0 30 1535 780])
% eegplot(pseudo_eeg{s}.trial{1}, 'srate', fs,'position',[0 30 1535 780])
% eegplot(tmp, 'srate', fs,'position',[0 30 1535 780])
% pseudo_eeg{s}.time{1}(pseudo_artf{s}.sample)
%% save
save('pseudo_data/dipole.mat','pseudo_dipole')
save('pseudo_data/source.mat','pseudo_source')
save('pseudo_data/artifact.mat','pseudo_artf')

eegplot(pseudo_eeg{s}.trial{1}, 'srate', fs,'position',[0 30 1535 780])

% eegplot(tmp, 'srate', fs,'position',[0 30 1535 780])
figure('position',[0 30 1000 300]);plot(0:1/fs:5-1/fs, tmp(:,100*fs:(105-1/fs)*fs))
xlabel('time (s)')
ylabel('Amplitude (µV)')
title(sprintf('Pseudo-source signals\nA'))
legend(atlas.tissuelabel(atlas.tissue(dipole_idx)))
ax = gca;
ax.TitleFontSizeMultiplier = 1.5;




% figure('position',[0 900 1100 1100]);imagesc(neighboring_matrix);colormap(gray)
% truesize([800,800])
% xticks(1:length(atlas.tissuelabel))
% xticklabels(atlas.tissuelabel)
% xtickangle(90)
% yticks(1:length(atlas.tissuelabel))
% yticklabels(atlas.tissuelabel)
% title('source neighboring matrix')
% ax = gca;
% ax.TitleFontSizeMultiplier = 2;

% fs = 2048
% s = 5
% eegplot(eeg.trial{1}, 'srate', fs,'position',[0 30 1535 780])
