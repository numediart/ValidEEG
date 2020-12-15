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
ft_defaults
cd(root)

%% create pseudo-EEG signal
fs = config.fsample;
n_dip = config.n_dipoles;
n_sess = config.n_sessions;
n_trial = config.n_trials;
n_artf = config.n_artifacts;
time = (0:config.pseudo_length*fs-1)/fs;

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
% figure;imagesc(neighboring_matrix);colormap(gray)

% prepare dipole simulation config
vol = load(config.headmodel);
vol = vol.(cell2mat(fieldnames(vol)));
elec = load(config.elec);
elec = elec.(cell2mat(fieldnames(elec)));
cfg      = [];
cfg.headmodel = vol;
cfg.elec = elec;
labels = upper(elec.label(str2num(config.channels)));
cfg.channel = labels;

% create event array: samples at which the template starts
duration = config.session_duration*60*fs;
event = randi(duration-length(time),n_sess,n_trial);

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
    all_artf{i} = all_artf{i}./max(abs(all_artf{i}),[],2);
    
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
    switch(config.type)
        case {'erp', 'ERP'}
            erp = load(config.ERP);
            erp = erp.(cell2mat(fieldnames(erp)));
            if ceil(erp.time(end)) ~= ceil(time(end))
                error('config.pseudo_length is different than the template ERP duration')
            end
            if erp.fsample ~= fs
                erp.erp = interp1(erp.time,erp.erp,time)';
                erp.erp(isnan(erp.erp)) = 0;
                erp.time = (0:fs-1)./fs;
                erp.fsample = fs;
            end
            
            signal = cell(n_trial,1);
            for trial = 1:n_trial
                erp_idx = randi(size(erp.erp,1),1,n_dip); %define each trial activity as one of the template ERPs
                tmp = erp.erp(erp_idx,:);
                % introduce causal interaction between the 2 first dipoles
                shift = ceil(fs*0.05); %shift by 50ms (for further causality study)
                tmp(1,end-shift+1:end) = 0;
                tmp(1,:) = circshift(tmp(1,:),shift);
                signal{trial} = tmp;

                % introduce frequency and phase shift on non-interacting sources
                if n_dip > 2
                    for dip = 3:n_dip
                        tmp = signal{trial}(dip,:);
                        freq_shift = rand(1)+0.5; %0.5 to 1.5 Hz frequency shift
                        ph_shift = rand(1)*2*pi;
                        signal{trial}(dip,:)= tmp.*cos(2*pi*freq_shift*time+ph_shift); %frequency modulation
                    end
                end
                tmp = awgn(signal{trial},config.snr_source,'measured'); %add white noise
                signal{trial} = tmp./max(tmp,[],2)*5; %normalize (max ampli=5µV)
                idx = event(s,trial):event(s,trial)+length(time)-1;
                pseudo_source{s}(:,idx) = pseudo_source{s}(:,idx) + signal{trial};
            end

        case {'oscil', 'OSCIL'}
            if isempty(config.freq)
                freq = (1:n_dip).*config.max_freq/n_dip;
                % impose 2 dipoles having the same frequency 
                % (for connectivity analysis purpose)
                rand_idx = randi(n_dip);
                new_idx = rand_idx+1;
                if new_idx > n_dip
                    new_idx = 1;
                end
                freq(new_idx) = freq(rand_idx);
            else
                freq = [];
                for f = fieldnames(config.freq)'
                    freq = [freq, config.freq.(f{1})];
                end
            end

            signal = cell(n_trial,1);
            for trial = 1:n_trial
                signal{trial} = zeros(n_dip,length(time));
                for i = 1:n_dip
                    if i > length(freq)
                        f = freq(mod(i-1,length(freq))+1);
                    else
                        f = freq(i);
                    end
                    rand_step = zeros(1,length(time));
                    step_idx = randi(length(time),1,floor(length(time)/10)); %a 10th of the samples will be frequency shifted
                    rand_step(step_idx) = randn(1,length(step_idx))*0.05; %random frequency shift(mean=0,std=0.05)
                    f_vec = filter(1,[1 -1], rand_step, f); %f(i) = f(i-1) + rand_step

                    % add abrupt change in freq to introduce causality afterwards
                    abrupt_times = [config.pseudo_length/3, 2*config.pseudo_length/3];
                    f_vec(and(time>abrupt_times(1), time<abrupt_times(2))) = ...
                        f_vec(and(time>abrupt_times(1), time<abrupt_times(2)))-(min(freq)/3);
                    f_vec(time>abrupt_times(2)) = f_vec(time>abrupt_times(2))+(min(freq)/3);

                    if i ~= new_idx
                        tmp = 5*sin(f_vec.*time.*2*pi);
                    else %interacting source
                        tmp = 5*sin(f_vec.*time.*2*pi - pi/3); %change phase of interacting source to make the difference with volume conduction
                        shift = ceil(fs*0.05); %shift by 50ms (for further causality study)
                        tmp(end-shift+1:end) = 0;
                        tmp = circshift(tmp,shift);
                    end
                    tmp = smooth(tmp,round(fs/50));
                    signal{trial}(i,:) = tmp./max(tmp)*5; %normalize (max ampli=5µV)
                end
                idx = event(s,trial):event(s,trial)+length(time)-1;
                pseudo_source{s}(:,idx) = pseudo_source{s}(:,idx) + signal{trial};
            end
    end
    
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
            else
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

    cfg.dip.pos = dipole_pos;
    cfg.dip.mom = dipole_mom;
    cfg.dip.signal = pseudo_source(s);
    cfg.fsample = fs;
    pseudo_eeg{s} = ft_dipolesimulation(cfg);
    tmp = pseudo_eeg{s}.trial{1};
    pseudo_eeg{s}.trial{1} = tmp/max(tmp(:)); %Normalization
    
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
eegplot(pseudo_eeg{s}.trial{1}, 'srate', fs,'position',[0 30 1535 780])
eegplot(tmp, 'srate', fs,'position',[0 30 1535 780])
pseudo_eeg{s}.time{1}(pseudo_artf{s}.sample)
%% save
save('pseudo_data/event.mat','event')
save('pseudo_data/dipole.mat','pseudo_dipole')
save('pseudo_data/source.mat','pseudo_source')
save('pseudo_data/artifact.mat','pseudo_artf')




