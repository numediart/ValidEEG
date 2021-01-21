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

evaluation = false;
if config.evaluation
    evaluation = true;
    % Load the reconstructed sources
    reconstr_source = load(config.reconstr_source);
    reconstr_source = reconstr_source.(cell2mat(fieldnames(reconstr_source)));
end

if config.benchmark
    fname = 'template/config_template.json';
    fid = fopen(fname);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    config = jsondecode(str);
end

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

% prepare dipole simulation config
vol = ft_read_headmodel(config.headmodel);
elec = ft_read_sens(config.elec);
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

switch(evaluation)
    case 0 % create pseudo-EEG signal
        % define the artifacts from the template dataset
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
            all_artf{i} = all_artf{i}./prctile(abs(all_artf{i}(:)),99);

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

        dipole_idx = zeros(n_sess,n_dip);
        pseudo_eeg = cell(n_sess,1);
        pseudo_dipole = cell(n_sess,1);
        pseudo_source = cell(n_sess,1);
        pseudo_source(:) = {zeros(n_dip,duration)};
        pseudo_artf = cell(n_sess,1);

        if config.benchmark
            pseudo_artf = load(config.benchmark_artf);
            pseudo_artf = pseudo_artf.(cell2mat(fieldnames(pseudo_artf)));
            pseudo_dipole = load(config.benchmark_dipole);
            pseudo_dipole = pseudo_dipole.(cell2mat(fieldnames(pseudo_dipole)));
        end

        if ~isempty(config.dipoles_selection)
            dipole_idx = load(config.dipoles_selection);
            dipole_idx = dipole_idx.(cell2mat(fieldnames(dipole_idx)));
        end

        for s = 1:n_sess
            if ~config.benchmark && isempty(config.dipoles_selection)
                % random dipole position
                rng('shuffle');
                neighb_cond = 0;
                non_zero_dip = find(atlas.tissue ~= 0);
                count = 0;
                clc;
                while(~neighb_cond)
                    dipole_idx(s,:) = non_zero_dip(randi(length(non_zero_dip),1,n_dip));
                    % avoid dipoles to be in neighbor regions (or neighbor of neighbor)
                    for dip = 1:n_dip
                        neighb_cond = 1;
                        [~,neighbors] = find(neighboring_matrix(atlas.tissue(dipole_idx(s,dip)),:));
                        if any(ismember(atlas.tissue(dipole_idx(s,:)),neighbors))
                            neighb_cond = 0;
                            break
                        elseif config.avoid_2nd_neighb
                            for neighb = neighbors
                                [~,second_neighb] = find(neighboring_matrix(neighb,:));
                                if sum(ismember(atlas.tissue(dipole_idx(s,:)),second_neighb))>1
                                    neighb_cond = 0;
                                    break
                                end
                            end
                        end
                        if ~neighb_cond
                            break
                        end
                    end
                    dipole_pos = atlas.pos(dipole_idx(s,:), :);
                end
                % random dipole moment
                dipole_mom = zeros(3*n_dip,1);
                for i = 1:n_dip
                    r = rand(3,1);
                    r = r/norm(r);
                    dipole_mom(3*(i-1)+1:3*i) = r';
                end
                pseudo_dipole{s}.pos = dipole_pos;
                pseudo_dipole{s}.region = atlas.tissue(dipole_idx(s,:));
                pseudo_dipole{s}.mom = dipole_mom;

            elseif ~config.benchmark
                dipole_pos = atlas.pos(dipole_idx(s,:),:);
                % random dipole moment
                dipole_mom = zeros(3*n_dip,1);
                for i = 1:n_dip
                    r = rand(3,1);
                    r = r/norm(r);
                    dipole_mom(3*(i-1)+1:3*i) = r';
                end
                pseudo_dipole{s}.pos = dipole_pos;
                pseudo_dipole{s}.region = atlas.tissue(dipole_idx(s,:));
                pseudo_dipole{s}.mom = dipole_mom;

            else
                dipole_pos = pseudo_dipole{s}.pos;
                dipole_mom = pseudo_dipole{s}.mom;
                if dipole_pos ~= atlas.pos(dipole_idx(s,:), :)
                    error("the dipole indices and positions does not match... Check template parameters")
                end
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
                            if atlas.pos(dipole_idx(s,dip), 1) > 0 %frontal region -> polarity inversion
                                tmp(dip,idx) = tmp(dip,idx) - signal;
                            else
                                tmp(dip,idx) = tmp(dip,idx) + signal;
                            end
                        end
                    end
                    % add pink noise
                    noise = pinknoise(size(tmp,1),size(tmp,2));
                    pRMS = rms(noise,2).^2;
                    tmp = tmp + noise./(pRMS.*config.snr_source);
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
                    % add pink noise
                    noise = pinknoise(size(tmp,1),size(tmp,2));
                    pRMS = rms(noise,2).^2;
                    tmp = tmp + noise./(pRMS.*config.snr_source);
                    pseudo_source{s} = tmp;            
            end

            cfg.dip.pos = dipole_pos;
            cfg.dip.mom = dipole_mom;
            cfg.dip.signal = pseudo_source(s);
            cfg.fsample = fs;
            cfg.relnoise = 1/config.snr_eeg;
            pseudo_eeg{s} = ft_dipolesimulation(cfg);
            tmp = pseudo_eeg{s}.trial{1};
            pseudo_eeg{s}.trial{1} = tmp./max(abs(tmp(:))); %Normalization

%             %(uncomment) test the reconstruction without artifact
%             eeg = pseudo_eeg{s};
%             test_path = '..\test';
%             save([test_path '\session_' num2str(s) '.mat'],'eeg')
%             continue

            if ~config.benchmark
                % add artifacts randomly selected from the templates
                pseudo_artf{s} = cell(n_artf,3);
                rand_artf = randi(length(all_artf),n_artf,1);
                rand_time = sort(randi(duration-length(artf.time),n_artf,1));

                pseudo_artf{s} = {rand_time, rand_artf, artf_type.name(artf_type.idx(rand_artf))};
                pseudo_artf{s} = cell2struct(pseudo_artf{s},{'sample' 'index' 'type'},2);
                pseudo_artf{s} = struct2table(pseudo_artf{s});
            else
                rand_artf = pseudo_artf{s}.rand_artf;
            end

            tmp = pseudo_eeg{s}.trial{1};
            for i = 1:n_artf
                tmp(:,rand_time(i):rand_time(i)+size(all_artf{rand_artf(i)},2)-1) = ...
                    tmp(:,rand_time(i):rand_time(i)+size(all_artf{rand_artf(i)},2)-1) + all_artf{rand_artf(i)}./4;
            end

            pseudo_eeg{s}.trial{1} = tmp;

            % save pseudo-eeg of each session (to avoid out of memory)
            eeg = pseudo_eeg{s};
            save(['pseudo_eeg/session_' num2str(s) '.mat'],'eeg')

            % save pseudo_data
            save('pseudo_data/dipole_idx.mat','dipole_idx')
            save('pseudo_data/dipole.mat','pseudo_dipole')
            save('pseudo_data/source.mat','pseudo_source')
            save('pseudo_data/artifact.mat','pseudo_artf')

        end
        % eegplot(pseudo_source{s}, 'srate', fs,'position',[0 30 1535 780])
        % eegplot(pseudo_eeg{s}.trial{1}, 'srate', fs,'position',[0 30 1535 780])
        % eegplot(tmp, 'srate', fs,'position',[0 30 1535 780])
        % pseudo_eeg{s}.time{1}(pseudo_artf{s}.sample)

    case 1
        %% Evaluation
        % The reconstructed sources must be given as a n_sessions structure which 
        % fields are "signal" (as n_regions*trial_length matrix) and "label" 
        % (n_regions*1 cell)
        
        % Compare with actual pseudo-source regions
        bins = [1,.5,0,-1];
        score = zeros(n_sess,n_dip);
        for s = 1:n_sess
            region_pow = rms(reconstr_source(s).signal,2);
            high_region = [];
            while(length(high_region) < n_dip)
                [~,max_idx] = max(region_pow);
                region_pow(max_idx) = 0;

                roi_idx = find(ismember(atlas.tissuelabel,reconstr_source(s).label(max_idx)));

                [~,neighbors] = find(neighboring_matrix(roi_idx,:));
                if ~any(ismember(high_region,neighbors))
                    high_region = [high_region, roi_idx];
                end
            end

            true_regions = pseudo_dipole{s}.region;
            for dip = 1:n_dip
                neighbors = find(neighboring_matrix(true_regions(dip),:));
                neighb_neighbors = [];
                for neighb = neighbors
                    [~,second_neighb] = find(neighboring_matrix(neighb,:));
                    neighb_neighbors = [neighb_neighbors, second_neighb];
                end
                if any(ismember(high_region,true_regions(dip)))
                    score(s,dip) = bins(1);
                elseif any(ismember(high_region,neighbors))
                    score(s,dip) = bins(2);
                elseif any(ismember(high_region,neighb_neighbors))
                    score(s,dip) = bins(3);
                else
                    score(s,dip) = bins(4);
                end
            end
        end

        bins = bins(end:-1:1);
        count = hist(score(:),bins);
        figure;bar(4:-1:1,count)
        xticklabels({"correct","neighbor","2nd neighbor","wrong"})
        ylabel('Number of regions')
        title('reconstructed regions distribution')
        ax = gca;
        ax.TitleFontSizeMultiplier = 1.5;

        %boxplot
        figure;boxplot(score(:))
        title(['mean score = ', num2str(mean(score(:)))])
        ax = gca;
        ax.TitleFontSizeMultiplier = 1.5;
        yticks(bins)
        xticklabels('')


%         % (uncomment) show ROIs on the cortex
%         s = 3; %set the desired session
%         m = zeros(size(atlas.tissue));
%         m(ismember(atlas.tissue,high_region)) = 2; %reconstructed region
%         m(ismember(atlas.tissue,atlas.tissue(dipole_idx(s,:)))) = 1; %true region
%         % m(ismember(atlas.tissue,63)) = 3; %overlapped region
%         figure;ft_plot_mesh(atlas,'vertexcolor',m,'facealpha',0.8);
%         view([90 90]); h = light; set(h, 'position', [0 0 0.2]); lighting gouraud; material dull
%         hold on
%         true_pos = mat2cell(atlas.pos(dipole_idx(s,:),:),3,[1,1,1]);
%         scatter3(true_pos{:},500,'r','filled')
%         colormap('jet')
%         title('reconstructed sources vs. ground truth')
%         ax = gca;
%         ax.TitleFontSizeMultiplier = 1.5;
end