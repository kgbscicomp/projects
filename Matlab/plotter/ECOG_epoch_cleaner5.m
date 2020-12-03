% ECOG_epoch_cleaner
% 
% RejThresh = 3;
% Induced_chan_var_array = zeros(size(hband_sig_morlet,1),1);
% for currchan = 1:size(hband_sig_morlet,1)
%         temp = squeeze(hband_sig_morlet(currchan,:,:));
%         a=var(temp(:));
%         Induced_chan_var_array(currchan,1) = a;
% end
% Induced_chan_bad_array = zeros(size(hband_sig_morlet,1),size(hband_sig_morlet,2));
% for currchan = 1:size(hband_sig_morlet,1)
%     mean_var = Induced_chan_var_array(currchan,1);
%     for currevent = 1:size(hband_sig_morlet,2)
%         a = var(squeeze(hband_sig_morlet(currchan,currevent,:)));
%         if a>mean_var*RejThresh
%             Induced_chan_bad_array(currchan,currevent) = 1;
%         end
%     end
% end
% 
% temp = mean2(Induced_chan_bad_array);
% fprintf('%3.1f%% rejected \n',temp*100)
% 


if isempty(SPECS.bad_trials)
    bad_trials = zeros(1,size(evoked_sig,3));
else
    bad_trials = SPECS.bad_trials;
end

% Make sure bad_trials matrix matches the number of trials in dataset

if size(evoked_sig,3) ~= length(bad_trials)
    error('bad_trials matrix does not match number of trials')
end

bad_trials = find(bad_trials);


% changed from epoch_cleaner2 to only search analyzed trials
% this way, more intuitive inference on how many trials rejected for the
% relevant conditions. still performs this search disregarding cond type
% (if more than one trial)
analyzed_trials = [];
for i=1:numevents
    if ismember(stim(i,2),SPECS.AnalyzedConds) && ~ismember(i,bad_trials)
        analyzed_trials = [analyzed_trials i];
    end
end

% for i=1:length(SPECS.plot_conds)
%     [temp] = find(SPECS.plot_conds(i)==stim(:,2));
%     analyzed_trials = [analyzed_trials temp'];
% end


% filter the data for the ERP but not the ERSP
try
    size(SPECS.filtfreq);
catch
    warning('SPECS.filtfreq not set. using defaults');
    SPECS.filtfreq=[60 120 180 240 300 360 420];
end

if ~exist('ersp_sig')
    ersp_sig=evoked_sig;
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[]; %Clear EEG before relooping
    eeglab
    EEG = pop_importdata('dataformat','array','nbchan',0,'data','evoked_sig','srate',srate,'pnts',0,'xmin',0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname','ecog','gui','off');
    
    for freqx = SPECS.filtfreq
        EEG  = pop_basicfilter( EEG, [ 1:size(evoked_sig,1)] , 'Cutoff',  freqx, 'Design', 'notch', 'Filter', 'PMnotch', 'Order',  freqx*3,'Boundary',[] );
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname','ECOG','overwrite','on','gui','off');
    end
    EEG = eeg_checkset( EEG);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    evoked_sig = EEG.data;
end

% Reject from std dev of max-min of epoch

time = epoch(1):1/srate:epoch(2);
epochtimes = dsearchn(time',CleanAxis');

Evoked_chan_var_array = zeros(size(evoked_sig,1),length(analyzed_trials));
Evoked_chan_bad_array = zeros(size(evoked_sig,1),length(analyzed_trials));

% Change to iterative process
for currchan = 1:size(evoked_sig,1)
    for currevent = 1:length(analyzed_trials)
        Evoked_chan_var_array(currchan,currevent) = std(squeeze(evoked_sig(currchan,epochtimes(1):epochtimes(2),analyzed_trials(currevent))));
    end
    
    tempbadevents = [];
    while 1
        breakloop = 0;
        tempgoodarray = 1:length(analyzed_trials);
        tempgoodarray(tempbadevents)=[];
        mean_var = mean(Evoked_chan_var_array(currchan,tempgoodarray));
        std_var = std(Evoked_chan_var_array(currchan,tempgoodarray));
        
        for currevent = tempgoodarray
            if Evoked_chan_var_array(currchan,currevent) > mean_var+(RejThreshEvoked*std_var)
                tempbadevents = [tempbadevents currevent];
            else
                breakloop=breakloop+1;
            end
        end
        if breakloop==length(tempgoodarray)
            break;
        end
    end
    Evoked_chan_bad_array(currchan,tempbadevents) = 1;
end

temp = mean2(Evoked_chan_bad_array);
fprintf('%3.1f%% rejected Evoked \n',temp*100)

close all


%then remove these events from the stim var



