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


% changed from epoch_cleaner2 to only search analyzed trials
% this way, more intuitive inference on how many trials rejected for the
% relevant conditions. still performs this search disregarding cond type
% (if more than one trial)
analyzed_trials = [];
for i=1:numevents
    if ismember(stim(i,2),SPECS.plot_conds)
        analyzed_trials = [analyzed_trials i];
    end
end

% for i=1:length(SPECS.plot_conds)
%     [temp] = find(SPECS.plot_conds(i)==stim(:,2));
%     analyzed_trials = [analyzed_trials temp'];
% end

% Reject from std dev of max-min of epoch

Evoked_chan_var_array = zeros(size(evoked_sig,1),1);
time = epoch(1):1/srate:epoch(2);
% epochtimes = length(time);%srate*1:srate*3;
epochtimes = dsearchn(time',SPECS.PLOTaxis');
for currchan = 1:size(evoked_sig,1)
        temp = squeeze(evoked_sig(currchan,epochtimes,analyzed_trials));
        a=std(temp(:));
        Evoked_chan_var_array(currchan,1) = a;
end
Evoked_chan_bad_array = zeros(size(evoked_sig,1),length(analyzed_trials));
for currchan = 1:size(evoked_sig,1)
    mean_var = Evoked_chan_var_array(currchan,1);
    for currevent = 1:length(analyzed_trials)
        a = std(squeeze(evoked_sig(currchan,epochtimes,currevent)));
        if a>mean_var*RejThreshEvoked
            Evoked_chan_bad_array(currchan,currevent) = 1;
        end
    end
end
temp = mean2(Evoked_chan_bad_array);
fprintf('%3.1f%% rejected Evoked \n',temp*100)

ERSP_chan_var_array = zeros(size(evoked_sig,1),1);
ERSPtime = [SPECS.PLOTaxis(1)-.5 SPECS.PLOTaxis(2)+.5];
epochtimes = dsearchn(time',ERSPtime');
for currchan = 1:size(evoked_sig,1)
        temp = squeeze(evoked_sig(currchan,epochtimes,analyzed_trials));
        a=std(temp(:));
        ERSP_chan_var_array(currchan,1) = a;
end
ERSP_chan_bad_array = zeros(size(evoked_sig,1),length(analyzed_trials));
for currchan = 1:size(evoked_sig,1)
    mean_var = ERSP_chan_var_array(currchan,1);
    for currevent = 1:length(analyzed_trials)
        a = std(squeeze(evoked_sig(currchan,epochtimes,currevent)));
        if a>mean_var*RejThreshEvoked
            ERSP_chan_bad_array(currchan,currevent) = 1;
        end
    end
end
temp = mean2(ERSP_chan_bad_array);
fprintf('%3.1f%% rejected ERSP \n',temp*100)




%then remove these events from the stim var



