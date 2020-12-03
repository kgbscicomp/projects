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

% Reject from std dev of max-min of epoch


Induced_chan_var_array = zeros(size(hband_data_all,1),1);
tempdata = squeeze(mean(hband_data_all,2));
for currchan = 1:size(hband_data_all,1)
    a = [];
    for currevent = 1:size(hband_data_all,4)
        a(currevent) = max(tempdata(currchan,:,currevent))-min(tempdata(currchan,:,currevent));
    end
    Induced_chan_var_array(currchan,1) = std(a);
end
Induced_chan_bad_array = zeros(size(hband_data_all,1),size(hband_data_all,4));
for currchan = 1:size(hband_data_all,1)
    mean_var = Induced_chan_var_array(currchan,1);
    for currevent = 1:size(hband_data_all,4)
        a = max(tempdata(currchan,:,currevent))-min(tempdata(currchan,:,currevent));
        if a>mean_var*RejThreshInduced
            Induced_chan_bad_array(currchan,currevent) = 1;
        end
    end
end

temp = mean2(Induced_chan_bad_array);
fprintf('%3.1f%% rejected Induced \n',temp*100)





Evoked_chan_var_array = zeros(size(evoked_sig,1),1);
for currchan = 1:size(evoked_sig,1)
        temp = squeeze(evoked_sig(currchan,:,:));
        a=std(temp(:));
        Evoked_chan_var_array(currchan,1) = a;
end
Evoked_chan_bad_array = zeros(size(evoked_sig,1),size(evoked_sig,3));
for currchan = 1:size(evoked_sig,1)
    mean_var = Evoked_chan_var_array(currchan,1);
    for currevent = 1:size(evoked_sig,3)
        a = std(squeeze(evoked_sig(currchan,:,currevent)));
        if a>mean_var*RejThreshEvoked
            Evoked_chan_bad_array(currchan,currevent) = 1;
        end
    end
end


temp = mean2(Evoked_chan_bad_array);
fprintf('%3.1f%% rejected Evoked \n',temp*100)





%then remove these events from the stim var



