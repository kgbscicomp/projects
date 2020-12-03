% ECOG_epoch_cleaner
% 
% RejThresh = 3;
% Induced_chan_var_array = zeros(size(hband_sig_buttord_wide,1),1);
% for currchan = 1:size(hband_sig_buttord_wide,1)
%         temp = squeeze(hband_sig_buttord_wide(currchan,:,:));
%         a=var(temp(:));
%         Induced_chan_var_array(currchan,1) = a;
% end
% Induced_chan_bad_array = zeros(size(hband_sig_buttord_wide,1),size(hband_sig_buttord_wide,2));
% for currchan = 1:size(hband_sig_buttord_wide,1)
%     mean_var = Induced_chan_var_array(currchan,1);
%     for currevent = 1:size(hband_sig_buttord_wide,2)
%         a = var(squeeze(hband_sig_buttord_wide(currchan,currevent,:)));
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


Induced_chan_var_array = zeros(size(hband_sig_buttord_wide,1),1);
for currchan = 1:size(hband_sig_buttord_wide,1)
    a = [];
    for currevent = 1:size(hband_sig_buttord_wide,2)
        a(currevent) = max(hband_sig_buttord_wide(currchan,currevent,:))-min(hband_sig_buttord_wide(currchan,currevent,:));
    end
    Induced_chan_var_array(currchan,1) = std(a);
end
Induced_chan_bad_array = zeros(size(hband_sig_buttord_wide,1),size(hband_sig_buttord_wide,2));
for currchan = 1:size(hband_sig_buttord_wide,1)
    mean_var = Induced_chan_var_array(currchan,1);
    for currevent = 1:size(hband_sig_buttord_wide,2)
        a = max(hband_sig_buttord_wide(currchan,currevent,:))-min(hband_sig_buttord_wide(currchan,currevent,:));
        if a>mean_var*RejThreshInduced
            Induced_chan_bad_array(currchan,currevent) = 1;
        end
    end
end

temp = mean2(Induced_chan_bad_array);
fprintf('%3.1f%% rejected Induced \n',temp*100)



%then remove these events from the stim var



