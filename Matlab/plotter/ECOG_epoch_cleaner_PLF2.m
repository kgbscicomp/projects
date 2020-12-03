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

temp=size(hband_sig_buttord_wide);
Induced_chan_var_array = zeros(temp([4 3]));
for currchan = 1:size(hband_sig_buttord_wide,4)
    for currfreq=1:size(hband_sig_buttord_wide,3)
        a = [];
        for currevent = 1:size(hband_sig_buttord_wide,1)
            a(currevent) = max(abs(hband_sig_buttord_wide(currevent,:,currfreq,currchan)))-min(abs(hband_sig_buttord_wide(currevent,:,currfreq,currchan)));
        end
        Induced_chan_var_array(currchan,currfreq) = std(a);
    end
end

Induced_chan_bad_array = zeros(temp([4 3 1]));
for currchan = 1:size(hband_sig_buttord_wide,4)
    for currfreq=1:size(hband_sig_buttord_wide,3)
        mean_var = Induced_chan_var_array(currchan,currfreq);
        for currevent = 1:size(hband_sig_buttord_wide,1)
            a = max(abs(hband_sig_buttord_wide(currevent,:,currfreq,currchan)))-min(abs(hband_sig_buttord_wide(currevent,:,currfreq,currchan)));
            if a>mean_var*RejThresh
                Induced_chan_bad_array(currchan,currfreq,currevent) = 1;
            end
        end
    end
end

temp = mean2(Induced_chan_bad_array);
fprintf('%3.1f%% rejected \n',temp*100)



%then remove these events from the stim var



