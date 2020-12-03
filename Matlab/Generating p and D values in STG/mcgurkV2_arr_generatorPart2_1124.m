
%% Cong compare
SPECS.PLOTaxis = [-.5 .5];
%SPECS.plot_conds = [1 2 3]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
SPECS.plot_conds = [102 202 302]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
%SPECS.plot_conds = [2]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
% SPECS.plot_chans = [20 28 30 31 37 38 39 40 42 43 44]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
%SPECS.plot_chans = [1:47]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
SPECS.plot_chans = [electrodes_to_analyze]';
SPECS.ERPbaseline = [-1 -.5]; %[-.05 0]
SPECS.Powerbaseline = [-1 -.5];
SPECS.Freq = [1 47.5];%[5 40];
SPECS.FreqSpacing = 0;
SPECS.FreqWidth = [1.5];%2.5;
SPECS.FreqCycles = [3 20];%[3];
SPECS.Normalization = 4; %0 = No Normalization, 1 = Log Norming, 2 = Z AVG, 3 = dB AVG, 4 = Z-Norm Indiv Trials/Freqs
SPECS.StatsType = 1; % 1 = Parametric. 2 = Non-Parametric. Parametric only works on No-Norm or Whitened Data
SPECS.NPerm=200; % Num permutations when needed
SPECS.PThresh = .01; % PVal threshold for stat comparisons. -1 for no stats
SPECS.PTail = 2; % 1 or 2 tailed tests
SPECS.ErrBars = 2; % 1 = SEM, 2 = One-Sample Conf Intervals
SPECS.MultCompType = 1; % 0 = None, 1 = FDR, 2 = MaxStat, 3 = Cluster, 4 = Bonferroni correction, 5 = Mean Amplitude (MeanA). 6 = Peak Amp (PKA)
%SPECS.Analysis_Windows = [-.5 .5]; % SPECS.Analysis_Windows = [.2 .3;.5 .6]; If empty, use full xaxis time-range. Otherwise perform MultCompType within limited time range. Req for MultCompType=5.
SPECS.PeakPolarity = 1; % Neg 1 (down) or Pos 1 (up). Only for Peak Amp analysis
SPECS.GridLines = 0; % Time in Secs between grid lines. 0 for no grid
SPECS.FigSize = 0; %1 for full screen, 0 for normal, -1 for don't display
SPECS.Replot =0; % 1 to simply replot last figure set. 0 to run analysis like normal
SPECS.types = {'ERP'};% SPECS.types = {'ERP','ERSP-line','ERSP','ERSP-individ'};
ECOG_plotter_v30c_testing 

close all
p_values = min(AVG.ERPpvalues_corr_vals,[],2);
[~, d_idx_max] = max(abs(AVG.ERP_EffectSize(:,an_window_start:an_window_end)),[],2);
d_values = arrayfun(@(n) AVG.ERP_EffectSize(n,an_window_start+d_idx_max(n)-1), 1:length(d_idx_max))';

%d_values = max(abs(AVG.ERP_EffectSize(:, 1501:2501)),[],2);

%d_values = arrayfun(@(n) max(abs(AVG.ERP_EffectSize(n,:))), 1:num_of_channels)';

arr(elec_locs,12) = num2cell(p_values(elec_locs));
arr(elec_locs,13) = num2cell(d_values(elec_locs));

%% Aud compare
SPECS.PLOTaxis = [-.5 .5];
%SPECS.plot_conds = [4 5 6]';
SPECS.plot_conds = [101 201 301]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
%SPECS.plot_conds = [2]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
% SPECS.plot_chans = [20 28 30 31 37 38 39 40 42 43 44]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
%SPECS.plot_chans = [1:47]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
SPECS.plot_chans = [electrodes_to_analyze]';
SPECS.ERPbaseline = [-1 -.5]; %[-.05 0]
SPECS.Powerbaseline = [-1 -.5];
SPECS.Freq = [1 47.5];%[5 40];
SPECS.FreqSpacing = 0;
SPECS.FreqWidth = [1.5];%2.5;
SPECS.FreqCycles = [3 20];%[3];
SPECS.Normalization = 4; %0 = No Normalization, 1 = Log Norming, 2 = Z AVG, 3 = dB AVG, 4 = Z-Norm Indiv Trials/Freqs
SPECS.StatsType = 1; % 1 = Parametric. 2 = Non-Parametric. Parametric only works on No-Norm or Whitened Data
SPECS.NPerm=200; % Num permutations when needed
SPECS.PThresh = .01; % PVal threshold for stat comparisons. -1 for no stats
SPECS.PTail = 2; % 1 or 2 tailed tests
SPECS.ErrBars = 2; % 1 = SEM, 2 = One-Sample Conf Intervals
SPECS.MultCompType = 1; % 0 = None, 1 = FDR, 2 = MaxStat, 3 = Cluster, 4 = Bonferroni correction, 5 = Mean Amplitude (MeanA). 6 = Peak Amp (PKA)
%SPECS.Analysis_Windows = [-.5 .5]; % SPECS.Analysis_Windows = [.2 .3;.5 .6]; If empty, use full xaxis time-range. Otherwise perform MultCompType within limited time range. Req for MultCompType=5.
SPECS.PeakPolarity = 1; % Neg 1 (down) or Pos 1 (up). Only for Peak Amp analysis
SPECS.GridLines = 0; % Time in Secs between grid lines. 0 for no grid
SPECS.FigSize = 0; %1 for full screen, 0 for normal, -1 for don't display
SPECS.Replot =0; % 1 to simply replot last figure set. 0 to run analysis like normal
SPECS.types = {'ERP'};% SPECS.types = {'ERP','ERSP-line','ERSP','ERSP-individ'};
ECOG_plotter_v30c_testing 

close all
p_values = min(AVG.ERPpvalues_corr_vals,[],2);
[~, d_idx_max] = max(abs(AVG.ERP_EffectSize(:,an_window_start:an_window_end)),[],2);
d_values = arrayfun(@(n) AVG.ERP_EffectSize(n,an_window_start+d_idx_max(n)-1), 1:length(d_idx_max))';

%d_values = max(abs(AVG.ERP_EffectSize(:, 1501:2501)),[],2);

%d_values = arrayfun(@(n) max(abs(AVG.ERP_EffectSize(n,:))), 1:num_of_channels)';

arr(elec_locs,14) = num2cell(p_values(elec_locs));
arr(elec_locs,15) = num2cell(d_values(elec_locs));

%% Vis compare
SPECS.PLOTaxis = [-.5 .5];
%SPECS.plot_conds = [7 8 9]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
SPECS.plot_conds = [104 204 304]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
%SPECS.plot_conds = [2]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
% SPECS.plot_chans = [20 28 30 31 37 38 39 40 42 43 44]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
%SPECS.plot_chans = [1:47]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
SPECS.plot_chans = [electrodes_to_analyze]';
SPECS.ERPbaseline = [-1 -.5]; %[-.05 0]
SPECS.Powerbaseline = [-1 -.5];
SPECS.Freq = [1 47.5];%[5 40];
SPECS.FreqSpacing = 0;
SPECS.FreqWidth = [1.5];%2.5;
SPECS.FreqCycles = [3 20];%[3];
SPECS.Normalization = 4; %0 = No Normalization, 1 = Log Norming, 2 = Z AVG, 3 = dB AVG, 4 = Z-Norm Indiv Trials/Freqs
SPECS.StatsType = 1; % 1 = Parametric. 2 = Non-Parametric. Parametric only works on No-Norm or Whitened Data
SPECS.NPerm=200; % Num permutations when needed
SPECS.PThresh = .01; % PVal threshold for stat comparisons. -1 for no stats
SPECS.PTail = 2; % 1 or 2 tailed tests
SPECS.ErrBars = 2; % 1 = SEM, 2 = One-Sample Conf Intervals
SPECS.MultCompType = 1; % 0 = None, 1 = FDR, 2 = MaxStat, 3 = Cluster, 4 = Bonferroni correction, 5 = Mean Amplitude (MeanA). 6 = Peak Amp (PKA)
%SPECS.Analysis_Windows = [-.5 .5]; % SPECS.Analysis_Windows = [.2 .3;.5 .6]; If empty, use full xaxis time-range. Otherwise perform MultCompType within limited time range. Req for MultCompType=5.
SPECS.PeakPolarity = 1; % Neg 1 (down) or Pos 1 (up). Only for Peak Amp analysis
SPECS.GridLines = 0; % Time in Secs between grid lines. 0 for no grid
SPECS.FigSize = 0; %1 for full screen, 0 for normal, -1 for don't display
SPECS.Replot =0; % 1 to simply replot last figure set. 0 to run analysis like normal
SPECS.types = {'ERP'};% SPECS.types = {'ERP','ERSP-line','ERSP','ERSP-individ'};
ECOG_plotter_v30c_testing 

close all
p_values = min(AVG.ERPpvalues_corr_vals,[],2);
[~, d_idx_max] = max(abs(AVG.ERP_EffectSize(:,an_window_start:an_window_end)),[],2);
d_values = arrayfun(@(n) AVG.ERP_EffectSize(n,an_window_start+d_idx_max(n)-1), 1:length(d_idx_max))';

%d_values = max(abs(AVG.ERP_EffectSize(:, 1501:2501)),[],2);

%d_values = arrayfun(@(n) max(abs(AVG.ERP_EffectSize(n,:))), 1:num_of_channels)';

arr(elec_locs,16) = num2cell(p_values(elec_locs));
arr(elec_locs,17) = num2cell(d_values(elec_locs));

electrodes_in_STG = arr;

elec_stg_dir = 'C:\Users\djbranglab-admin\Desktop\Karthik\IMRF\Electrodes in STG';
save_dir = [elec_stg_dir,'\',subjid];
if(strcmp(analysis_window,'aud-onset'))
    mkdir(save_dir,'audonset')
    save([save_dir '\audonset\' 'Elec_STG.mat'],'electrodes_in_STG')
elseif(strcmp(analysis_window,'pre-stim'))
    mkdir(save_dir,'\prestim')
    save([save_dir '\prestim\' 'Elec_STG.mat'],'electrodes_in_STG')
elseif(strcmp(analysis_window,'face-onset'))
    mkdir(save_dir,'faceonset')
    save([save_dir '\faceonset\' 'Elec_STG.mat'],'electrodes_in_STG')
elseif(strcmp(analysis_window,'face-move'))    
    mkdir(save_dir,'facemove')
    save([save_dir '\facemove\' 'Elec_STG.mat'],'electrodes_in_STG')
end

save([save_dir '\Electrodes\' 'Electrode_Labels.mat'],'Electrode_Labels')

