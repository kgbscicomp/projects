clear all
clc
subjid = '1187HF';
load('HF1187_McGurk-bipolar.mat');
load('C:\Users\djbranglab-admin\Desktop\Karthik\Electrode_Registration\1187HF\Electrodes\Electrode_Labels.mat');

num_of_channels = length(OUTPUT.chan_names);
% analysis_window = 'aud-onset';%options: 'pre-stim', 'face-onset', 'face-move','aud-onset'
% analysis_window = 'pre-stim';
% analysis_window = 'face-onset';
analysis_window = 'face-move';
if(strcmp(analysis_window,'aud-onset'))
    SPECS.Analysis_Windows = [0 .25]; % SPECS.Analysis_Windows = [.2 .3;.5 .6]; If empty, use full xaxis time-range. Otherwise perform MultCompType within limited time range. Req for MultCompType=5.
elseif(strcmp(analysis_window,'pre-stim'))
    SPECS.Analysis_Windows = [-1.25 -.75];
elseif(strcmp(analysis_window,'face-onset'))
    SPECS.Analysis_Windows = [-.75 -.5];
elseif(strcmp(analysis_window,'face-move'))    
    SPECS.Analysis_Windows = [-.25 0];    
end

SPECS.PLOTaxis = [-.5 .5];
SPECS.plot_conds = [7;8;9]';
%SPECS.plot_conds = [101 201 301 401]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
%SPECS.plot_conds = [1]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
% SPECS.plot_chans = [20 28 30 31 37 38 39 40 42 43 44]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
SPECS.plot_chans = [1:num_of_channels]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
%SPECS.plot_chans = [electrodes_to_analyze]';
SPECS.ERPbaseline = [-1.25 -.75]; %[-.05 0]
SPECS.Powerbaseline = [-1.25 -.75];
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
SPECS.PeakPolarity = 1; % Neg 1 (down) or Pos 1 (up). Only for Peak Amp analysis
SPECS.GridLines = 0; % Time in Secs between grid lines. 0 for no grid
SPECS.FigSize = 0; %1 for full screen, 0 for normal, -1 for don't display
SPECS.Replot =0; % 1 to simply replot last figure set. 0 to run analysis like normal
SPECS.types = {'ERP'};% SPECS.types = {'ERP','ERSP-line','ERSP','ERSP-individ'};
ECOG_plotter_v30c_testing

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%% Extracts labels of electrodes that are present in the temporal region
%% shortlists electrodes that have a p-value of less than 0.05 and are 
%% present in the temporal region 

Electrode_labels_backup = Electrode_Labels.Parc.wmparc.Bipolar.Volume; 
region_labels = Electrode_Labels.Parc.wmparc.Bipolar.Volume;
region_names = fieldnames(region_labels);
counter = 1;
for r = 1:size(region_names,1)
    region_values = extractfield(region_labels,region_names{r});
    region_values = region_values{1,1};
    stg_presence = cell(size(region_values));
    pattern = ["superiortemporal","middletemporal","transversetemporal","bankssts"];
    for i = 1:size(region_values,1)
        for j = 1:2:size(region_values,2)
            str = region_values(i,:);
            if(find(~cellfun('isempty', str(j))))
                stg_presence{i,j} = contains(str(j),pattern,'IgnoreCase',true);
                stg_presence(i,j+1) = region_values(i,j+1);
            end
        end
    end

    roi_labels = [];
    M = cellfun(@(x) double(x),stg_presence,'UniformOutput',false);
    for i = 1:size(region_values,1)
        sum_probs = 0;
        a = cell2mat(M(i,:));
        a(numel(region_values(1,:))) = 0;
        for j = 1:2:size(region_values,2)-1
            if a(j) ~= 0
                sum_probs = sum_probs + a(j+1);
            end
        end
        if(sum_probs > 0.5)
            roi_labels(size(roi_labels,1)+1,1) = i;
        end
    end
    Electrode_labels_backup.(region_names{r}) = roi_labels;
    for k = 1:size(roi_labels,1)
        elec_names{counter,1} = [region_names{r},num2str(Electrode_labels_backup.(region_names{r})(k),'%d')];
        counter = counter + 1;
    end
    
end

%% Extracts electrode number of channels that are present only in temporal region
all_chan_names = OUTPUT.chan_names';
str = all_chan_names;
pattern = strcat(elec_names,'_b');
for i = 1:size(str,1)
    a(i) = contains(str(i),pattern,'IgnoreCase',true);
end
electrodes_to_analyze = find(a);

elec_locs = electrodes_to_analyze;

if(length(pattern) ~= length(electrodes_to_analyze))
    error('Error. Non-auditory electrodes included in analysis')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Run time_lapse 2. run extract_roi after loading electrode labels
%% 3. run construct_xcel till activations 4. run phoneme_compare 5. run rest of construct_xcel 

an_window_start =(2000-(abs(SPECS.Analysis_Windows(1))*1000))+1;
an_window_end = an_window_start + abs(SPECS.Analysis_Windows(1)-SPECS.Analysis_Windows(2))*1000;

p_values = min(AVG.ERPpvalues_corr_vals,[],2);
[~, d_idx_max] = max(abs(AVG.ERP_EffectSize(:,an_window_start:an_window_end)),[],2);
d_values = arrayfun(@(n) AVG.ERP_EffectSize(n,an_window_start+d_idx_max(n)-1), 1:length(d_idx_max))';


%d_values = arrayfun(@(n) max(abs(AVG.ERP_EffectSize(n,~isnan(AVG.ERPpvalues_corr_vals(1,:))))),...
%                    1:num_of_channels)';


thresh_ind = elec_locs;

arr = cell(size(chan_names,2),4);
arr(:,1) = all_chan_names;
arr(:,2) = num2cell(p_values);
arr(:,3) = num2cell(d_values);
arr(thresh_ind,4) = num2cell(p_values(thresh_ind));  %%Electrodes with sig p values
arr(thresh_ind,5) = num2cell(d_values(thresh_ind));

arr(elec_locs,6) = num2cell(p_values(elec_locs));    %% Electrodes at STG
arr(elec_locs,7) = num2cell(d_values(elec_locs));



%% audio
SPECS.PLOTaxis = [-.5 .5];
SPECS.plot_conds = [1;2;3]';
%SPECS.plot_conds = [101 201 301 401]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
%SPECS.plot_conds = [4]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
% SPECS.plot_chans = [20 28 30 31 37 38 39 40 42 43 44]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
%SPECS.plot_chans = [1:47]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
SPECS.plot_chans = [electrodes_to_analyze]';
SPECS.ERPbaseline = [-1.25 -.75]; %[-.05 0]
SPECS.Powerbaseline = [-1.25 -.75];
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
%d_values = arrayfun(@(n) max(abs(AVG.ERP_EffectSize(n, 1501:2501))), 1:num_of_channels)';

arr(elec_locs,8) = num2cell(p_values(elec_locs));
arr(elec_locs,9) = num2cell(d_values(elec_locs));



%% Visual
SPECS.PLOTaxis = [-.5 .5];
SPECS.plot_conds = [4;5;6]';
%SPECS.plot_conds = [101 201 301 401]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
%SPECS.plot_conds = [2]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
% SPECS.plot_chans = [20 28 30 31 37 38 39 40 42 43 44]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
%SPECS.plot_chans = [1:47]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
SPECS.plot_chans = [electrodes_to_analyze]';
SPECS.ERPbaseline = [-1.25 -.75]; %[-.05 0]
SPECS.Powerbaseline = [-1.25 -.75];
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

%d_values = arrayfun(@(n) max(abs(AVG.ERP_EffectSize(n,~isnan(AVG.ERPpvalues_corr_vals(1,:))))),...
%                    1:num_of_channels)';

arr(elec_locs,10) = num2cell(p_values(elec_locs));
arr(elec_locs,11) = num2cell(d_values(elec_locs));



%% Cong compare
SPECS.PLOTaxis = [-.5 .5];
SPECS.plot_conds = [7 8 9]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
%SPECS.plot_conds = [101 201 301 401]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
%SPECS.plot_conds = [2]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
% SPECS.plot_chans = [20 28 30 31 37 38 39 40 42 43 44]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
%SPECS.plot_chans = [1:47]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
SPECS.plot_chans = [electrodes_to_analyze]';
SPECS.ERPbaseline = [-1.25 -.75]; %[-.05 0]
SPECS.Powerbaseline = [-1.25 -.75];
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
SPECS.plot_conds = [1 2 3]';
%SPECS.plot_conds = [104 204 304 404]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
%SPECS.plot_conds = [2]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
% SPECS.plot_chans = [20 28 30 31 37 38 39 40 42 43 44]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
%SPECS.plot_chans = [1:47]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
SPECS.plot_chans = [electrodes_to_analyze]';
SPECS.ERPbaseline = [-1.25 -.75]; %[-.05 0]
SPECS.Powerbaseline = [-1.25 -.75];
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
SPECS.plot_conds = [4 5 6]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
%SPECS.plot_conds = [102 202 302 402]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
%SPECS.plot_conds = [2]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
% SPECS.plot_chans = [20 28 30 31 37 38 39 40 42 43 44]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
%SPECS.plot_chans = [1:47]'; % if more than one column, channels concatenated, e.g., [1 2;5 7] combines channels 1 and 2. Put a 0 to combine uneven combinations of channels (e.g., 1 and 2, and 5 by itself)
SPECS.plot_chans = [electrodes_to_analyze]';
SPECS.ERPbaseline = [-1.25 -.75]; %[-.05 0]
SPECS.Powerbaseline = [-1.25 -.75];
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

