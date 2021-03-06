% % % % % function ECOG_plotter_v30h_testing
%% Example
% tic
% SPECS.plot_chans =[1]';
% SPECS.PLOTaxis = [-1 1];
% SPECS.plot_conds = [1 2]'; % if more than one column, stim concatenated, e.g., [100 101;500 500] combines stim from 100 and 101. Put a 0 to test uneven combinations of conditions (e.g., Aud & Vis vs. Aud-Vis)
% SPECS.ERPbaseline = [-1 0]; %[-.05 0]
% SPECS.Powerbaseline = [-1 -.1];
% SPECS.Freq = [70 150];
% SPECS.FreqWidth = 5; % Distance in frequencies between wavelet center frequencies
% SPECS.FreqSpacing = 0; % 0 = Linear spacing. 1 = Log spacing.
% SPECS.FreqCycles = [20];%[3]; Will always do a min of 3 cycles
% SPECS.Normalization = 2; % 0 = No Normalization, 1 = dB change relative to baseline, 2 = 1/f norm (RECOMMENDED), 3 = Percent Change (uses median)
% SPECS.StatsType = 1; % 1 = Parametric. 2 = Non-Parametric. Parametric only works on No-Norm or Whitened Data
% SPECS.NPerm=200; % Num permutations when needed
% SPECS.PThresh = .05; % PVal threshold for stat comparisons. -1 for no stats
% SPECS.PTail = 2; % 1 or 2 tailed tests
% SPECS.ErrBars = 2; % 1 = SEM, 2 = One-Sample Conf Intervals
% SPECS.MultCompType = 1; % 0 = None, 1 = FDR, 2 = MaxStat, 3 = Cluster, 4 = Bonferroni correction, 5 = Mean time range
% SPECS.Analysis_Windows = []; % SPECS.Analysis_Windows = [.2 .3;.5 .6]; If empty, use full xaxis time-range. Otherwise perform MultCompType within limited time range. Req for MultCompType=5.
% SPECS.RemoveOutliers = 1; % 0 = No removal (Only those rejected during pre-processing). 1 = Replace w/ Clipped Values (RECOMMENDED). 2 = Replace w/ NaN. 3 = Replace w/ Median. 4. Include all data, including previously rejected (use with classification)
% SPECS.OutlierMethod = 1; % 1 = 3 scaled MAD (RECOMMENDED). 2 = 3 SD. 3 = 1.5 IQR. 4 = Grubbs. 5 = gesd. 6 = 3 x cooks. 7 = sn. 8 = qn
% SPECS.OrthogonalPreSelection = 0; % 0 = No (Default), 1 = Yes. Only analyze electrodes that have significant activity regardless of condition.
% SPECS.OrthogonalPreSelection_PThresh = .01; % .05 (Default). Mult comparison corrected threshold for inclusion.
% SPECS.OrthogonalPreSelection_Polarity = 1; % 0 = Both sides. 1=Positive sig only. 2=Negative Sig Only. HGp is usually going to be positive sig only.
% SPECS.GridLines = 0; % Time in Secs between grid lines. 0 for no grid
% SPECS.FigSize = 0; %1 for full screen, 0 for normal, -1 for don't display
% SPECS.Replot =0; % 1 to simply replot last figure set. 0 to run analysis like normal
% SPECS.types = {'ERSP-line'};% SPECS.types = {'ERP','ERSP-line','ERSP','ERSP-individ'};
% ECOG_plotter_v30h_testing









%%

% This function takes the output of the function ECOG_filter_analysis7
% along with the following inputs to plot the processed ECoG data:
%
% SPECS.PLOTaxis = [-1 1];
% SPECS.plot_conds = [1 2;3 3]; % Plot 1-4 conditions, with ordered colors of red, green, blue, orange
% % if more than one column, columns concatenated, e.g., [100 101;500 500] combines stim from 100 and 101
% RejThresh = 3;
% ECOG_epoch_cleaner3

% SPECS.bad_trials = []; % Bad trials by trial number (e.g. for blinks)
% SPECS.ERPbaseline = [-.5 0]; %[-.05 0]
% SPECS.Powerbaseline = [-1 -.5];
% SPECS.xticksize = .25;
% SPECS.plot_chans = [33];% 34 46:47];%1:length(chan_names);
% SPECS.Freq = [3 40];%[5 40];
% SPECS.FreqWidth = [2];%2.5;
% SPECS.FreqCycles = [3 10];%[3];
% SPECS.Normalization = 1; %0 = No Normalization, 1 = Whitening, 2 = Z, 3 = dB
% SPECS.StatsType = 1; % 1 = Parametric. 2 = Non-Parametric. Parametric only works on No-Norm or Whitened Data
% SPECS.NPerm=100; % Num permutations when needed
% SPECS.PThresh = .05; % PVal threshold for stat comparisons
% SPECS.PTail = 2; % 1 or 2 tailed tests
% SPECS.MultCompType = 1; % FDR, MaxStat, Cluster
% SPECS.MultCompPoints = 1; % Number of datapoints averaged in each bin prior to conducting statistical correction
% SPECS.FigSize = 1; %1 for full screen, 0 for normal, -1 for don't display
% SPECS.types = {'ERP','ERSP'};% SPECS.types = {'ERP','ERSP-line','ERSP','ERSP-individ','ITPC-line','ITPC','R2-line','R2','Distribution'};
% ECOG_plotter

% Oct 2016
% Added in feature to plot distribution of activity in a range across conditions. Useful for AudLoc_v2.
% Also increased the number of channels that could be compared. If more than 4 skipping stats

% April 2017
% Rewriting of most of the code.
%
% Added new normalization techniques. Whitening and No Norm will allow
% parametric analysis of the data and is similar to how Parvizi's group's
% preprocessing pipeline.

% 9/15/2017 - DB
% Removed detrend function since should be done as a high-pass filter in
% intial pre-processing stage.
% Add check for whether functions were changed in the pre-processing script
% relative to this version

% 10/11/2017 - DB
% Added in 4th option for normalization
% znorms each freq bin across the entire timerange
% piloting shows less variability that log norming

% 10/18/2017 - DB
% Fixed the fixed epoch variable from -2 2. Now takes input
% Updated the ECOG_plotter_dep_parametric_statistics_v8 to perform analyses
% based on time-windows

% 10/27/2017 - DB - Not implemented yet
% ROI Analysis
% Only works with one-sample tests to prevent hunting for sig between conds
% Will generate .mat file with results for easy access for future plots and documentation
% Give it 2 params for each input of interest as search params
% SPECS.Freq, will narrow down (e.g., 60-250 will become 70-150). Will search in step sizes set by SPECS.FreqWidth
% SPECS.FreqCycles will be harder to optimize, particularly if already have 2 inputs. Omit for now.
% SPECS.Analysis_Windows could eventually be tested but not implemented yet.
% Will also save array of sig versus non-sig electrodes

% 11/6/2017 - DB
% log-spacing of frequencies
% when given 2 inputs to SPECS.FreqCycles, cycles will be rounded.

% 2/27/19 - DB
% Fixed ERSP processing for non-symetrical epochs

% 12/22/19 - DB
% 1. Added percent change option
% 2. Added trial rejection
% 3. Added orthogonal option in which if multiple conditions are present,
%    the script will first preselect electrodes for analysis based on them
%    having significant activity regardless of condition type


% 2/9/20 - DB
% Updated the way that outliers and processed and added additional outlier
% detection methods



%% Updated with changes of the pre-processing scripts as needed
% matching versions of the pre-processing scripts
% if changes there don't affect changes here, can have mult options
Expected_Versions = [10 11]; 

% if OUTPUT.reference==2 && OUTPUT.Script_Version==9
%     error('There was a problem with bipolar referencing in Filter Analysis v9. Reprocess with newer script.')
% end

%% Check version of pre-processing
if isfield(OUTPUT, 'Script_Version')==1
    if ~ismember(OUTPUT.Script_Version,Expected_Versions)
        error(['This script was designed to be used with ECOG_filter_analysis_v' num2str(Expected_Versions(end))]);
    end
else
    error(['This script was designed to be used with ECOG_filter_analysis_v' num2str(Expected_Versions(end))]);
end


%% Prep all the plot information
if SPECS.Replot==0
    plot_conds_matrix=SPECS.plot_conds';
    plot_conds = plot_conds_matrix(1,:);
    plot_chans=SPECS.plot_chans';
    SPECS.Options.types = {'ERP','ERSP-line','ERSP','ERSP-individ','ITPC-line','ITPC','R2-line','R2','Distribution','RawSpect-line','RawSpect','ITPA','Variance','Variance-line'};
    SPECS.Options.plots = [1 1 .5 1 1 .5 1 .5 1 1 .5 .5 .5 1]; % how many of two rows will each plot take in the subplot
    SPECS.Options.ylabel = {'microVolts','Normalized Amplitude','Normalized Amplitude','Normalized Amplitude','Phase-Consistency (0-1)','Phase-Consistency (0-1)','Phase-Consistency (0-1)','Phase-Consistency (0-1)','Mean Normalized Amplitude','Normalized Amplitude','Normalized Amplitude','Phase-Consistency (0-1)','Coefficient of Variation \sigma/\mu','Coefficient of Variation \sigma/\mu'};
    tempPLOTtypes = SPECS.types;
    numanalyses = length(tempPLOTtypes);
    LineThickness = 2.5;
    
    bad_trials = OUTPUT.bad_trials;
    time = OUTPUT.time;
    srate = OUTPUT.srate;
    trial_length = size(OUTPUT.evoked_sig,2);
    evoked_sig = OUTPUT.evoked_sig;
    ersp_sig = OUTPUT.ersp_sig;
    CondNums = OUTPUT.CondNums;
    Conditions = OUTPUT.Conditions;
    stim = OUTPUT.stim;
    chan_names = OUTPUT.chan_names;
    Evoked_chan_bad_array = OUTPUT.Evoked_chan_bad_array;
    oneoverf_array = OUTPUT.oneoverf;
    
    if ~ismember(SPECS.StatsType,[1 2])
        error('SPECS.StatsType must be 1 or 2.')
    end
    
    
    % Setup outlier rej options
    % SPECS.RemoveOutliers
    % 0 = No removal
    % 1 = Replace w/ Clipped Values (RECOMMENDED)
    % 2 = Replace w/ NaN
    % 3 = Replace w/ Median
    % 4 = No removal and clear previous removal from pre-processing
    RemoveOutliers = SPECS.RemoveOutliers;
    if RemoveOutliers==0 || RemoveOutliers==4
        RemoveOutliers = 'NONE';
    elseif RemoveOutliers==1
        RemoveOutliers = 'clip';
    elseif RemoveOutliers==2
        RemoveOutliers = 'NaN';
    elseif RemoveOutliers==3
        RemoveOutliers = 'median';
    end
    
    % SPECS.OutlierMethod
    % 1 = 3 scaled MAD (RECOMMENDED)
    % 2 = 3 SD
    % 3 = 1.5 IQR
    % 4 = Grubbs
    % 5 = gesd
    % 6 = 3 x cooks
    % 7 = sn
    % 8 = qn
    OutlierMethod = SPECS.OutlierMethod;
    if OutlierMethod==1
        OutlierMethod = 'median';
    elseif OutlierMethod==2
        OutlierMethod = 'mean';
    elseif OutlierMethod==3
        OutlierMethod = 'quartiles';
    elseif OutlierMethod==4
        OutlierMethod = 'grubbs';
    elseif OutlierMethod==5
        OutlierMethod = 'gesd';
    elseif OutlierMethod==6
        OutlierMethod = 'cooks';
    elseif OutlierMethod==7
        OutlierMethod = 'sn';
    elseif OutlierMethod==8
        OutlierMethod = 'qn';
    end
    
    % Zero out Evoked_chan_bad_array and check on indiv conds
    if ~strcmp(RemoveOutliers,'NONE') || SPECS.RemoveOutliers==4
       Evoked_chan_bad_array = Evoked_chan_bad_array*0;
    end
    
    % % % % % Parametric = SPECS.Parametric;
    
    PLOT_Specs.plotnum = 0;
    PLOT_Specs.rowcount = size(plot_conds_matrix,2);
    PLOT_Specs.colcount = 0;
    PLOT_Specs.halfwhole = [];
    for counter=1:numanalyses
        PLOT_Specs.colcount=PLOT_Specs.colcount+1;
        temp = SPECS.Options.plots(find(strcmp(tempPLOTtypes(counter),SPECS.Options.types)==1));
        if isempty(temp), error('One or more SPECS.types are invalid.');end
        if temp == .5, PLOT_Specs.plotnum = PLOT_Specs.plotnum+PLOT_Specs.rowcount;PLOT_Specs.halfwhole(counter)=.5;
        else PLOT_Specs.plotnum = PLOT_Specs.plotnum+1;PLOT_Specs.halfwhole(counter)=1;end
    end
    PLOT_Specs.arraycount = PLOT_Specs.colcount*PLOT_Specs.rowcount;
    PLOT_Specs.plotarray = zeros(PLOT_Specs.rowcount,PLOT_Specs.colcount);
    counter=1;
    for j=1:PLOT_Specs.rowcount
        for i=1:PLOT_Specs.colcount
            PLOT_Specs.plotarray(j,i)=counter;
            counter=counter+1;
        end
    end
    PLOT_Specs.plotarray=PLOT_Specs.plotarray(:)';
    PLOT_Specs.plotlocs = zeros(PLOT_Specs.plotnum,PLOT_Specs.rowcount);
    PLOT_Specs.types = {};
    counter1=1;
    counter2=1;
    for i=1:PLOT_Specs.colcount
        if PLOT_Specs.halfwhole(i)==.5
            for j=1:PLOT_Specs.rowcount
                PLOT_Specs.types(counter2)=tempPLOTtypes(i);
                PLOT_Specs.plotlocs(counter2,:) = ones(1,PLOT_Specs.rowcount)*PLOT_Specs.plotarray(counter1);
                counter1=counter1+1;
                counter2=counter2+1;
            end
        else
            temp = counter1:1:counter1+PLOT_Specs.rowcount-1;
            PLOT_Specs.types(counter2)=tempPLOTtypes(i);
            PLOT_Specs.plotlocs(counter2,:) =PLOT_Specs.plotarray(temp);
            counter1=counter1+PLOT_Specs.rowcount;
            counter2=counter2+1;
        end
    end
    
    
    % YLIM = SPECS.Yaxis;
    SPECS.Options.ERP = 1; % All processes need ERP Preprocessing
    
    % If anything but ERP, will need spectral
    if length(PLOT_Specs.types)>1 || isempty(strmatch('ERP', SPECS.types))==1
        SPECS.Options.ERSP = 1;
    else
        SPECS.Options.ERSP=0; 
    end
        
    if SPECS.Options.ERSP==1
        if isempty(SPECS.Powerbaseline)
             disp('Warning. ERSP Baseline not being used')
        else
            if SPECS.Powerbaseline(1)<SPECS.PLOTaxis(1) || SPECS.Powerbaseline(2)>SPECS.PLOTaxis(2)
                disp('Warning. ERSP Baseline is outside of plotted range')
            end
        end
    end
    
    % Enter condition numbers interested in plotting
    xaxis = dsearchn(time',SPECS.PLOTaxis');
    xaxis = xaxis(1):xaxis(2);
    zerotime=dsearchn(time',0);
    
    % Create the baselines by default    
    % ERP Baseline    
    plotbaseline = dsearchn(time',SPECS.ERPbaseline');
    if plotbaseline(1)>=plotbaseline(2)
        error('First value in SPECS.ERPbaseline needs to be before second value')
    end
    plotbaseline = plotbaseline(1):plotbaseline(2);

    % ERSP Baseline
    if isempty(SPECS.Powerbaseline)
        zbaseline = [];
    else
        zbaseline = dsearchn(time',SPECS.Powerbaseline');
        zbaseline = zbaseline(1):zbaseline(2);
    end
    
    % If SPECS.Analysis_Windows is empty or doesn't exist, use xaxis
    if isfield(SPECS, 'Analysis_Windows')==0
        SPECS.Analysis_Windows = [time(xaxis(1)) time(xaxis(end))];
    elseif isempty(SPECS.Analysis_Windows)
        SPECS.Analysis_Windows = [time(xaxis(1)) time(xaxis(end))];
    end
    
    % Iterate through each analysis window
    clear analysis_windows
    for i=1:size(SPECS.Analysis_Windows,1)
        temp = dsearchn(time',SPECS.Analysis_Windows(i,:)');
        analysis_windows.(genvarname(char(['TIMES',num2str(i)]))) = temp(1):temp(2);
    end
    
    % Include warning if the SPECS.Analysis_Windows overlap
    if size(SPECS.Analysis_Windows,1)>1
    for i=1:size(SPECS.Analysis_Windows,1)-1
        for j=i+1:size(SPECS.Analysis_Windows,1)
            if isempty(intersect(analysis_windows.(genvarname(char(['TIMES',num2str(i)]))),analysis_windows.(genvarname(char(['TIMES',num2str(j)])))))==0
                warndlg('SPECS.Analysis_Windows Overlap. Overlapping time periods will be overwritten by subsequent window analysis','!! Warning !!')
            end
        end
    end
    end
    
    % Check if SPECS.PeakPolarity exists
    if isfield(SPECS, 'PeakPolarity')==0
        SPECS.PeakPolarity = [];
    end
    if ismember(SPECS.MultCompType,6:7)==1 && isempty(SPECS.PeakPolarity)
        error('SPECS.PeakPolarity not set. Need to set the polarity of the peak of interest');
    end
    PeakPolarity = SPECS.PeakPolarity;
    
    PThresh = SPECS.PThresh;
    PTail = SPECS.PTail;
    if ~ismember(PTail,1:2)
        error('SPECS.PTail requires entry of 1 or 2 tails');
    end
    
    if SPECS.FigSize==1,FullSrnFig=1;else FullSrnFig=0;end
    NPerm=SPECS.NPerm;
    NumConds = PLOT_Specs.rowcount;
    
    % colorarray = {'r','g','b','c'};%old
    % colorarray=[1 0 0;0 0 1;0 1 0;0 1 1; 1 1 0; 1 0 1; 1 .5 0; 0 1 .5;0 .5 1;1 0 .5;.5 1 0;.5 0 1;.5 .5 .5; 0 0 0]; % red, green, blue, cyan, yellow, magenta, orange, blue-green, light blue,fuscia,lime-green,purple, grey,black
    colorarray=[190/255 32/255 38/255;39/255 126/255 183/255;59/255 181/255 74/255; 251/255 191/255 21/255; .188 .188 .188;1 0 1; 1 .5 0; 0 1 .5;0 .5 1;1 0 .5;.5 1 0;.5 0 1;.5 .5 .5; 0 0 0]; % red, green, blue, cyan, yellow, magenta, orange, blue-green, light blue,fuscia,lime-green,purple, grey,black
    %colorarrayShade = {'r','b','g','c'};
        
    % Log or linear spacing of freqs
    if isfield(SPECS,'FreqSpacing')
        if SPECS.FreqSpacing==1
            loglin_freq = 1; % log
        else
            loglin_freq = 0; % linear
        end
    else
        loglin_freq = 0; % linear
    end
        
    % define wavelet parameters
    if SPECS.Options.ERSP==1
        Freq = SPECS.Freq;
        FreqWidth = SPECS.FreqWidth;
        FreqCycles = SPECS.FreqCycles;
        if length(Freq)==1
            frex = Freq(1);
        elseif length(Freq)==2 && loglin_freq==0
            frex = Freq(1):FreqWidth:Freq(2);
            % frex = linspace(Freq(1),Freq(2),FreqBins);
        elseif length(Freq)==2 && loglin_freq==1
            frex = logspace(log10(Freq(1)),log10(Freq(2)),length(Freq(1):FreqWidth:Freq(2)));
        else
            error('frequencies not properly defined')
        end
            
        num_frex = length(frex);
        if length(FreqCycles)==1 % variable number of cycles to equate time length
            % cyclex = [FreqCycles(1) round(FreqCycles(1).*(frex(2:end)./frex(1)))];
            cyclex = [FreqCycles(1) (FreqCycles(1).*(frex(2:end)./frex(1)))];
        elseif length(FreqCycles)==2
            if loglin_freq==0
                % cyclex = round(linspace(FreqCycles(1),FreqCycles(2),num_frex));
                cyclex = (linspace(FreqCycles(1),FreqCycles(2),num_frex));
            elseif loglin_freq==1
                % cyclex  = round(logspace(log10(FreqCycles(1)),log10(FreqCycles(2)),num_frex));
                cyclex  = (logspace(log10(FreqCycles(1)),log10(FreqCycles(2)),num_frex));
            end
        end
        min_cycle_3 = find(cyclex<3);
        if ~isempty(min_cycle_3)
            cyclex(min_cycle_3) = 3;
            warning('Less than 3 cycles found for 1 or more frequencies. Corrected')
        end
    end
    
    % create legend
    legend_array = {};
    for j=1:NumConds
        temp = [];
        tempCounter = 1;
        for i=1:size(plot_conds_matrix,1)
            if plot_conds_matrix(i,j)~=0
                temp(tempCounter) = find(CondNums==plot_conds_matrix(i,j));
                tempCounter = tempCounter+1;
            end
        end
        tempname = Conditions(temp(1));
        if length(unique(temp))>1
            tempCounter = 1;
            for i=2:length(temp)%size(plot_conds_matrix,1)
                %if plot_conds_matrix(tempCounter,j)~=0
                tempname = strcat(tempname,'+',Conditions(temp(i)));
            end
        end
        legend_array(length(legend_array)+1) = tempname;
    end
    
    %% Create groups of plot types
    % Create an array of line plots, then feed iteratively to funct
    % Need to have axes titles up at top of script
    
    line_plot_conds = {}; % What line plots are there in current plots?
    line_plot_index = []; % What's the plot number out of current plots?
    line_plot_index_full = []; % What's the plot number out of all possible plots?
    for lineplots = 1:length(PLOT_Specs.types)
        temp = char(PLOT_Specs.types(lineplots));
        k = strfind(temp,'-');
        if strcmp(temp,'ERP') || ~isempty(k)
            % Add the line plot to the array
            line_plot_conds = [line_plot_conds temp];
            
            % find the index of each line plot in the PLOT_Specs.types variable
            line_plot_index = [line_plot_index lineplots];
        end
        
        for i=1:length(SPECS.Options.types)
            k = strfind(SPECS.Options.types(i),temp);
            if ~isempty(cell2mat(k))
                line_plot_index_full = [line_plot_index_full i];
            end
        end
        
    end
    
    spect_plot_conds = {}; % What spect plots are there in current plots?
    spect_plot_index = []; % What's the plot number out of current plots?
    spect_plot_index_full = []; % What's the plot number out of all possible plots?
    for spectplots = 1:length(PLOT_Specs.types)
        temp = char(PLOT_Specs.types(spectplots));
        k = strfind(temp,'-');
        if isempty(k)
            if ~strcmp(temp,'ERP') && ~strcmp(temp,'Distribution')
                % Add the spect plot to the array
                spect_plot_conds = [spect_plot_conds temp];
                
                % find the index of each spect plot in the PLOT_Specs.types variable
                spect_plot_index = [spect_plot_index spectplots];
            end
        end
        for i=1:length(SPECS.Options.types)
            %k = strfind(SPECS.Options.types(i),temp);
            %if ~isempty(cell2mat(k))
            if strcmp(SPECS.Options.types(i),temp)
                spect_plot_index_full = [spect_plot_index_full i];
            end
        end
        
    end
            
    %% Check if we need to combine data across channels
    % If we do:
    % Update bad_trial array
    % Update channel labels
    % Change data in evoked_sig and ersp_sig variables
    
    chan_num_counter = length(chan_names)+1;
    plot_chans_virtual = [];
    
    if size(plot_chans,1)>1
        for i=1:size(plot_chans,2)
            
            % ignore 0s
            merge_chans = plot_chans(:,i)';
            merge_chans(find(merge_chans==0))=[];
            
            if length(merge_chans)>1
                
                % renumber chan numbers
                chan_names{chan_num_counter} = OUTPUT.chan_names{plot_chans(1,i)};
                for j=2:size(plot_chans,1)
                    chan_names{chan_num_counter} = [chan_names{chan_num_counter},'+',OUTPUT.chan_names{plot_chans(j,i)}];
                end
                
                % Identify any bad trial in a chan as bad trial in all chan
                Evoked_chan_bad_array(chan_num_counter,:) = max(OUTPUT.Evoked_chan_bad_array(merge_chans,:));
                
                % Update data variables
                evoked_sig(chan_num_counter,:,:) = nanmean(OUTPUT.evoked_sig(merge_chans,:,:),1);
                ersp_sig(chan_num_counter,:,:) = nanmean(OUTPUT.ersp_sig(merge_chans,:,:),1);
                oneoverf_array(chan_num_counter,:) = nanmean(oneoverf_array(merge_chans,:),1);
                
                % need to update plot_chans
                plot_chans_virtual = [plot_chans_virtual chan_num_counter];
                
                % update counter
                chan_num_counter = chan_num_counter+1;
            else
                plot_chans_virtual = [plot_chans_virtual merge_chans];
            end
            
        end
    else
        % need to update plot_chans
        plot_chans_virtual = plot_chans;
    end
        
    
    
    %% Orthogonal Pre-Selection
    % If active, only analyze electrodes that have significant activity regardless of condition
    
    if isfield(SPECS,'OrthogonalPreSelection')==0
        OrthogonalPreSelection=0;
    else
        OrthogonalPreSelection = SPECS.OrthogonalPreSelection;
    end
    if ~ismember(OrthogonalPreSelection,[0 1])
        error('SPECS.OrthogonalPreSelection should be 0 (No) or 1 (Yes)');
    end
    
    if isfield(SPECS,'OrthogonalPreSelection_PThresh')==0
        OrthogonalPreSelection_PThresh=0.05;
    else
        OrthogonalPreSelection_PThresh = SPECS.OrthogonalPreSelection_PThresh;
    end
    if OrthogonalPreSelection_PThresh<0 || OrthogonalPreSelection_PThresh>1
        error('SPECS.OrthogonalPreSelection_PThresh should be between 0 and 1');
    end    
    
    if isfield(SPECS,'OrthogonalPreSelection_Polarity')==0
        OrthogonalPreSelection_Polarity=0;
    else
        OrthogonalPreSelection_Polarity = SPECS.OrthogonalPreSelection_Polarity;
    end
    if ~ismember(OrthogonalPreSelection_Polarity,[-1 0 1])
        error('SPECS.OrthogonalPreSelection_Polarity should be -1 (Neg Sig), 0 (Both), or 1 (Pos Sig)');
    end
    
    for OrthogRun = 1:OrthogonalPreSelection+1
        clear AVG PLOT

        % Check that only running one type of analysis (e.g., ERP or ERSP)
        if OrthogonalPreSelection==1
           if size(SPECS.types,2)>1
               error('Only one analysis can be run at a time with OrthogonalPreSelection. Update SPECS.types.')
           end
        end
        
        if OrthogonalPreSelection==1 && OrthogRun==1
            % Update conditions
            NumConds_Orig = NumConds;
            plot_conds_matrix_Orig = plot_conds_matrix;
            plot_conds_matrix = plot_conds_matrix(:);
            NumConds = 1;
            
            % Update pval
            PThresh_Orig = PThresh;
            PThresh = OrthogonalPreSelection_PThresh;
            
        elseif OrthogonalPreSelection==1 && OrthogRun==2
            % Update conditions
            NumConds = NumConds_Orig;
            plot_conds_matrix = plot_conds_matrix_Orig;
            
            % Update pval
            PThresh = PThresh_Orig;
            
        end
        
        % Check if we need to combine conditions
        % Updated to ignore zeros in the columns
        stim_new = stim;
        if size(plot_conds_matrix,1)>1
            for i=1:size(plot_conds_matrix,2)
                % renumber additional stim numbers
                for j=2:size(plot_conds_matrix,1)
                    if plot_conds_matrix(j,i)~=0
                        temp = find(stim(:,2)==plot_conds_matrix(j,i));
                        stim_new(temp,2) = plot_conds_matrix(1,i);
                    end
                end
            end
        end
        
        % if OrthogonalPreSelection active, need more than one cond
        if OrthogonalPreSelection==1 && OrthogRun==1
            if NumConds_Orig==1
                error('Need more than one condition (excluding averaged conditions) for OrthogonalPreSelection');
            end
        end
        
        %% Iterate through all channels
        counter = 1;
        for chanx = plot_chans_virtual %[29 30 37 38 43:46 86:91]
            
            clear PLOT
            disp(['Channel ' int2str(counter) ' of ' int2str(size(plot_chans_virtual,2))])
            
            
            %         error('a');
            for j=1:NumConds
                
                % Pull  data from the current channel for ERP and ERSP
                PLOT.(genvarname(['ERP',num2str(plot_conds(j))])) = squeeze(evoked_sig(chanx,:,:));
                PLOT.(genvarname(['ERSP_ERP',num2str(plot_conds(j))])) = squeeze(ersp_sig(chanx,:,:));
                
                % Only include trials for this condition (others nan)
                exclude_trials = find(stim_new(:,2)~=plot_conds(j));
                PLOT.(genvarname(['ERP',num2str(plot_conds(j))]))(:,exclude_trials)=nan;
                PLOT.(genvarname(['ERSP_ERP',num2str(plot_conds(j))]))(:,exclude_trials)=nan;
                
                % Change bad trials for this chan into nan
                exclude_trials = find(Evoked_chan_bad_array(chanx,:)==1);
                PLOT.(genvarname(['ERP',num2str(plot_conds(j))]))(:,exclude_trials)=nan;
                PLOT.(genvarname(['ERSP_ERP',num2str(plot_conds(j))]))(:,exclude_trials)=nan;
                
                % Remove baseline from ERP
                PLOT.(genvarname(['ERP',num2str(plot_conds(j))])) = ...
                    PLOT.(genvarname(['ERP',num2str(plot_conds(j))])) - ...
                    nanmean(PLOT.(genvarname(['ERP',num2str(plot_conds(j))]))(plotbaseline,:));
                
                % Remove baseline from ERSP
                % This is the initial baseline removal from the pre-filtered
                % data. It is necessary to mean center data before filtering
                % After filtering we'll apply another baseline removal
                PLOT.(genvarname(['ERSP_ERP',num2str(plot_conds(j))])) = ...
                    PLOT.(genvarname(['ERSP_ERP',num2str(plot_conds(j))])) - ...
                    nanmean(PLOT.(genvarname(['ERSP_ERP',num2str(plot_conds(j))]))(plotbaseline,:));
                
                % Remove outliers from ERP
                % Use Median Absolute Deviation
                % Leys, C.; et al. (2013). "Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median" (PDF). Journal of Experimental Social Psychology. 49 (4): 764?766. doi:10.1016/j.jesp.2013.03.013.
                % Need different measure if we want to consider the time-series
                
                
                % a = [PLOT.(genvarname(['ERP',num2str(plot_conds(j))]))(1,:)];
                if ~strcmp(RemoveOutliers,'NONE')
                    PLOT.(genvarname(['ERP',num2str(plot_conds(j))])) = ECOG_filloutliers(PLOT.(genvarname(['ERP',num2str(plot_conds(j))])),OutlierMethod,RemoveOutliers);
                end
                
                %figure;histogram(PLOT.(genvarname(['ERP',num2str(plot_conds(j))]))(1,:))
                %a = [a;PLOT.(genvarname(['ERP',num2str(plot_conds(j))]))(1,:)]';
                %error('a')
                
                % Create average ERP
                AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,:) = nanmean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)';
                
                % Pull individual ERP
                AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_indiv'))))(chanx,:,:) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))));
                
                % Check if RawSpect is one of the plots, then add variable.
                if sum(strcmp(PLOT_Specs.types,'RawSpect-line')>0) || sum(strcmp(PLOT_Specs.types,'RawSpect')>0) ||  sum(strcmp(PLOT_Specs.types,'RawSpect-individ')>0)
                    PLOT.(genvarname(char(strcat('RawSpect_ERP',num2str(plot_conds(j)))))) = PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j))))));
                    
                    if SPECS.Normalization==1
                        error('RawSpect measures cannot use SPECS.Normalization=1 becuase the raw filter data yields complex values when logged. Use SPECS.Normalization 0, 2, or 3.')
                    end
                    
                end
                
            end
            
            PHASEcond = {};
            RawSpectcond = {};
            Variancecond = {};
            if SPECS.Options.ERSP==1
                
                ERSP_and_Induced = {};
                ERSPcond = {};
                
                if sum(strcmp(PLOT_Specs.types,'ERSP-line')>0) || sum(strcmp(PLOT_Specs.types,'ERSP')>0) ||  sum(strcmp(PLOT_Specs.types,'ERSP-individ')>0) ||  sum(strcmp(PLOT_Specs.types,'Distribution')>0)
                    ERSPcond = [ERSPcond 'ERSP'];
                    ERSP_and_Induced = [ERSP_and_Induced 'ERSP'];
                end
                
                if sum(strcmp(PLOT_Specs.types,'RawSpect-line')>0) || sum(strcmp(PLOT_Specs.types,'RawSpect')>0)
                    RawSpectcond = [RawSpectcond 'RawSpect'];
                end
                
                if sum(strcmp(PLOT_Specs.types,'ITPC-line')>0) || sum(strcmp(PLOT_Specs.types,'ITPC')>0)
                    PHASEcond = [PHASEcond 'ITPC'];
                end
                if sum(strcmp(PLOT_Specs.types,'R2-line')>0) || sum(strcmp(PLOT_Specs.types,'R2')>0)
                    PHASEcond = [PHASEcond 'R2'];
                end
                if sum(strcmp(PLOT_Specs.types,'ITPA-line')>0) || sum(strcmp(PLOT_Specs.types,'ITPA')>0)
                    PHASEcond = [PHASEcond 'ITPA'];
                end
                if sum(strcmp(PLOT_Specs.types,'Variance-line')>0) || sum(strcmp(PLOT_Specs.types,'Variance')>0)
                    Variancecond = [Variancecond 'Variance'];
                end
                
                if isempty(ERSPcond) && ~isempty(PHASEcond)
                    ERSPcond = [ERSPcond 'ERSP'];
                    ERSP_and_Induced = [ERSP_and_Induced 'ERSP'];
                end
                
                if isempty(ERSPcond) && isempty(PHASEcond) && ~isempty(Variancecond)
                    ERSPcond = [ERSPcond 'ERSP'];
                    ERSP_and_Induced = [ERSP_and_Induced 'ERSP'];
                end
                
                % if isempty(ERSPcond) && ~isempty(RawSpectcond)
                if ~isempty(RawSpectcond)
                    ERSPcond = [ERSPcond 'RawSpect'];
                    ERSP_and_Induced = [ERSP_and_Induced 'RawSpect'];
                end
                
                % Check if we need to do ERSP, Induced, or both
                ERSP_and_Induced_Array = {};
                if sum([strcmp(ERSP_and_Induced,'ERSP'),strcmp(ERSP_and_Induced,'Variance')])>0
                    ERSP_and_Induced_Array = [ERSP_and_Induced_Array 'ERSP'];
                end
                if sum(strcmp(ERSP_and_Induced,'Induced'))>0
                    ERSP_and_Induced_Array = [ERSP_and_Induced_Array 'Induced'];
                end
                if sum(strcmp(ERSP_and_Induced,'RawSpect'))>0
                    ERSP_and_Induced_Array = [ERSP_and_Induced_Array 'RawSpect'];
                end
                
                for ERSPcLength=1:length(ERSP_and_Induced_Array)
                    % Power calculated on freq averages for ERSP
                    for j=1:NumConds
                        
                        % Wavelet analyses require time to be symmetrical
                        % Pad the data if not symmetrical around zero
                        
                        wave_time = time-mean(time);
                        temp = PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'_ERP',num2str(plot_conds(j)))))); % time x trials format
                        
                        % For ease we'll want to remove the nans
                        % before filtering, and them put the nans back
                        nan_hit = find(isnan(temp(1,:)));
                        nan_miss = find(~isnan(temp(1,:)));
                        temp(:,nan_hit)=[];
                        
                        EEGpnts = size(temp,1);
                        EEGtrials = size(temp,2);
                        
                        
                        % definte convolution parameters
                        n_wavelet            = length(time);
                        n_data               = EEGpnts*EEGtrials;
                        n_convolution        = n_wavelet+n_data-1;
                        n_conv_pow2          = pow2(nextpow2(n_convolution));
                        half_of_wavelet_size = (n_wavelet-1)/2;
                        
                        % Get FFT of data
                        eegfft = fft(reshape(temp,1,EEGpnts*EEGtrials),n_conv_pow2);
                        
                        % initialize
                        eegpower = zeros(num_frex,EEGpnts); % frequencies X time X trials
                        
                        % pull current 1/f
                        oneoverf = oneoverf_array(chanx,:);%OUTPUT.oneoverf(chanx,:);
                        oneoverf_mean = zeros(1,num_frex);
                        
                        % loop through frequencies and filter
                        s = cyclex./(2*pi*frex);
                        
                        % For 1/f calculation
                        hz = linspace(0,srate/2,floor(length(time)/2)+1);
                        
                        for fi=1:num_frex
                            
                            fwave = exp(2*1i*pi*frex(fi).*wave_time) .* exp(-wave_time.^2./(2*(s(fi)^2))) ;
                            wavelet = fft( fwave, n_conv_pow2 );
                            
                            % convolution
                            eegconv = ifft(wavelet.*eegfft);
                            eegconv = eegconv(1:n_convolution);
                            eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
                            sig1 = reshape(eegconv,EEGpnts,EEGtrials);
                            
                            % test norm by 1/f
                            fwave_filt = fft(fwave);
                            fwave_filt = abs(fwave_filt(1:length(hz)))*2;
                            oneoverf_mean(fi) = trapz(hz,oneoverf.*fwave_filt); % area under curve
                            
                            % Put nans back
                            temp = sig1;
                            sig1 = nan(length(time),length(nan_miss)+length(nan_hit));
                            sig1(:,nan_miss)=temp;
                            
                            % Store values
                            PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j))))))(fi,:,:) = sig1.*conj(sig1); % same as mean(abs(sig1).^2,2);
                            % PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j))))))(fi,:,:) = abs(sig1); % without transforming into power
                            
                            % Remove outliers from ERSP
                            % Can't do this here since not baselined
                            % Use Median Absolute Deviation
                            %                         PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j))))))(fi,:,:) = filloutliers(PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j))))))(fi,:,:),'clip',3);
                            %                         PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j))))))(fi,:,:) = filloutliers(PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j))))))(fi,:,:),nan,3);
                            %                         PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j))))))(fi,:,:) = filloutliers(PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j))))))(fi,:,:),'center',3);
                            
                            if length(PHASEcond)>0
                                PLOT.(genvarname(char(strcat('Phase',num2str(plot_conds(j))))))(fi,:,:) = angle(sig1); % hilbert isn't needed because filtered response is already complex
                            end
                            if length(RawSpectcond)>0
                                PLOT.(genvarname(char(strcat('RawSpect',num2str(plot_conds(j))))))(fi,:,:) = real(sig1);
                            end
                            
                        end
                    end
                    
                    for j=1:NumConds
                        
                        %Frequency Normalization
                        curr_plot_type = [];
                        if sum(strcmp(PLOT_Specs.types,strcat(ERSP_and_Induced_Array(ERSPcLength),'-line'))>0),curr_plot_type=[curr_plot_type 1];end
                        if sum(strcmp(PLOT_Specs.types,strcat(ERSP_and_Induced_Array(ERSPcLength)))>0),curr_plot_type=[curr_plot_type 2];end
                        if sum(strcmp(PLOT_Specs.types,strcat(ERSP_and_Induced_Array(ERSPcLength),'-individ'))>0),curr_plot_type=[curr_plot_type 3];end
                        [tempDatIndiv,tempDat,tempDatSpectral,tempDatSpectral_nobase]=ECOG_plotter_dep_normalize_v12(PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j)))))),zbaseline,xaxis,SPECS.Normalization,curr_plot_type,oneoverf_mean,ERSP_and_Induced_Array{ERSPcLength},OutlierMethod,RemoveOutliers);
                        
                        % save unprocessed data to AVG
                        % AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(chanx,:,:,:) = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))));
                        
                        % Save data to variables post-normalization
                        PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j)),'_spect_nobase')))) = tempDatSpectral_nobase;
                        if ismember(1,curr_plot_type),AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j)),'_avg'))))(chanx,:) = tempDat;end
                        if ismember(2,curr_plot_type)
                            PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j)),'_spect')))) = tempDatSpectral;
                            AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j)),'_spect'))))(chanx,:,:) = nanmean(tempDatSpectral,1);
                        end
                        if sum(ismember([1 3],curr_plot_type))>0
                            PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j)),'_indiv')))) = tempDatIndiv';
                            
                            % Also save the individual responses into an array in the AVG variable
                            % Have nans in empty trials to keep number consistent across electrodes
                            %                         AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j)),'_indiv'))))(chanx,:,:) = nan(size(tempDatIndiv,2),length(find(stim_new(:,2)==plot_conds(j))));
                            %                         temp_good = find(OUTPUT.Evoked_chan_bad_array(chanx,find(stim_new(:,2)==plot_conds(j)))==0);
                            %                         AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j)),'_indiv'))))(chanx,:,temp_good) = tempDatIndiv';
                            AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j)),'_indiv'))))(chanx,:,:) = tempDatIndiv';
                            
                        end
                    end
                end
            end
            
            PThresh_adj = PThresh*2/PTail;
            
            %% Statistics
            if SPECS.StatsType == 1 % Parametric
                
                % *************************************************************
                % Wrap this into a for loop to do on all conditions, like plots
                % Need to push raw and corr pvals into correct matrices. E.g.,
                % ERSP-line_p and ERSP-line_pcorr
                % *************************************************************
                
                % Spectral-line
                if SPECS.Options.ERSP==1
                    for ERSPcLength=1:length(ERSP_and_Induced_Array)
                        if sum(strcmp(PLOT_Specs.types,strcat(ERSP_and_Induced_Array(ERSPcLength),'-line'))>0)
                            clear input_data
                            for Condx = 1:NumConds
                                input_data.(genvarname(['data',num2str(Condx)])) = squeeze(PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(Condx)),'_indiv')))))';
                            end
                            STATS = ECOG_plotter_dep_parametric_statistics_v8a(input_data,PThresh_adj,SPECS.ErrBars,SPECS.MultCompType,NPerm,analysis_windows);
                            AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'pvalues'))))(chanx,:) = STATS.p;
                            AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'pvalues_corr'))))(chanx,:) = STATS.multcomp;
                            
                            % Only pull if not stats 0
                            %                         if ismember(SPECS.MultCompType,0)==0
                            AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'pvalues_corr_values'))))(chanx,:) = STATS.multcompvals;
                            %                         end
                            
                            % AVG.ERSPpvalues_corr(chanx,:) = STATS.multcomp;
                            for Condx = 1:NumConds
                                AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(Condx)),'_err'))))(chanx,:) = STATS.ci(Condx,:);
                            end
                        end
                    end
                end
                
                % ERP
                if sum(strcmp(PLOT_Specs.types,'ERP'))>0
                    clear input_data
                    for Condx = 1:NumConds
                        input_data.(genvarname(['data',num2str(Condx)])) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(Condx))))))';
                    end
                    STATS = ECOG_plotter_dep_parametric_statistics_v8a(input_data,PThresh_adj,SPECS.ErrBars,SPECS.MultCompType,NPerm,analysis_windows);
                    AVG.ERPpvalues(chanx,:) = STATS.p;
                    AVG.ERPpvalues_corr(chanx,:) = STATS.multcomp;
                    AVG.ERP_EffectSize(chanx,:) = STATS.CohenD;
                    
                    % Only pull if not stats 0
                    %                 if ismember(SPECS.MultCompType,0)==0
                    AVG.ERPpvalues_corr_vals(chanx,:) = STATS.multcompvals;
                    %                 end
                    
                    for Condx = 1:NumConds
                        AVG.(genvarname(char(['ERP',num2str(plot_conds(Condx)),'_err'])))(chanx,:) = STATS.ci(Condx,:);
                    end
                end
                
                % Spectral
                if SPECS.Options.ERSP==1
                    for ERSPcLength=1:length(ERSP_and_Induced_Array)
                        if sum(strcmp(PLOT_Specs.types,ERSP_and_Induced_Array(ERSPcLength)))>0
                            clear input_data
                            for Condx = 1:NumConds
                                tempdata = PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(Condx)),'_spect'))));
                                input_data.(genvarname(['data',num2str(Condx)])) = squeeze(tempdata);
                            end
                            STATS = ECOG_plotter_dep_parametric_statistics_v8a(input_data,PThresh_adj,SPECS.ErrBars,SPECS.MultCompType,NPerm,analysis_windows);
                            AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'_spect_pvalues'))))(chanx,:,:) = STATS.p;
                            AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'_spect_pvalues_corr'))))(chanx,:,:) = STATS.multcomp;
                            
                            % Only pull if not stats 0
                            %                         if ismember(SPECS.MultCompType,0)==0
                            AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'_spect_pvalues_corr_values'))))(chanx,:,:) = STATS.multcompvals;
                            %                         end
                            
                        end
                    end
                end
                
                % ITPC & R2
                % Stats calculated in the same function for both line and
                % spectral computations since ITPC calculated independently at
                % each frequency bin, and averaged across freqs when
                % necessary.
                if length(PHASEcond)>0
                    for PHASEtype = 1:length(PHASEcond)
                        clear input_data
                        for Condx = 1:NumConds
                            tempdata = PLOT.(genvarname(char(strcat('Phase',num2str(plot_conds(Condx))))));
                            input_data.(genvarname(['data',num2str(Condx)])) = squeeze(tempdata);
                        end
                        
                        if sum(strcmp(PLOT_Specs.types,strcat(PHASEcond(PHASEtype),'-line'))>0)
                            line_or_spect=1; % Line
                        else
                            line_or_spect=2; % Spect
                        end
                        
                        if strcmp(PHASEcond(PHASEtype),'ITPC')
                            
                            % Check how many tails being used
                            % In general, ITPC should increase relative to baseline.
                            % I'm not sure mechanistically why ITPC would decrease below baseline
                            % unless the baseline had high ITPC values
                            if PTail==2
                                warning('Parametric ITPC calculations are 1-tailed by default')
                            end
                            [STATS,Output] = ECOG_plotter_dep_parametric_statistics_ITPC_v4i(input_data,PThresh,SPECS.ErrBars,SPECS.MultCompType,zbaseline,line_or_spect,analysis_windows);
                            
                            % Save mean permuted ITPC values
                            % AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'_spect_perm'))))(chanx,:,:) = STATS.perm_diff;
                            
                        elseif strcmp(PHASEcond(PHASEtype),'ITPA')
                            [STATS,Output] = ECOG_plotter_dep_parametric_statistics_ITPA_v1(input_data,PThresh,PTail,SPECS.ErrBars,SPECS.MultCompType,NPerm,zbaseline,line_or_spect,analysis_windows);
                        end
                        
                        
                        if line_or_spect==1
                            for Condx = 1:NumConds
                                AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),num2str(plot_conds(Condx)),'_avg'))))(chanx,:) = Output.(genvarname(char(strcat('data_Cluster',num2str(Condx)))));
                                AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),num2str(plot_conds(Condx)),'_err'))))(chanx,:) = STATS.ci(Condx,:);
                                AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'pvalues_corr'))))(chanx,:,:) = STATS.multcomp;
                                
                                % Only pull if not stats 0, 2 and 3.
                                %                             if ismember(SPECS.MultCompType,[0 2 3])==0
                                AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'pvalues_corr_values'))))(chanx,:,:) = STATS.multcompvals;
                                %                             end
                            end
                        else
                            for Condx = 1:NumConds
                                AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),num2str(plot_conds(Condx)),'_spect'))))(chanx,:,:) = Output.(genvarname(char(strcat('data_Cluster',num2str(Condx)))))';
                            end
                            AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'_spect_pvalues'))))(chanx,:,:) = STATS.p;
                            AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'_spect_pvalues_corr'))))(chanx,:,:) = STATS.multcomp;
                            
                            % % % % %                         % Only pull if not stats 0, 2 and 3.
                            % % % % %                         if ismember(SPECS.MultCompType,[0 2 3])==0
                            AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'_spect_pvalues_corr_values'))))(chanx,:,:) = STATS.multcompvals;
                            % % % % %                         end
                        end
                    end
                end
                
            else
                
                % Need to implement non-parametrics
                % Initial will be just for ITPC
                
                %             error('a')
                
                % Spectral-line
                if isempty(Variancecond)==0
                    
                    % Variance measures can only be calculated on raw power
                    % values, so normalization needs to be either:
                    % 0 = No norm
                    % 2 = 1/f norm
                    % dB, Z, or percent all are calculated relative to a
                    %  pre-stimulus baseline, meaning they have negative values
                    
                    if ismember(SPECS.Normalization,[0 2])==0
                        error('Variance measures require normalization options 0 or 2.')
                    end
                    
                    clear input_data
                    for Condx = 1:NumConds
                        input_data.(genvarname(['data',num2str(Condx)])) = squeeze(PLOT.(genvarname(['ERSP' num2str(plot_conds(Condx)),'_spect_nobase'])));
                    end
                    
                    if sum(strcmp(PLOT_Specs.types,['Variance','-line'])>0)
                        line_or_spect=1; % Line
                    else
                        line_or_spect=2; % Spect
                    end
                    
                    [STATS,Output] = ECOG_plotter_dep_nonparametric_statistics_Variance_v1a(input_data,PThresh_adj,SPECS.ErrBars,SPECS.MultCompType,NPerm,zbaseline,line_or_spect,analysis_windows,xaxis);
                    
                    % pull average
                    for Condx = 1:NumConds
                        
                        if line_or_spect==1
                            AVG.(genvarname(char(strcat('Variance',num2str(plot_conds(Condx)),'_avg'))))(chanx,:) = Output.(genvarname(char(strcat('data_Cluster',num2str((Condx))))));
                            AVG.Variancepvalues(chanx,:) = STATS.p;
                            AVG.Variancepvalues_corr(chanx,:) = STATS.multcomp;
                            AVG.Variancepvalues_corr_values(chanx,:) = STATS.multcompvals;
                            
                            % CI
                            AVG.(genvarname(char(strcat('Variance',num2str(plot_conds(Condx)),'_err'))))(chanx,:) = STATS.ci(Condx,:);
                        else
                            AVG.(genvarname(char(strcat('Variance',num2str(plot_conds(Condx)),'_spect'))))(chanx,:,:) = permute(Output.(genvarname(char(strcat('data_Cluster',num2str((Condx)))))),[2 1]);
                            %                             AVG.Variance_spect_pvalues(chanx,:,:) = permute(STATS.p,[2 1]);
                            %                             AVG.Variance_spect_pvalues_corr(chanx,:,:) = permute(STATS.multcomp,[2 1]);
                            %                             AVG.Variance_spect_pvalues_corr_values(chanx,:,:) = permute(STATS.multcompvals,[2 1]);
                            AVG.Variance_spect_pvalues(chanx,:,:) = STATS.p;
                            AVG.Variance_spect_pvalues_corr(chanx,:,:) = STATS.multcomp;
                            AVG.Variance_spect_pvalues_corr_values(chanx,:,:) = STATS.multcompvals;
                            
                        end
                    end
                    
                    
                    % AVG.ERSPpvalues_corr(chanx,:) = STATS.multcomp;
                end
                
                % % % % %             % Spectral
                % % % % %             if SPECS.Options.ERSP==1
                % % % % %                 for ERSPcLength=1:length(ERSP_and_Induced_Array)
                % % % % %                     if sum(strcmp(PLOT_Specs.types,ERSP_and_Induced_Array(ERSPcLength)))>0
                % % % % %                         clear input_data
                % % % % %                         for Condx = 1:NumConds
                % % % % %                             tempdata = PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(Condx)),'_spect'))));
                % % % % %                             input_data.(genvarname(['data',num2str(Condx)])) = squeeze(tempdata);
                % % % % %                         end
                % % % % %                         STATS = ECOG_plotter_dep_parametric_statistics_v8a(input_data,PThresh_adj,SPECS.ErrBars,SPECS.MultCompType,NPerm,analysis_windows);
                % % % % %                         AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'_spect_pvalues'))))(chanx,:,:) = STATS.p;
                % % % % %                         AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'_spect_pvalues_corr'))))(chanx,:,:) = STATS.multcomp;
                % % % % %
                % % % % %                         % Only pull if not stats 0
                % % % % % %                         if ismember(SPECS.MultCompType,0)==0
                % % % % %                             AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'_spect_pvalues_corr_values'))))(chanx,:,:) = STATS.multcompvals;
                % % % % % %                         end
                % % % % %
                % % % % %                     end
                % % % % %                 end
                % % % % %             end
                
                
                
                
                % ITPC & R2
                % Stats calculated in the same function for both line and
                % spectral computations since ITPC calculated independently at
                % each frequency bin, and averaged across freqs when
                % necessary.
                if length(PHASEcond)>0
                    for PHASEtype = 1:length(PHASEcond)
                        clear input_data
                        for Condx = 1:NumConds
                            tempdata = PLOT.(genvarname(char(strcat('Phase',num2str(plot_conds(Condx))))));
                            input_data.(genvarname(['data',num2str(Condx)])) = squeeze(tempdata);
                        end
                        
                        if sum(strcmp(PLOT_Specs.types,strcat(PHASEcond(PHASEtype),'-line'))>0)
                            line_or_spect=1; % Line
                        else
                            line_or_spect=2; % Spect
                        end
                        
                        if strcmp(PHASEcond(PHASEtype),'ITPC')
                            
                            % Check how many tails being used
                            % In general, ITPC should increase relative to baseline.
                            % I'm not sure mechanistically why ITPC would decrease below baseline
                            % unless the baseline had high ITPC values
                            if PTail==2 && NumConds==1
                                warning('Using Two-Tailed Option with only one condition. Consider using One-Tailed tests, as single-condition ITPC will be expected to increase relative to baseline in most contexts.')
                            end
                            
                            % Time required for cluster stats scale non-linearly with increasing number of perms.
                            % To counter-act this, first test with a lower number of perms to see if p-val is in conf int
                            
                            if ~SPECS.MultCompType==3 || NPerm <=200
                                [STATS,Output] = ECOG_plotter_dep_nonparametric_statistics_ITPC_v1a(input_data,PThresh_adj,SPECS.ErrBars,SPECS.MultCompType,NPerm,zbaseline,line_or_spect,analysis_windows,xaxis);
                            else
                                % adjust the p-val to be at the upper limit of
                                % the confidence interval
                                % https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Randomise/Theory
                                NPerm_temp = 100;
                                PThresh_temp = PThresh+2*sqrt(PThresh*(1-PThresh)/NPerm_temp);
                                [STATS,Output] = ECOG_plotter_dep_nonparametric_statistics_ITPC_v1a(input_data,PThresh_adj,SPECS.ErrBars,SPECS.MultCompType,NPerm_temp,zbaseline,line_or_spect,analysis_windows,xaxis);
                                
                                % Check if anything sig; if so, run all perms
                                if isempty(find(STATS.multcomp==1))==0
                                    [STATS,Output] = ECOG_plotter_dep_nonparametric_statistics_ITPC_v1a(input_data,PThresh_adj,SPECS.ErrBars,SPECS.MultCompType,NPerm,zbaseline,line_or_spect,analysis_windows,xaxis);
                                end
                            end
                            
                            % Save mean permuted ITPC values
                            % AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'_spect_perm'))))(chanx,:,:) = STATS.perm_diff;
                            
                        end
                        
                        
                        if line_or_spect==1
                            for Condx = 1:NumConds
                                AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),num2str(plot_conds(Condx)),'_avg'))))(chanx,:) = Output.(genvarname(char(strcat('data_Cluster',num2str(Condx)))));
                                AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),num2str(plot_conds(Condx)),'_err'))))(chanx,:) = STATS.ci(Condx,:);
                                AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'pvalues_corr'))))(chanx,:,:) = STATS.multcomp;
                                
                                % Only pull if not stats 0, 2 and 3.
                                %                             if ismember(SPECS.MultCompType,[0 2 3])==0
                                AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'pvalues_corr_values'))))(chanx,:,:) = STATS.multcompvals;
                                %                             end
                            end
                        else
                            for Condx = 1:NumConds
                                AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),num2str(plot_conds(Condx)),'_spect'))))(chanx,:,:) = Output.(genvarname(char(strcat('data_Cluster',num2str(Condx)))))';
                            end
                            AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'_spect_pvalues'))))(chanx,:,:) = STATS.p;
                            AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'_spect_pvalues_corr'))))(chanx,:,:) = STATS.multcomp;
                            
                            % % % % %                         % Only pull if not stats 0, 2 and 3.
                            % % % % %                         if ismember(SPECS.MultCompType,[0 2 3])==0
                            AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'_spect_pvalues_corr_values'))))(chanx,:,:) = STATS.multcompvals;
                            % % % % %                         end
                        end
                    end
                end
                
                
                
                
                
            end
            counter = counter+1;
        end
        
        % error('a')
        if OrthogonalPreSelection==1 && OrthogRun==1
            
            % Find fieldnames
            Index = find(contains(fieldnames(AVG),'pvalues_corr'),1);
            temp = fieldnames(AVG);
            pvalues_name = temp{Index};
            analysis_name = [pvalues_name(1:end-12) num2str(plot_conds(1)) '_avg'];

            % Combine all analysis windows into 1
            windows_fields = fieldnames(analysis_windows);
            
            % Create single array with all window time-points for some multcomps
            analysis_windows_all = [];
            for windowx=1:length(windows_fields)
                analysis_windows_all = [ analysis_windows_all analysis_windows.(genvarname(char(['TIMES',num2str(windowx)])))];
            end

            % identify sig channels
            temp = AVG.(genvarname(pvalues_name))(:,analysis_windows_all);
            
            % Get polarity
            temp_polarity = sign(AVG.(genvarname(analysis_name))(:,analysis_windows_all));
            
            % Combine
            temp = temp.*temp_polarity;
            
            if OrthogonalPreSelection_Polarity==-1
                temp((temp==1)) = 0;
                temp = abs(temp);
            elseif OrthogonalPreSelection_Polarity==1
                temp((temp==-1)) = 0;
            else
                temp = abs(temp);
            end
            
            Sig_Chans = find(nansum(temp,2)>0);
            
            % update plot_chans_virtual
            plot_chans_virtual = Sig_Chans';
            
        end
    end
    SPECS.plot_chans_virtual = plot_chans_virtual;
end







%% Plotting
disp('Plotting Channels')
clear f1
if SPECS.FigSize ~=-1
    % ECOG_plotter_dep_showFigures
    counter=1;
    
    % Check if Replot used incorrectly
    if ~exist('plot_chans_virtual') && SPECS.Replot==1
       error('SPECS.Replot=1, but nothing to replot') 
    end
    
    for chanx = fliplr(plot_chans_virtual) %[29 30 37 38 43:46 86:91]
        
        f1(counter)=figure(chanx);
        counter = counter+1;
        
        if FullSrnFig==1
            set(0,'Units','pixels')
            screen_size = get(0,'ScreenSize');
            set(f1,'Position',[0 0 screen_size(3) screen_size(4)]);
        end
        
        % First plot the lines
        for plotX = 1:length(line_plot_conds)
            subplot(PLOT_Specs.rowcount,PLOT_Specs.colcount,PLOT_Specs.plotlocs(line_plot_index(plotX),:))
            ECOG_plotter_dep_figures_line_v3(zerotime,NumConds,AVG,chanx,xaxis,plot_conds,colorarray,LineThickness,PLOT_Specs.types(line_plot_index(plotX)),time,SPECS.Options.ylabel(line_plot_index_full(line_plot_index(plotX))),SPECS,legend_array)
        end
        
        % Then plot the spectrals
        temp = spect_plot_conds;
        temp = unique(temp);
        for i = 1:length(temp)
            plotX = 1+((i-1)*NumConds);
            plotXtemp = [];
            for Condx = 1:NumConds
                plotXtemp = [plotXtemp plotX-1+Condx];
            end
            ECOG_plotter_dep_figures_spectral_v3(NumConds,AVG,chanx,xaxis,plot_conds,PLOT_Specs.types(spect_plot_index(plotX)),SPECS.Options.ylabel(spect_plot_index_full(spect_plot_index(plotX))),legend_array,PLOT_Specs,time,frex,spect_plot_index([plotXtemp]),loglin_freq)
        end
        
        
        % Individual plots    
        
        % Distribution plots - mainly for tonotopy
        
        
        
        suptitle(strcat('Channel:',{' '}, chan_names(chanx)));
        
        
    end
end


% % % % % end

