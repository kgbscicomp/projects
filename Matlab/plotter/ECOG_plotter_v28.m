% % % % % function ECOG_plotter_v21

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
% SPECS.types = {'ERP','ERSP'};% SPECS.types = {'ERP','ERSP-line','ERSP','ERSP-individ','Induced-line','Induced','Induced-individ','ITPC-line','ITPC','R2-line','R2','Distribution'};
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



%% Updated with changes of the pre-processing scripts as needed
% matching versions of the pre-processing scripts
% if changes there don't affect changes here, can have mult options
Expected_Versions = [8 9 10]; 

if OUTPUT.reference==2 && OUTPUT.Script_Version==9
    error('There was a problem with bipolar referencing in Filter Analysis v9. Reprocess with newer script.')
end

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
    clear AVG PLOT Analysis
    plot_conds_matrix=SPECS.plot_conds';
    plot_conds = plot_conds_matrix(1,:);
    plot_chans=SPECS.plot_chans';
    SPECS.Options.types = {'ERP','ERSP-line','ERSP','ERSP-individ','Induced-line','Induced','Induced-individ','ITPC-line','ITPC','R2-line','R2','Distribution'};
    SPECS.Options.plots = [1 1 .5 1 1 .5 1 1 .5 1 .5 1]; % how many of two rows will each plot take in the subplot
    SPECS.Options.ylabel = {'microVolts','Normalized Amplitude','Normalized Amplitude','Normalized Amplitude','Normalized Amplitude','Normalized Amplitude','Normalized Amplitude','Phase-Consistency (0-1)','Phase-Consistency (0-1)','Phase-Consistency (0-1)','Phase-Consistency (0-1)','Mean Normalized Amplitude'};
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
    if length(PLOT_Specs.types)>1 || isempty(strmatch('ERP', SPECS.types))==1, SPECS.Options.ERSP = 1;else SPECS.Options.ERSP=0; end
    
    if SPECS.Options.ERSP==1
        if SPECS.Powerbaseline(1)<SPECS.PLOTaxis(1) || SPECS.Powerbaseline(2)>SPECS.PLOTaxis(2)
            disp('Warning. ERSP Baseline is outside of plotted range')
        end
    end
    
    % Enter condition numbers interested in plotting
    xaxis = dsearchn(time',SPECS.PLOTaxis');
    xaxis = xaxis(1):xaxis(2);
    zerotime=dsearchn(time',0);
    plotbaseline = dsearchn(time',SPECS.ERPbaseline');
    plotbaseline = plotbaseline(1):plotbaseline(2);
    zbaseline = dsearchn(time',SPECS.Powerbaseline');
    zbaseline = zbaseline(1):zbaseline(2);
    
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
            cyclex = [FreqCycles(1) round(FreqCycles(1).*(frex(2:end)./frex(1)))];
        elseif length(FreqCycles)==2
            if loglin_freq==0
                cyclex = round(linspace(FreqCycles(1),FreqCycles(2),num_frex));
            elseif loglin_freq==1
                cyclex  = round(logspace(log10(FreqCycles(1)),log10(FreqCycles(2)),num_frex));
            end
        end
    end
    
    % create legend
    legend_array = {};
    for j=1:NumConds
        temp = [];
        tempCounter = 1;
        for i=1:size(plot_conds_matrix,1)
            if plot_conds_matrix(i,j)~=0
                temp(tempCounter) = find(CondNums==plot_conds_matrix(tempCounter,j));
                tempCounter = tempCounter+1;
            end
        end
        tempname = Conditions(temp(1));
        if length(unique(temp))>1
            tempCounter = 1;
            for i=2:size(plot_conds_matrix,1)
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
    
    
    
    
    %% Check if we need to combine conditions
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
                evoked_sig(chan_num_counter,:,:) = mean(OUTPUT.evoked_sig(merge_chans,:,:),1);
                ersp_sig(chan_num_counter,:,:) = mean(OUTPUT.ersp_sig(merge_chans,:,:),1);
                
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
    
    %% Preallocate up to the max size of any var, then make sure to pull only
    %relevant info
    NTrials_ERP=[];
    NTrials_ERSP=[];
    for chanx = plot_chans_virtual
        for j=1:NumConds
            tempsizeERP=0;
            tempsizeERSP=0;
            for i=1:size(Evoked_chan_bad_array,2)
                if Evoked_chan_bad_array(chanx,i)==0
                    if stim_new(i,2)==plot_conds(j)
                        tempsizeERP=tempsizeERP+1;
                    end
                end
                if Evoked_chan_bad_array(chanx,i)==0
                    if stim_new(i,2)==plot_conds(j)
                        tempsizeERSP=tempsizeERSP+1;
                    end
                end
            end
            NTrials_ERP(chanx,j) = tempsizeERP;
            NTrials_ERSP(chanx,j) = tempsizeERSP;
        end
    end
    MaxChans = size(NTrials_ERP,1);
    MaxTrials_ERP = max(max(NTrials_ERP));
    MaxTrials_ERSP = max(max(NTrials_ERSP));
    
    % if sum(strcmp(PLOT_Specs.types,'ERP')>0),AVG.ERPpvalues_corr= nan(MaxChans,trial_length);AVG.ERPpvalues= nan(MaxChans,trial_length);end
    % if sum(strcmp(PLOT_Specs.types,'ERSP-line')>0),AVG.ERSPpvalues_corr= nan(MaxChans,trial_length);AVG.ERSPpvalues= nan(MaxChans,trial_length);end
    
    %% Iterate through all channels
    counter = 1;
    for chanx = plot_chans_virtual %[29 30 37 38 43:46 86:91]
        
        clear PLOT
        
        disp(['Channel ' int2str(counter) ' of ' int2str(size(plot_chans_virtual,2))])
        
        NTrials_ERP = [];
        NTrials_ERSP = [];
        for j=1:NumConds
            
            tempsizeERP=1;
            tempsizeERSP=1;
            for i=1:size(Evoked_chan_bad_array,2)
                if Evoked_chan_bad_array(chanx,i)==0
                    if stim_new(i,2)==plot_conds(j)
                        temp = (squeeze(evoked_sig(chanx,:,i)));
                        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,tempsizeERP) = temp;
                        tempsizeERP=tempsizeERP+1;
                    end
                end
                if Evoked_chan_bad_array(chanx,i)==0
                    if stim_new(i,2)==plot_conds(j)
                        temp = (squeeze(ersp_sig(chanx,:,i)));
                        PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j))))))(:,tempsizeERSP) = temp;
                        tempsizeERSP=tempsizeERSP+1;
                    end
                end
            end
            NTrials_ERP(j) = tempsizeERP-1;
            NTrials_ERSP(j) = tempsizeERSP-1;
            
            for i=1:NTrials_ERP(j)
                % Baseline data - Removes mean
                PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(plotbaseline,i),1));
            end
            for i=1:NTrials_ERSP(j)
                % Baseline data - Removes mean
                PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j))))))(plotbaseline,i),1));
            end
            %AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(chanx,:,1:NTrials_ERP(j)) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))));
            AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,:) = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)';
            
             % Check if Induced is one of the plots, then add variable.
            if sum(strcmp(PLOT_Specs.types,'Induced-line')>0) || sum(strcmp(PLOT_Specs.types,'Induced')>0) ||  sum(strcmp(PLOT_Specs.types,'Induced-individ')>0)
                for i=1:NTrials_ERSP(j)
                    PLOT.(genvarname(char(strcat('Induced_ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j))))))(:,i)'-mean(PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j)))))),2)';
                end
            end
        end
        
        PHASEcond = {};
        if SPECS.Options.ERSP==1
            
            ERSP_and_Induced = {};
            ERSPcond = {};
            
            if sum(strcmp(PLOT_Specs.types,'ERSP-line')>0) || sum(strcmp(PLOT_Specs.types,'ERSP')>0) ||  sum(strcmp(PLOT_Specs.types,'ERSP-individ')>0) ||  sum(strcmp(PLOT_Specs.types,'Distribution')>0)
                ERSPcond = [ERSPcond 'ERSP'];
                ERSP_and_Induced = [ERSP_and_Induced 'ERSP'];
            end
            
            % Check if Induced is one of the plots, then add variable.
            if sum(strcmp(PLOT_Specs.types,'Induced-line')>0) || sum(strcmp(PLOT_Specs.types,'Induced')>0) ||  sum(strcmp(PLOT_Specs.types,'Induced-individ')>0)
                ERSPcond = [ERSPcond 'Induced'];
                ERSP_and_Induced = [ERSP_and_Induced 'Induced'];
            end
            
            
            if sum(strcmp(PLOT_Specs.types,'ITPC-line')>0) || sum(strcmp(PLOT_Specs.types,'ITPC')>0)
                PHASEcond = [PHASEcond 'ITPC'];
            end
            if sum(strcmp(PLOT_Specs.types,'R2-line')>0) || sum(strcmp(PLOT_Specs.types,'R2')>0)
                PHASEcond = [PHASEcond 'R2'];
            end
            if isempty(ERSPcond) && ~isempty(PHASEcond)
                ERSPcond = [ERSPcond 'ERSP'];
                ERSP_and_Induced = [ERSP_and_Induced 'ERSP'];
            end
            
            % Check if we need to do ERSP, Induced, or both
            ERSP_and_Induced_Array = {};
            if sum(strcmp(ERSP_and_Induced,'ERSP'))>0
                ERSP_and_Induced_Array = [ERSP_and_Induced_Array 'ERSP'];
            end
            if sum(strcmp(ERSP_and_Induced,'Induced'))>0
                ERSP_and_Induced_Array = [ERSP_and_Induced_Array 'Induced'];
            end
            
            for ERSPcLength=1:length(ERSP_and_Induced_Array)
                % Power calculated on freq averages for ERSP
                for j=1:NumConds
                    
                    % this wavelet analysis requires time to be symmetrical
                    % not 100% sure why
                    % center time for now
                    
                    wave_time = time-mean(time);
                    temp = PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'_ERP',num2str(plot_conds(j)))))); % time x trials format
                    
                    EEGpnts = size(temp,1);
                    EEGtrials = size(temp,2);
                    
                    % definte convolution parameters
                    n_wavelet            = length(wave_time);
                    n_data               = EEGpnts*EEGtrials;
                    n_convolution        = n_wavelet+n_data-1;
                    half_of_wavelet_size = (n_wavelet-1)/2;
                    
                    % get FFT of data
                    eegfft = fft(reshape(temp,1,EEGpnts*EEGtrials),n_convolution);
                    
                    % initialize
                    eegpower = zeros(num_frex,EEGpnts); % frequencies X time X trials
                    
                    % loop through frequencies and compute synchronization
                    s_val = cyclex./(2*pi*frex);
                    
                    for fi=1:num_frex
                        wavelet_fft = fft( exp(2*1i*pi*frex(fi).*wave_time) .* exp(-wave_time.^2./(2*(s_val(fi)^2))) ,n_convolution);
                        
                        convolution_result_fft = ifft(wavelet_fft.*eegfft,n_convolution);
                        convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
                        sig1 = reshape(convolution_result_fft,EEGpnts,EEGtrials);
                        PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j))))))(fi,:,:) = sig1.*conj(sig1); % same as mean(abs(sig1).^2,2);
                        if length(PHASEcond)>0
                            PLOT.(genvarname(char(strcat('Phase',num2str(plot_conds(j))))))(fi,:,:) = angle(sig1); % hilbert isn't needed because filtered response is already complex
                        end
                    end
                end
                
                for j=1:NumConds
                    
                    %Frequency Normalization
                    curr_plot_type = [];
                    if sum(strcmp(PLOT_Specs.types,strcat(ERSP_and_Induced_Array(ERSPcLength),'-line'))>0),curr_plot_type=[curr_plot_type 1];end
                    if sum(strcmp(PLOT_Specs.types,strcat(ERSP_and_Induced_Array(ERSPcLength)))>0),curr_plot_type=[curr_plot_type 2];end
                    if sum(strcmp(PLOT_Specs.types,strcat(ERSP_and_Induced_Array(ERSPcLength),'-individ'))>0),curr_plot_type=[curr_plot_type 3];end
                    [tempDatIndiv,tempDat,tempDatSpectral]=ECOG_plotter_dep_normalize_v5(PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j)))))),zbaseline,xaxis,SPECS.Normalization,curr_plot_type);
                    
                    % save unprocessed data to AVG
                    % AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(chanx,:,:,:) = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))));
                    
                    % Save data to variables post-normalization
                    if ismember(1,curr_plot_type),AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j)),'_avg'))))(chanx,:) = tempDat;end
                    if ismember(2,curr_plot_type)
                        PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j)),'_spect')))) = tempDatSpectral;
                        AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j)),'_spect'))))(chanx,:,:) = mean(tempDatSpectral,1);
                    end
                    if sum(ismember([1 3],curr_plot_type))>0,PLOT.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),num2str(plot_conds(j)),'_indiv')))) = tempDatIndiv';end
                end
            end
        end
        
        %% Statistics
        if SPECS.StatsType == 1 % Parametric
            
            % Parametric only works on Whitened, No Normalization, and ERP plots
            if ismember(SPECS.Normalization, [2 3])
                error('Parametric Statistics can only be run on Whitened and No Normalization options (SPECS.Normalization = 2 or 3)')
            end
            
            PThresh_adj = PThresh*2/PTail;
            
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
                        STATS = ECOG_plotter_dep_parametric_statistics_v8(input_data,PThresh_adj,SPECS.ErrBars,SPECS.MultCompType,NPerm,analysis_windows,PeakPolarity);
                        AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'pvalues'))))(chanx,:) = STATS.p;
                        AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'pvalues_corr'))))(chanx,:) = STATS.multcomp;
                        
                        % Only pull if not stats 0, 2 and 3.
                        if ismember(SPECS.MultCompType,[0 2 3])==0
                            AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'pvalues_corr_values'))))(chanx,:) = STATS.multcompvals;
                        end
                        
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
                STATS = ECOG_plotter_dep_parametric_statistics_v8(input_data,PThresh_adj,SPECS.ErrBars,SPECS.MultCompType,NPerm,analysis_windows,PeakPolarity);
                AVG.ERPpvalues(chanx,:) = STATS.p;
                AVG.ERPpvalues_corr(chanx,:) = STATS.multcomp;
                
                % Only pull if not stats 0, 2 and 3.
                if ismember(SPECS.MultCompType,[0 2 3])==0
                    AVG.ERPpvalues_corr_vals(chanx,:) = STATS.multcompvals;
                end
                
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
                        STATS = ECOG_plotter_dep_parametric_statistics_v8(input_data,PThresh_adj,SPECS.ErrBars,SPECS.MultCompType,NPerm,analysis_windows,PeakPolarity);
                        AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'_spect_pvalues'))))(chanx,:,:) = STATS.p;
                        AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'_spect_pvalues_corr'))))(chanx,:,:) = STATS.multcomp;
                        
                        % Only pull if not stats 0, 2 and 3.
                        if ismember(SPECS.MultCompType,[0 2 3])==0
                            AVG.(genvarname(char(strcat(ERSP_and_Induced_Array(ERSPcLength),'_spect_pvalues_corr_values'))))(chanx,:,:) = STATS.multcompvals;
                        end
                        
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
                    
                    [STATS,Output] = ECOG_plotter_dep_parametric_statistics_ITPC_v1(input_data,PHASEcond(PHASEtype),PThresh,PTail,SPECS.ErrBars,SPECS.MultCompType,xaxis,NPerm,zbaseline,line_or_spect);
                    
                    if line_or_spect==1
                        for Condx = 1:NumConds
                            AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),num2str(plot_conds(Condx)),'_avg'))))(chanx,:) = Output.(genvarname(char(strcat('data_Cluster',num2str(Condx)))));
                            AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),num2str(plot_conds(Condx)),'_err'))))(chanx,:) = STATS.ci(Condx,:);
                            AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'pvalues_corr'))))(chanx,:,:) = STATS.multcomp;
                            
                            % Only pull if not stats 0, 2 and 3.
                            if ismember(SPECS.MultCompType,[0 2 3])==0
                                AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'pvalues_corr_values'))))(chanx,:,:) = STATS.multcompvals;
                            end
                        end
                    else
                        for Condx = 1:NumConds
                            AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),num2str(plot_conds(Condx)),'_spect'))))(chanx,:,:) = Output.(genvarname(char(strcat('data_Cluster',num2str(Condx)))))';
                        end
                        AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'_spect_pvalues'))))(chanx,:,:) = STATS.p;
                        AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'_spect_pvalues_corr'))))(chanx,:,:) = STATS.multcomp;
                        
% % % % %                         % Only pull if not stats 0, 2 and 3.
% % % % %                         if ismember(SPECS.MultCompType,[0 2 3])==0
% % % % %                             AVG.(genvarname(char(strcat(PHASEcond(PHASEtype),'_spect_pvalues_corr_values'))))(chanx,:,:) = STATS.multcompvals;
% % % % %                         end
                    end
                end
            end
            
        end
        counter = counter+1;
    end
end










%% Plotting
disp('Plotting Channels')
clear f1
if SPECS.FigSize ~=-1
    % ECOG_plotter_dep_showFigures
    counter=1;
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

