% ECoG_plotter
% This function takes the output of the function ECOG_filter_analysis7
% along with the following inputs to plot the processed ECoG data:
%
% SPECS.PLOTaxis = [-1 1];
% SPECS.plot_conds = [1 2;3 3]; % Plot 1-4 conditions, with ordered colors of red, green, blue, orange
% % if more than one column, columns concatenated, e.g., [100 101;500 500] combines stim from 100 and 101
% RejThreshEvoked = 3;
% ECOG_epoch_cleaner3
% SPECS.bad_trials = []; % Bad trials by trial number (e.g. for blinks)
% SPECS.ERPbaseline = [-.5 0]; %[-.05 0]
% SPECS.Powerbaseline = [-1 -.5];
% SPECS.xticksize = .25;
% SPECS.Freq = [3 40];%[5 40];
% SPECS.FreqWidth = [2];%2.5;
% SPECS.FreqCycles = [3 10];%[3];
% SPECS.db_on = 0; % Z if 0, dB if 1
% SPECS.plot_chans = [33];% 34 46:47];%1:length(chan_names);
% SPECS.NPerm=100; % Num permutations when needed
% SPECS.PThresh = .05; % PVal threshold for stat comparisons
% SPECS.PTail = 2; % 1 or 2 tailed tests
% SPECS.FigSize = 1; %1 for full screen, 0 for normal, -1 for don't display
% SPECS.types = {'ERP','ERSP'};% SPECS.types = {'ERP','ERSP-line','ERSP','ERSP-individ','Induced-line','Induced','Induced-individ','ITPC-line','ITPC','R2-line','R2','Distribution'};
% SPECS.DoStats = 1; %1 to get significant timepoints, pvalues (in the 
% SPECS.PrintFigs =1; %1 to print figures to a file (make sure the SPECS.imgdir is right!)
% SPECS.imgdir = '/Users/jacobzweig/Dropbox/Research/ECOG_Jacob_Analysis/UC1106/images';
% SPECS.imglabels = 'UC1106_TEST'; %Add anything you want to be a part of the filename that will be printed (chan name will be at end)
% SPECS.interpreter = 'latex' or 'tex'(default) or 'none'
% ECOG_plotter



% Oct 2016
% Added in feature to plot distribution of activity in a range across conditions. Useful for AudLoc_v2.
% Also increased the number of channels that could be compared. If more than 4 skipping stats




clear AVG PLOT Analysis
plot_conds_matrix=SPECS.plot_conds';
plot_conds = plot_conds_matrix(1,:);
plot_chans=SPECS.plot_chans;
SPECS.Options.types = {'ERP','ERSP-line','ERSP','ERSP-individ','Induced-line','Induced','Induced-individ','ITPC-line','ITPC','R2-line','R2','Distribution'};
SPECS.Options.plots = [1 1 .5 1 1 .5 1 1 .5 1 .5 1]; % how many of two rows will each plot take in the subplot
tempPLOTtypes = SPECS.types;
numanalyses = length(tempPLOTtypes);
LineThickness = 2;
RejThreshERSP = RejThreshEvoked;

% % % % % Parametric = SPECS.Parametric;


%Set up plot font interpreter
if ~exist('SPECS.interpreter','var')
    SPECS.interpreter = 'tex'; %default to tex
end
if strcmp(SPECS.interpreter,'latex')
    interpreter_string = '"FontUnits","points", "interpreter","latex", "FontSize", 13, "FontName", "Times","Location","NorthEast"';
elseif strcmp(SPECS.interpreter,'none')
    interpreter_string = '"FontUnits","points", "interpreter","none", "FontSize", 13, "FontName", "Times","Location","NorthEast"';
elseif strcmp(SPECS.interpreter,'tex')
    interpreter_string = '"FontUnits","points", "interpreter","tex", "FontSize", 13, "FontName", "Times","Location","NorthEast"';
end

PLOT.plotnum = 0;
PLOT.rowcount = size(plot_conds_matrix,2);
PLOT.colcount = 0;
PLOT.halfwhole = [];
for counter=1:numanalyses
    PLOT.colcount=PLOT.colcount+1;
    temp = SPECS.Options.plots(find(strcmp(tempPLOTtypes(counter),SPECS.Options.types)==1));
    if isempty(temp), error('One or more SPECS.types are invalid.');end
    if temp == .5, PLOT.plotnum = PLOT.plotnum+PLOT.rowcount;PLOT.halfwhole(counter)=.5;
    else PLOT.plotnum = PLOT.plotnum+1;PLOT.halfwhole(counter)=1;end
end
PLOT.arraycount = PLOT.colcount*PLOT.rowcount;
PLOT.plotarray = zeros(PLOT.rowcount,PLOT.colcount);
counter=1;
for j=1:PLOT.rowcount
    for i=1:PLOT.colcount
        PLOT.plotarray(j,i)=counter;
        counter=counter+1;
    end
end
PLOT.plotarray=PLOT.plotarray(:)';
PLOT.plotlocs = zeros(PLOT.plotnum,PLOT.rowcount);
PLOT.types = {};
counter1=1;
counter2=1;
for i=1:PLOT.colcount
    if PLOT.halfwhole(i)==.5
        for j=1:PLOT.rowcount
            PLOT.types(counter2)=tempPLOTtypes(i);
            PLOT.plotlocs(counter2,:) = ones(1,PLOT.rowcount)*PLOT.plotarray(counter1);
            counter1=counter1+1;
            counter2=counter2+1;
        end    
    else
        temp = counter1:1:counter1+PLOT.rowcount-1;
        PLOT.types(counter2)=tempPLOTtypes(i);
        PLOT.plotlocs(counter2,:) =PLOT.plotarray(temp);
        counter1=counter1+PLOT.rowcount;
        counter2=counter2+1;
    end
end


% YLIM = SPECS.Yaxis;
SPECS.Options.ERP = 1; % All processes need ERP Preprocessing
if length(PLOT.types)>1 || isempty(strmatch('ERP', SPECS.types))==1, SPECS.Options.ERSP = 1;else SPECS.Options.ERSP=0; end

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
xticksize=SPECS.xticksize;

PThresh = SPECS.PThresh;
PTail = SPECS.PTail;
if ismember(PTail,1:2)
%     PThresh = PThresh/2;
else
    error('SPECS.PTail requires entry of 1 or 2 tails');
end

if SPECS.FigSize==1,FullSrnFig=1;else FullSrnFig=0;end
db_on=SPECS.db_on;
NPerm=SPECS.NPerm;
NumConds = PLOT.rowcount;
xticklabel = SPECS.PLOTaxis(1):xticksize:SPECS.PLOTaxis(2);
for i=1:size(xticklabel,2)
    if floor(xticklabel(i)*10)==0
        xticklabel(i)=0;
    end
end

xtickrange = 1:srate*xticksize:trial_length;
% colorarray = {'r','g','b','c'};%old
colorarray=[1 0 0;0 1 0;0 0 1;0 1 1; 1 1 0; 1 0 1; 1 .5 0; 0 1 .5;0 .5 1;1 0 .5;.5 1 0;.5 0 1;.5 .5 .5; 0 0 0]; % red, green, blue, cyan, yellow, magenta, orange, blue-green, light blue,fuscia,lime-green,purple, grey,black
colorarrayShade = {'r','g','b','c'};


% if isfield(SPECS,'bad_trials') == 1, bad_trials = [];else bad_trials = SPECS.bad_trials; end



% define wavelet parameters
if SPECS.Options.ERSP==1
    Freq = SPECS.Freq;
    FreqWidth = SPECS.FreqWidth;
    FreqCycles = SPECS.FreqCycles;
    if length(Freq)==1
        frex = Freq(1);
    else
        frex = Freq(1):FreqWidth:Freq(2);
    end
    num_frex = length(frex);
    if length(FreqCycles)==1 % variable number of cycles to equate time length
        numcycles_arry = [];
        FWHM_arry = [];
        numcycles = FreqCycles(1);
        wavelet = exp(2*1i*pi*frex(1).*time) .* exp(-time.^2./(2*(numcycles/(2*pi*frex(1)))^2));
        fft_wav = 2*abs(fft(wavelet));
        hz_wav  = linspace(0,srate/2,round(length(wavelet)/2)+1);
        test = fft_wav(1:length(hz_wav));
        test = (test - min(test)) / ( max(test) - min(test) );
        [a maxval] =max((test));
        [a b] =min(abs(test(1:maxval-1)-.5));
        lefthz = hz_wav(b);
        [a b] =min(abs(test(maxval:end)-.5));
        righthz = hz_wav(b+maxval-1);
        refFWHM = righthz-lefthz;
        numcycles_arry(1) = numcycles;
        FWHM_arry(1) = refFWHM;
        
        refFWHM = numcycles_arry(1);
        for j=2:1:num_frex
            test_cycle_arry = zeros(1,50);
            for k=1:length(test_cycle_arry) % test cycles
                wavelet = exp(2*1i*pi*frex(j).*time) .* exp(-time.^2./(2*((numcycles+k-1)/(2*pi*frex(j)))^2));
                fft_wav = 2*abs(fft(wavelet));
                hz_wav  = linspace(0,srate/2,round(length(wavelet)/2)+1);
                test = fft_wav(1:length(hz_wav));
                test = (test - min(test)) / ( max(test) - min(test) );
                [a maxval] =max((test));
                [a b] =min(abs(test(1:maxval-1)-.5));
                lefthz = hz_wav(b);
                
                [a b] =min(abs(test(maxval:end)-.5));
                righthz = hz_wav(b+maxval-1);
                test_cycle_arry(k) = righthz-lefthz;
            end
            [a b] =min(abs(test_cycle_arry-FWHM_arry(1)));
            numcycles_arry(j) = numcycles_arry(1)+b-1;
        end
        cyclex = numcycles_arry;
    elseif length(FreqCycles)==2
        cyclex = linspace(FreqCycles(1),FreqCycles(2),num_frex);
    end
end

% create legend
legend_array = {};
for j=1:NumConds
    temp = zeros(1,size(plot_conds_matrix,1));
    for i=1:size(plot_conds_matrix,1)
        temp(i) = find(CondNums==plot_conds_matrix(i,j));
    end
    tempname = Conditions(temp(1));
    if length(unique(temp))>1
        for i=2:size(plot_conds_matrix,1)
            tempname = strcat(tempname,'+',Conditions(temp(i)));
        end
    end
    legend_array(length(legend_array)+1) = tempname;
end

% Check if we need to combine conditions
stim_new = stim;
if size(plot_conds_matrix,1)>1
    for i=1:size(plot_conds_matrix,2)
        % renumber additional stim numbers
        for j=2:size(plot_conds_matrix,1)
            temp = find(stim(:,2)==plot_conds_matrix(j,i));
            stim_new(temp,2) = plot_conds_matrix(1,i);
        end
    end
end

%Preallocate up to the max size of any var, then make sure to pull only
%relevant info
NTrials_ERP=[];
NTrials_ERSP=[];
for chanx = plot_chans
    for j=1:NumConds
        tempsizeERP=0;
        tempsizeERSP=0;
        for i=1:length(analyzed_trials);
            if Evoked_chan_bad_array(chanx,i)==0 && ~ismember(analyzed_trials(i),bad_trials)
                if stim_new(analyzed_trials(i),2)==plot_conds(j)
                    tempsizeERP=tempsizeERP+1;
                end
            end
            if Evoked_chan_bad_array(chanx,i)==0 && ~ismember(analyzed_trials(i),bad_trials)
                if stim_new(analyzed_trials(i),2)==plot_conds(j)
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

for j=1:NumConds
    AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))) = nan(MaxChans,trial_length,MaxTrials_ERP);
    AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))) = nan(MaxChans,trial_length);
    AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg')))) = nan(MaxChans,trial_length);
    if SPECS.Options.ERSP==1
        AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)))))) = nan(MaxChans,num_frex,trial_length,MaxTrials_ERSP);
        AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_err')))) = nan(MaxChans,trial_length);
        AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_avg')))) = nan(MaxChans,trial_length);
        if sum(strcmp(PLOT.types,'ERSP')>0),AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_spect')))) = nan(MaxChans,num_frex,trial_length);end
        if sum(strcmp(PLOT.types,'ERSP-individ')>0),AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_indiv')))) = nan(MaxChans,trial_length,MaxTrials_ERSP);end
        if sum(strcmp(PLOT.types,'Induced-line')>0) || sum(strcmp(PLOT.types,'Induced')>0) ||  sum(strcmp(PLOT.types,'Induced-individ')>0)
            AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)))))) = nan(MaxChans,num_frex,trial_length,MaxTrials_ERSP);
            AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_err')))) = nan(MaxChans,trial_length);
            AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_avg')))) = nan(MaxChans,trial_length);
        end
        if sum(strcmp(PLOT.types,'Induced')>0),AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_spect')))) = nan(MaxChans,num_frex,trial_length);end
        if sum(strcmp(PLOT.types,'Induced-individ')>0),AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_indiv')))) = nan(MaxChans,trial_length,MaxTrials_ERSP);end
        
        if sum(strcmp(PLOT.types,'ITPC-line')>0) || sum(strcmp(PLOT.types,'ITPC')>0)
            AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)))))) = nan(MaxChans,num_frex,trial_length,MaxTrials_ERSP);
            AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_err')))) = nan(MaxChans,trial_length);
            AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_avg')))) = nan(MaxChans,trial_length);
        end
        if sum(strcmp(PLOT.types,'ITPC')>0),AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_spect')))) = nan(MaxChans,num_frex,trial_length);end
        if sum(strcmp(PLOT.types,'R2-line')>0) || sum(strcmp(PLOT.types,'R2')>0)
            AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)))))) = nan(MaxChans,num_frex,trial_length,MaxTrials_ERSP);
            AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_err')))) = nan(MaxChans,trial_length);
            AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_avg')))) = nan(MaxChans,trial_length);
        end
        if sum(strcmp(PLOT.types,'R2')>0),AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_spect')))) = nan(MaxChans,num_frex,trial_length);end
    end
end
if sum(strcmp(PLOT.types,'ERP')>0),AVG.ERPstats= nan(MaxChans,trial_length);AVG.ERPpstats= nan(MaxChans,trial_length);end
if sum(strcmp(PLOT.types,'ERSP-line')>0),AVG.ERSPstats= nan(MaxChans,trial_length);end
if sum(strcmp(PLOT.types,'Induced-line')>0),AVG.Inducedstats= nan(MaxChans,trial_length);end
if sum(strcmp(PLOT.types,'ITPC-line')>0),AVG.ITPCstats= nan(MaxChans,trial_length);end
if sum(strcmp(PLOT.types,'R2-line')>0),AVG.R2stats= nan(MaxChans,trial_length);end

counter = 1;
for chanx = plot_chans %[29 30 37 38 43:46 86:91]
    disp(['Channel ' int2str(counter) ' of ' int2str(length(plot_chans))])
    for j=1:NumConds
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))) = [];
        if SPECS.Options.ERSP==1
            PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j)))))) = []; % Diff ERP in case diff trials excluded in ERP from longer time window
            PLOT.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)))))) = [];
            if sum(strcmp(PLOT.types,'Induced-line')>0) || sum(strcmp(PLOT.types,'Induced')>0) ||  sum(strcmp(PLOT.types,'Induced-individ')>0)
                PLOT.(genvarname(char(strcat('Induced_ERP',num2str(plot_conds(j)))))) = [];
                PLOT.(genvarname(char(strcat('Induced',num2str(plot_conds(j)))))) = [];
            end
            if sum(strcmp(PLOT.types,'ITPC-line')>0) || sum(strcmp(PLOT.types,'ITPC')>0)
                PLOT.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)))))) = [];
            end
            if sum(strcmp(PLOT.types,'R2-line')>0) || sum(strcmp(PLOT.types,'R2')>0)
                PLOT.(genvarname(char(strcat('R2',num2str(plot_conds(j)))))) = [];
            end
        end
    end
    
    NTrials_ERP = [];
    NTrials_ERSP = [];
    for j=1:NumConds
        
        tempsizeERP=1;
        tempsizeERSP=1;
        for i=1:length(analyzed_trials);
            if Evoked_chan_bad_array(chanx,i)==0 && ~ismember(analyzed_trials(i),bad_trials)
                if stim_new(analyzed_trials(i),2)==plot_conds(j)
                    temp = (squeeze(evoked_sig(chanx,:,analyzed_trials(i))));
                    PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,tempsizeERP) = temp;
                    tempsizeERP=tempsizeERP+1;
                end
            end
            if Evoked_chan_bad_array(chanx,i)==0 && ~ismember(analyzed_trials(i),bad_trials)
                if stim_new(analyzed_trials(i),2)==plot_conds(j)
                    temp = (squeeze(ersp_sig(chanx,:,analyzed_trials(i))));
                    PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j))))))(:,tempsizeERSP) = temp;
                    tempsizeERSP=tempsizeERSP+1;
                end
            end
        end
        NTrials_ERP(j) = tempsizeERP-1;
        NTrials_ERSP(j) = tempsizeERSP-1;
        
        for i=1:NTrials_ERP(j)
            PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(plotbaseline,i),1));
        end
        for i=1:NTrials_ERSP(j)
            PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j))))))(plotbaseline,i),1));
        end
        AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(chanx,:,1:NTrials_ERP(j)) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))));
        AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err'))))(chanx,:) = std(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))')/sqrt(NTrials_ERP(j));
        AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,:) = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)';
        if sum(strcmp(PLOT.types,'Induced-line')>0) || sum(strcmp(PLOT.types,'Induced')>0) ||  sum(strcmp(PLOT.types,'Induced-individ')>0)
            AVG.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j)),'_avg'))))(chanx,:) = mean(PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j)))))),2)';
            AVG.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j)),'_err'))))(chanx,:) = std(PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j))))))')/sqrt(NTrials_ERP(j));
            AVG.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j)),'_avg'))))(chanx,:) = mean(PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j)))))),2)';
        end
        
        %NOT SURE IF WE NEED ANYMORE
        %subtract ERP from individual trials
        if sum(strcmp(PLOT.types,'Induced-line')>0) || sum(strcmp(PLOT.types,'Induced')>0) ||  sum(strcmp(PLOT.types,'Induced-individ')>0)
            for i=1:NTrials_ERSP(j)
                PLOT.(genvarname(char(strcat('Induced_ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j))))))(:,i)-AVG.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j)),'_avg'))))(chanx,:)';
            end
        end
    end
    
    % conduct stats on the ERP
    if PThresh~=-1
    if NumConds==1 % 95% confidence interval
        [h,p,ci,stats] = ttest(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))))',0,PThresh);
        tempMean = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2)';
        AVG.ERPstats(chanx,:) = abs(ci(1,:)-tempMean);
        AVG.ERPpstats(chanx,:) = p;
    elseif NumConds==2 % ttest
        tempdata1 = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))))';
        tempdata2 = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(2))))))';
        diff = mean(tempdata1,1)-mean(tempdata2,1);
        [h,p,ci,stats] = ttest2(tempdata1, tempdata2,PThresh);
        AVG.ERPstats(chanx,:) = abs(ci(2,:)-diff)/2;
        AVG.ERPpstats(chanx,:) = p;
    elseif NumConds==3 % ANOVA
        tempgroup = [ones(1,NTrials_ERP(1)),ones(1,NTrials_ERP(2))+1,ones(1,NTrials_ERP(3))+2];
        for timex = 1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),1)
            tempdata1 = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))))(timex,:);
            tempdata2 = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(2))))))(timex,:);
            tempdata3 = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(3))))))(timex,:);
            tempstat=anova1([tempdata1,tempdata2,tempdata3],tempgroup,'off');
            if tempstat< PThresh
                AVG.ERPstats(chanx,timex)=1;
            else
                AVG.ERPstats(chanx,timex)=0;
            end
        end
    elseif NumConds==4 % ANOVA
        tempgroup = [ones(1,NTrials_ERP(1)),ones(1,NTrials_ERP(2))+1,ones(1,NTrials_ERP(3))+2,ones(1,NTrials_ERP(4))+3];
        for timex = 1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),1)
            tempdata1 = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))))(timex,:);
            tempdata2 = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(2))))))(timex,:);
            tempdata3 = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(3))))))(timex,:);
            tempdata4 = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(4))))))(timex,:);
            tempstat=anova1([tempdata1,tempdata2,tempdata3,tempdata4],tempgroup,'off');
            if tempstat< PThresh
                AVG.ERPstats(chanx,timex)=1;
            else
                AVG.ERPstats(chanx,timex)=0;
            end
        end
    else % more than 4 conds
        for timex = 1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),1)
        AVG.ERPstats(chanx,timex)=0;
        end
    end
    end
    
    
    if SPECS.Options.ERSP==1
        ERSPcond = {};
        if sum(strcmp(PLOT.types,'ERSP-line')>0) || sum(strcmp(PLOT.types,'ERSP')>0) ||  sum(strcmp(PLOT.types,'ERSP-individ')>0) ||  sum(strcmp(PLOT.types,'Distribution')>0)
            ERSPcond = [ERSPcond 'ERSP'];
        end
        if sum(strcmp(PLOT.types,'Induced-line')>0) || sum(strcmp(PLOT.types,'Induced')>0) ||  sum(strcmp(PLOT.types,'Induced-individ')>0)
            ERSPcond = [ERSPcond 'Induced'];
        end
        PHASEcond = {};
        if sum(strcmp(PLOT.types,'ITPC-line')>0) || sum(strcmp(PLOT.types,'ITPC')>0)
            PHASEcond = [PHASEcond 'ITPC'];
        end
        if sum(strcmp(PLOT.types,'R2-line')>0) || sum(strcmp(PLOT.types,'R2')>0)
            PHASEcond = [PHASEcond 'R2'];
        end
        if isempty(ERSPcond) && ~isempty(PHASEcond)
            ERSPcond = [ERSPcond 'ERSP'];
        end
        
        for ERSPcLength=1:length(ERSPcond)
            % Power calculated on freq averages for ERSP
            for j=1:NumConds
                
                % this wavelet analysis requires time to be symmetrical
                % not 100% sure why
                % center time for now
                
                wave_time = time-mean(time);
                
                temp = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),'_ERP',num2str(plot_conds(j)))))); % time x trials format
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
                    power1 = sig1.*conj(sig1); % same as mean(abs(sig1).^2,2);
                    
                    % baseline norm the power to enable sign switching perm analyses
                    %                 power1 = power1 - repmat(mean(power1(zbaseline,:),1),size(power1,1),1);
                    PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(fi,:,:) = power1;
                end
            end
            
            
            STDthresh = []; % in case we wanted to do a second level of trial rejs just on the induced. currently off (effectively)
            Zgroup = [];
            for fi=1:num_frex
                tempFreqData = [];
                for j=1:NumConds
                    tempFreqData = [tempFreqData squeeze(PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(fi,zbaseline(1):xaxis(end),:))];
                end
                Zgroup(:,:,fi) = zscore(tempFreqData);
            end
            Zgroup = mean(Zgroup,3);
            STDthresh = std(std(Zgroup))*RejThreshERSP;
            MEANthresh = mean(std(Zgroup));
            
            clear GoodTrials
            for j=1:NumConds
                C_bad_trials = [];
                Zgroup = [];
                for fi=1:num_frex
                    tempFreqData = squeeze(PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(fi,zbaseline(1):xaxis(end),:));
                    Zgroup(:,:,fi) = zscore(tempFreqData);
                end
                Zgroup = mean(Zgroup,3);
                Zgroup = std(Zgroup);
                for i=1:length(Zgroup)
                    if Zgroup(i)>(STDthresh+MEANthresh) || Zgroup(i)<(MEANthresh-STDthresh)
                        C_bad_trials(fi,i) = 1;
                    else
                        C_bad_trials(fi,i) = 0;
                    end
                end
                [I,J] = find(C_bad_trials);
                C_good_trials = 1:size(C_bad_trials,2);
                C_good_trials(J)=[];
                GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j)))))) = C_good_trials;
                GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j)),'total')))) = size(PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j)))))),3);
                
                tempDatIndiv = [];
                tempDat = [];
                for fi=1:num_frex
                    power1 = squeeze(PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(fi,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))));
                    meanDAT = mean(power1,2)';
                    meanC = mean(meanDAT(zbaseline));
                    meanIndiv = mean(power1(zbaseline,:),1);
                    stdx = std(meanDAT(zbaseline));
                    stdIndiv = std(power1(zbaseline,:));
                    if db_on==1
                        tempDB = [10*log10(meanDAT/meanC)]';
                        tempDat(fi,:) = tempDB-mean(tempDB(plotbaseline));
                        
                        tempIndivDB=[];
                        for i=1:size(power1,2)
                            tempIndivDB = [10*log10(power1(:,i)/meanIndiv(i))]';
                            tempDatIndiv(i,:,fi) = tempIndivDB-mean(tempIndivDB(plotbaseline));
                        end
                        
                    else
                        tempDat(fi,:) = (meanDAT-meanC)/stdx;
                        tempDat(fi,:) = tempDat(fi,:)-mean(tempDat(fi,plotbaseline));
                        
                        tempIndivDB=[];
                        for i=1:size(power1,2)
                            tempIndivDB = [(power1(:,i)-meanIndiv(i))/stdIndiv(i)]';
                            tempDatIndiv(i,:,fi) = tempIndivDB-mean(tempIndivDB(plotbaseline));
                        end
                    end
                end
                
                AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j)),'_avg'))))(chanx,:) = mean(tempDat,1);
                AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(chanx,:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))) = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j)))))));
                if sum(strcmp(PLOT.types,strcat(ERSPcond(ERSPcLength),'-individ'))>0),AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j)),'_indiv'))))(chanx,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))) = mean(tempDatIndiv,3)';end
                if sum(strcmp(PLOT.types,ERSPcond(ERSPcLength))>0),AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j)),'_spect'))))(chanx,:,:) = tempDat;end
                
                if sum(strcmp(PLOT.types,strcat(ERSPcond(ERSPcLength),'-line')))>0
                    % run perms for stdev for the ERSP-line analysis
                    if ismember(NumConds,[1 3:14]) && PThresh~=-1
                        tempPerm = [];
                        for k=1:NPerm
                            permAvg = [];
                            permX = randi(length(GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))),[1 length(GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j)))))))]);
                            permX = GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(permX);
                            tempDat = [];
                            for fi=1:num_frex % calculate power for each freq
                                % pull the rand freqs many times then repeat the avging above
                                meanDAT = mean(PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(fi,:,permX),3);
                                meanC = mean(meanDAT(zbaseline));
                                stdx = std(meanDAT(zbaseline));
                                if db_on==1
                                    tempDB = [10*log10(meanDAT/meanC)]';
                                    tempDat(fi,:) = tempDB-mean(tempDB(plotbaseline));
                                else
                                    tempDat(fi,:) = (meanDAT-meanC)/stdx;
                                    tempDat(fi,:) = tempDat(fi,:)-mean(tempDat(fi,plotbaseline));
                                end
                            end
                            tempPerm(k,:) = mean(tempDat,1);
                            %             disp(['Perm: ' int2str(k) ' of ' int2str(NPerm)])
                        end
                        AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j)),'_err'))))(chanx,:) = std(tempPerm);
                    end
                end
            end
            
            if sum(strcmp(PLOT.types,strcat(ERSPcond(ERSPcLength),'-line'))>0)
                % Run stats on Power
                if PThresh~=-1
                if NumConds==1 % 1-PVal confidence interval
                    % One-Sample Permuation using Sign-Switching (w replacement, i.e., not 50/50)
                    % The p-value for the test of the hypothesis that the mean of x differs from 0
                    % Good, P. 2000. Permutation Tests. Springer, New York.
                    % Same as onetPermutation in R
                    
                    % for now still use the 95% bootstrapped CIs above
                    % Bootstrapped gives similar values and easier to compute CIs
                    % % % % %
                            plotSize = size(PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1))))))));
                            realMean = AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1)),'_avg'))))(chanx,:);
                            Analysis.PermTest = zeros(1,plotSize(2));
                            Analysis.PermValues = zeros(NPerm,plotSize(2));
                            for k=1:NPerm
                                signx = randi(2,1,plotSize(3));
                                signx = reshape(signx,[1 1 plotSize(3)]);
                                signx = repmat(signx,plotSize(1),plotSize(2));
                                signx( signx==2 )=-1;
                                tempPerm = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1))))))).*signx;
                                tempDat = [];
                                for fi=1:num_frex % calculate power for each freq
                                    meanDAT = mean(tempPerm(fi,:,:),3);
                                    meanC = mean(meanDAT(zbaseline));
                                    stdx = std(meanDAT(zbaseline));
                                    if db_on==1
                                        tempDB = [10*log10(meanDAT/meanC)]';
                                        tempDat(fi,:) = tempDB-mean(tempDB(plotbaseline));
                                    else
                                        tempDat(fi,:) = (meanDAT-meanC)/stdx;
                                        tempDat(fi,:) = tempDat(fi,:)-mean(tempDat(fi,plotbaseline));
                                    end
                                end
                                Analysis.PermValues(k,:) = mean(tempDat,1);
                                for j = 1:plotSize(2)
                                    if realMean(j) > 0
                                        if realMean(j)-Analysis.PermValues(k,j)> 0
                                            Analysis.PermTest(j) = Analysis.PermTest(j)+1;
                                        end
                                    elseif realMean(j) < 0
                                        if realMean(j)-Analysis.PermValues(k,j)< 0
                                            Analysis.PermTest(j) = Analysis.PermTest(j)-1;
                                        end
                                    end
                                end
                            end
                            Analysis.PermTest = Analysis.PermTest/NPerm;
                            Analysis.PermTest = Analysis.PermTest';
                            AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'p'))))(chanx,:) = Analysis.PermTest;
                            for j = 1:plotSize(2)
                                if realMean(j) < 0, Analysis.PermTest(j,1) = -1-Analysis.PermTest(j,1);end
                                if realMean(j) > 0, Analysis.PermTest(j,1) = 1-Analysis.PermTest(j,1);end
                            end
                            CutoffValmax = ceil(quantile(1:NPerm,1-(PThresh/PTail)));
                            if NPerm*(PThresh/PTail) < 1
                               error('Not enough permutations to calculate this p-value');
                            end
                            tempA = sort(Analysis.PermValues,1);
                            % AVG.ERSPstats(chanx,:) = tempA(CutoffValmax',:)/2;
                            AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'stats'))))(chanx,:) = tempA(CutoffValmax',:);
                            AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'P'))))(chanx,:) = Analysis.PermTest;
                            
                elseif NumConds==2 % ttest
                    
% % % % %                     if Parametric==0
                        % randomly switch condition labels, resampling without replacement
                        % permutation analysis
                        
                        tempdata1 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1)))))));
                        tempdata2 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2)))))));
                        diff = AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1)),'_avg'))))(chanx,:)-AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2)),'_avg'))))(chanx,:);
                        plotSize = size(diff);
                        plotSize1 = size(tempdata1);
                        plotSize2 = size(tempdata2);
                        Analysis.PermTest = zeros(1,plotSize(2));
                        Analysis.PermValues = zeros(NPerm,plotSize(2));
                        
                        Analysis.test1n = 1:size(tempdata1,3);
                        Analysis.test2n = plotSize1(3)+1:plotSize1(3)+plotSize2(3);
                        Analysis.All = tempdata1;
                        Analysis.All(:,:,Analysis.test2n) = tempdata2;
                        for k=1:NPerm
                            Analysis.orderX = randperm(size(Analysis.All,3));
                            Analysis.permdiff  = mean(Analysis.All(:,:,Analysis.orderX(Analysis.test1n)),3)-mean(Analysis.All(:,:,Analysis.orderX(Analysis.test2n)),3);
                            tempDat = [];
                            for fi=1:num_frex % calculate power for each freq
                                meanDAT = Analysis.permdiff(fi,:,:);
                                meanC = mean(meanDAT(zbaseline));
                                stdx = std(meanDAT(zbaseline));
                                if db_on==1
                                    tempDB = [10*log10(meanDAT/meanC)]';
                                    tempDat(fi,:) = tempDB-mean(tempDB(plotbaseline));
                                else
                                    tempDat(fi,:) = (meanDAT-meanC)/stdx;
                                    tempDat(fi,:) = tempDat(fi,:)-mean(tempDat(fi,plotbaseline));
                                end
                            end
                            Analysis.PermValues(k,:) = mean(tempDat,1);
                            for j = 1:plotSize(2)
                                
                                if diff(j) > 0
                                    if Analysis.PermValues(k,j) > diff(j)
                                        Analysis.PermTest(j) = Analysis.PermTest(j)+1;
                                    end
                                elseif diff(j) < 0
                                    if Analysis.PermValues(k,j) < diff(j)
                                        Analysis.PermTest(j) = Analysis.PermTest(j)-1;
                                    end
                                end
                            end
                            %             if ismember(k,round(linspace(1,NPerm,20)))
                            %                 (k/NPerm)
                            %             end
                        end
                        Analysis.PermTest = (Analysis.PermTest'/NPerm)*PTail;
                        Analysis.(genvarname(char(strcat(ERSPcond(ERSPcLength),'_PermTest')))) = Analysis.PermTest;
                        Analysis.(genvarname(char(strcat(ERSPcond(ERSPcLength),'_PermValues')))) = Analysis.PermValues;
                        for j = 1:plotSize(2)
                            if diff(j) < 0, Analysis.PermTest(j,1) = -1-Analysis.PermTest(j,1);end
                            if diff(j) > 0, Analysis.PermTest(j,1) = 1-Analysis.PermTest(j,1);end
                        end
                        CutoffValmax = ceil(quantile(1:NPerm,1-(PThresh/PTail)));
                        if NPerm*(PThresh/PTail) < 1
                            error('Not enough permutations to calculate this p-value');
                        end
                        tempA = sort(Analysis.PermValues,1);
                        % AVG.ERSPstats(chanx,:) = tempA(CutoffValmax',:)/2;
                        AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'stats'))))(chanx,:) = tempA(CutoffValmax',:)/2;
                        AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'P'))))(chanx,:) = Analysis.PermTest;
% % % % %                     else
% % % % %                         
% % % % %                         tempdata1 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1)))))));
% % % % %                         tempdata2 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2)))))));
% % % % %                         
% % % % %                         % znorm each trial
% % % % %                         meanDat1 = [];meanDat2 = [];
% % % % %                         meanC1 = [];meanC2 = [];
% % % % %                         stdx1 = [];stdx2 = [];
% % % % %                         tempDat1 = [];tempDat2 = [];
% % % % %                         for fi=1:num_frex % calculate power for each freq
% % % % %                             meanDat1 = squeeze(tempdata1(fi,:,:));
% % % % %                             meanC1(fi,:) = squeeze(mean(tempdata1(fi,zbaseline,:),2))';
% % % % %                             stdx1(fi,:) = std(squeeze(tempdata1(fi,zbaseline,:)));
% % % % %                             meanDat2 = squeeze(tempdata2(fi,:,:));
% % % % %                             meanC2(fi,:) = squeeze(mean(tempdata2(fi,zbaseline,:),2))';
% % % % %                             stdx2(fi,:) = std(squeeze(tempdata2(fi,zbaseline,:)));
% % % % %                             tempDB1 = [];
% % % % %                             tempDB2 = [];
% % % % %                             for trialX = 1:size(meanDat1,2)
% % % % %                                 if db_on==1
% % % % %                                     tempDB1 = [10*log10(meanDat1(:,trialX)/meanC1(trialX))]';
% % % % %                                     tempDat1(fi,trialX,:) = tempDB1-mean(tempDB1(plotbaseline));
% % % % %                                 else
% % % % %                                     tempDat1(fi,trialX,:) = (meanDat1(:,trialX)-meanC1(trialX))/stdx1(fi,trialX);
% % % % %                                     tempDat1(fi,trialX,:) = tempDat1(fi,trialX,:) - mean(tempDat1(fi,trialX,plotbaseline));
% % % % %                                 end
% % % % %                             end
% % % % %                             for trialX = 1:size(meanDat2,2)
% % % % %                                 if db_on==1
% % % % %                                     tempDB2 = [10*log10(meanDat2(:,trialX)/meanC2(trialX))]';
% % % % %                                     tempDat2(fi,trialX,:) = tempDB2-mean(tempDB2(plotbaseline));
% % % % %                                 else
% % % % %                                     tempDat2(fi,trialX,:) = (meanDat2(:,trialX)-meanC2(trialX))/stdx2(fi,trialX);
% % % % %                                     tempDat2(fi,trialX,:) = tempDat2(fi,trialX,:) - mean(tempDat2(fi,trialX,plotbaseline));
% % % % %                                 end
% % % % %                             end
% % % % %                         end
% % % % %                         tempDat1 = squeeze(mean(tempDat1,1));
% % % % %                         tempDat2 = squeeze(mean(tempDat2,1));
% % % % %                         [h,p,ci,stats] = ttest2(tempDat1, tempDat2,PThresh);
% % % % %                         
% % % % %                     end
                elseif NumConds==3 % ANOVA
                    
                    % Calculate the sum of squared deviations of group means from the grand mean
                    % https://www.uvm.edu/~dhowell/StatPages/Resampling/Bootst1way/bootstrapping_oneway.html
                    % 1. find real: sum((mean of all Z - mean of each Z)^2)
                    
                    
                    tempdata1 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1)))))));
                    tempdata2 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2)))))));
                    tempdata3 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(3))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(3)))))));
                    tempMean1 = AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1)),'_avg'))))(chanx,:);
                    tempMean2 = AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2)),'_avg'))))(chanx,:);
                    tempMean3 = AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(3)),'_avg'))))(chanx,:);
                    tempMu = mean([tempMean1;tempMean2;tempMean3],1);
                    Analysis.SS = ((tempMean1-tempMu).^2)+((tempMean2-tempMu).^2)+((tempMean3-tempMu).^2);
                    plotSize = size(Analysis.SS);
                    plotSize1 = size(tempdata1);
                    plotSize2 = size(tempdata2);
                    plotSize3 = size(tempdata3);
                    Analysis.PermTest = zeros(1,plotSize(2));
                    Analysis.PermValues = zeros(NPerm,plotSize(2));
                    
                    Analysis.test1n = 1:size(tempdata1,3);
                    Analysis.test2n = plotSize1(3)+1:plotSize1(3)+plotSize2(3);
                    Analysis.test3n = plotSize1(3)+plotSize2(3)+1:plotSize1(3)+plotSize2(3)+plotSize3(3);
                    Analysis.All = tempdata1;
                    Analysis.All(:,:,Analysis.test2n) = tempdata2;
                    Analysis.All(:,:,Analysis.test3n) = tempdata3;
                    
                    for k=1:NPerm
                        Analysis.orderX = randperm(size(Analysis.All,3));
                        Analysis.Perm1n = Analysis.All(:,:,Analysis.orderX(Analysis.test1n));
                        Analysis.Perm2n = Analysis.All(:,:,Analysis.orderX(Analysis.test2n));
                        Analysis.Perm3n = Analysis.All(:,:,Analysis.orderX(Analysis.test3n));
                        
                        tempDat1 = [];
                        tempDat2 = [];
                        tempDat3 = [];
                        for fi=1:num_frex % calculate power for each freq
                            meanDAT1 = mean(Analysis.Perm1n(fi,:,:),3);
                            meanDAT2 = mean(Analysis.Perm2n(fi,:,:),3);
                            meanDAT3 = mean(Analysis.Perm3n(fi,:,:),3);
                            meanC1 = mean(meanDAT1(zbaseline));
                            meanC2 = mean(meanDAT2(zbaseline));
                            meanC3 = mean(meanDAT3(zbaseline));
                            stdx1 = std(meanDAT1(zbaseline));
                            stdx2 = std(meanDAT2(zbaseline));
                            stdx3 = std(meanDAT3(zbaseline));
                            if db_on==1
                                tempDB1 = [10*log10(meanDAT1/meanC1)]';
                                tempDat1(fi,:) = tempDB1-mean(tempDB1(plotbaseline));
                                tempDB2 = [10*log10(meanDAT2/meanC2)]';
                                tempDat2(fi,:) = tempDB2-mean(tempDB2(plotbaseline));
                                tempDB3 = [10*log10(meanDAT3/meanC3)]';
                                tempDat3(fi,:) = tempDB3-mean(tempDB3(plotbaseline));
                            else
                                tempDat1(fi,:) = (meanDAT1-meanC1)/stdx1;
                                tempDat1(fi,:) = tempDat1(fi,:)-mean(tempDat1(fi,plotbaseline));
                                tempDat2(fi,:) = (meanDAT2-meanC2)/stdx2;
                                tempDat2(fi,:) = tempDat2(fi,:)-mean(tempDat2(fi,plotbaseline));
                                tempDat3(fi,:) = (meanDAT3-meanC3)/stdx3;
                                tempDat3(fi,:) = tempDat3(fi,:)-mean(tempDat3(fi,plotbaseline));
                            end
                        end
                        tempDat1 = mean(tempDat1,1);
                        tempDat2 = mean(tempDat2,1);
                        tempDat3 = mean(tempDat3,1);
                        permMu = mean([tempDat1;tempDat2;tempDat3],1);
                        Analysis.PermValues(k,:) = ((tempDat1-permMu).^2)+((tempDat2-permMu).^2)+((tempDat3-permMu).^2);
                        
                        for j = 1:plotSize(2)
                            
                            if Analysis.SS(j) > 0
                                if Analysis.PermValues(k,j) > Analysis.SS(j)
                                    Analysis.PermTest(j) = Analysis.PermTest(j)+1;
                                end
                            elseif Analysis.SS(j) < 0
                                if Analysis.PermValues(k,j) < Analysis.SS(j)
                                    Analysis.PermTest(j) = Analysis.PermTest(j)+1;
                                end
                            end
                        end
                        %             if ismember(k,round(linspace(1,NPerm,20)))
                        %                 (k/NPerm)
                        %             end
                    end
                    Analysis.PermTest = (Analysis.PermTest'/NPerm);
                    Analysis.(genvarname(char(strcat(ERSPcond(ERSPcLength),'_PermTest')))) = Analysis.PermTest;
                    Analysis.(genvarname(char(strcat(ERSPcond(ERSPcLength),'_PermValues')))) = Analysis.PermValues;
                    for timex=1:plotSize(2)
                        if Analysis.PermTest(timex)< PThresh
                            AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'stats'))))(chanx,timex) = 1;
                        else
                            AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'stats'))))(chanx,timex) = 0;
                        end
                    end
                    AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'P'))))(chanx,:) = Analysis.PermTest;
                    
                elseif NumConds==4 % ANOVA
                    
                    
                    tempdata1 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1)))))));
                    tempdata2 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2)))))));
                    tempdata3 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(3))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(3)))))));
                    tempdata4 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(4))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(4)))))));
                    tempMean1 = AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1)),'_avg'))))(chanx,:);
                    tempMean2 = AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2)),'_avg'))))(chanx,:);
                    tempMean3 = AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(3)),'_avg'))))(chanx,:);
                    tempMean4 = AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(4)),'_avg'))))(chanx,:);
                    tempMu = mean([tempMean1;tempMean2;tempMean3;tempMean4],1);
                    Analysis.SS = ((tempMean1-tempMu).^2)+((tempMean2-tempMu).^2)+((tempMean3-tempMu).^2)+((tempMean4-tempMu).^2);
                    plotSize = size(Analysis.SS);
                    plotSize1 = size(tempdata1);
                    plotSize2 = size(tempdata2);
                    plotSize3 = size(tempdata3);
                    plotSize4 = size(tempdata4);
                    Analysis.PermTest = zeros(1,plotSize(2));
                    Analysis.PermValues = zeros(NPerm,plotSize(2));
                    
                    Analysis.test1n = 1:size(tempdata1,3);
                    Analysis.test2n = plotSize1(3)+1:plotSize1(3)+plotSize2(3);
                    Analysis.test3n = plotSize1(3)+plotSize2(3)+1:plotSize1(3)+plotSize2(3)+plotSize3(3);
                    Analysis.test4n = plotSize1(3)+plotSize2(3)+plotSize3(3)+1:plotSize1(3)+plotSize2(3)+plotSize3(3)+plotSize4(3);
                    Analysis.All = tempdata1;
                    Analysis.All(:,:,Analysis.test2n) = tempdata2;
                    Analysis.All(:,:,Analysis.test3n) = tempdata3;
                    Analysis.All(:,:,Analysis.test4n) = tempdata4;
                    for k=1:NPerm
                        Analysis.orderX = randperm(size(Analysis.All,3));
                        Analysis.Perm1n = Analysis.All(:,:,Analysis.orderX(Analysis.test1n));
                        Analysis.Perm2n = Analysis.All(:,:,Analysis.orderX(Analysis.test2n));
                        Analysis.Perm3n = Analysis.All(:,:,Analysis.orderX(Analysis.test3n));
                        Analysis.Perm4n = Analysis.All(:,:,Analysis.orderX(Analysis.test4n));
                        
                        tempDat1 = [];
                        tempDat2 = [];
                        tempDat3 = [];
                        tempDat4 = [];
                        for fi=1:num_frex % calculate power for each freq
                            meanDAT1 = mean(Analysis.Perm1n(fi,:,:),3);
                            meanDAT2 = mean(Analysis.Perm2n(fi,:,:),3);
                            meanDAT3 = mean(Analysis.Perm3n(fi,:,:),3);
                            meanDAT4 = mean(Analysis.Perm4n(fi,:,:),3);
                            meanC1 = mean(meanDAT1(zbaseline));
                            meanC2 = mean(meanDAT2(zbaseline));
                            meanC3 = mean(meanDAT3(zbaseline));
                            meanC4 = mean(meanDAT4(zbaseline));
                            stdx1 = std(meanDAT1(zbaseline));
                            stdx2 = std(meanDAT2(zbaseline));
                            stdx3 = std(meanDAT3(zbaseline));
                            stdx4 = std(meanDAT4(zbaseline));
                            if db_on==1
                                tempDB1 = [10*log10(meanDAT1/meanC1)]';
                                tempDat1(fi,:) = tempDB1-mean(tempDB1(plotbaseline));
                                tempDB2 = [10*log10(meanDAT2/meanC2)]';
                                tempDat2(fi,:) = tempDB2-mean(tempDB2(plotbaseline));
                                tempDB3 = [10*log10(meanDAT3/meanC3)]';
                                tempDat3(fi,:) = tempDB3-mean(tempDB3(plotbaseline));
                                tempDB4 = [10*log10(meanDAT4/meanC4)]';
                                tempDat4(fi,:) = tempDB4-mean(tempDB4(plotbaseline));
                            else
                                tempDat1(fi,:) = (meanDAT1-meanC1)/stdx1;
                                tempDat1(fi,:) = tempDat1(fi,:)-mean(tempDat1(fi,plotbaseline));
                                tempDat2(fi,:) = (meanDAT2-meanC2)/stdx2;
                                tempDat2(fi,:) = tempDat2(fi,:)-mean(tempDat2(fi,plotbaseline));
                                tempDat3(fi,:) = (meanDAT3-meanC3)/stdx3;
                                tempDat3(fi,:) = tempDat3(fi,:)-mean(tempDat3(fi,plotbaseline));
                                tempDat4(fi,:) = (meanDAT4-meanC4)/stdx4;
                                tempDat4(fi,:) = tempDat4(fi,:)-mean(tempDat4(fi,plotbaseline));
                            end
                        end
                        tempDat1 = mean(tempDat1,1);
                        tempDat2 = mean(tempDat2,1);
                        tempDat3 = mean(tempDat3,1);
                        tempDat4 = mean(tempDat4,1);
                        permMu = mean([tempDat1;tempDat2;tempDat3;tempDat4],1);
                        Analysis.PermValues(k,:) = ((tempDat1-permMu).^2)+((tempDat2-permMu).^2)+((tempDat3-permMu).^2)+((tempDat4-permMu).^2);
                        
                        for j = 1:plotSize(2)
                            
                            if Analysis.SS(j) > 0
                                if Analysis.PermValues(k,j) > Analysis.SS(j)
                                    Analysis.PermTest(j) = Analysis.PermTest(j)+1;
                                end
                            elseif Analysis.SS(j) < 0
                                if Analysis.PermValues(k,j) < Analysis.SS(j)
                                    Analysis.PermTest(j) = Analysis.PermTest(j)+1;
                                end
                            end
                        end
                        %             if ismember(k,round(linspace(1,NPerm,20)))
                        %                 (k/NPerm)
                        %             end
                    end
                    Analysis.PermTest = (Analysis.PermTest'/NPerm);
                    Analysis.(genvarname(char(strcat(ERSPcond(ERSPcLength),'_PermTest')))) = Analysis.PermTest;
                    Analysis.(genvarname(char(strcat(ERSPcond(ERSPcLength),'_PermValues')))) = Analysis.PermValues;
                    for timex=1:plotSize(2)
                        if Analysis.PermTest(timex)< PThresh
                            AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'stats'))))(chanx,timex) =1;
                        else
                            AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'stats'))))(chanx,timex) =0;
                        end
                    end
                    AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'P'))))(chanx,:) = Analysis.PermTest;
                    
                end
            end
        end
        end
        
        
        
        
        for PHASEcLength=1:length(PHASEcond)
            % Power calculated on freq averages for ERSP
            for j=1:NumConds
                temp = PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j)))))); % time x trials format
                EEGpnts = size(temp,1);
                EEGtrials = size(temp,2);
                
                % definte convolution parameters
                n_wavelet            = length(time);
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
                    wavelet_fft = fft( exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s_val(fi)^2))) ,n_convolution);
                    
                    convolution_result_fft = ifft(wavelet_fft.*eegfft,n_convolution);
                    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
                    sig1 = reshape(convolution_result_fft,EEGpnts,EEGtrials);
                    % power1 = sig1.*conj(sig1); % same as mean(abs(sig1).^2,2);
                    angle1 = angle(sig1);
                    
                    % baseline norm the power to enable sign switching perm analyses
                    % inactive because caused probs with dB analysis not
                    % being able to take negative values
                    % power1 = power1 - repmat(mean(power1(zbaseline,:),1),size(power1,1),1);
                    PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(j))))))(fi,:,:) = angle1;
                end
            end
            
            
            % Good and bad trials should be based on variability in power,
            % not phase. Even with high ITC=.7, stdev was stil large
            
            
%             STDthresh = []; % in case we wanted to do a second level of trial rejs just on the induced. currently off (effectively)
%             for fi=1:num_frex
%                 tempFreqData = [];
%                 for j=1:NumConds
%                     tempFreqData = [tempFreqData squeeze(PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(j))))))(fi,:,:))];
%                 end
%                 tempdat = abs(mean(exp(1i*tempFreqData),2))';
%                 STDthresh(fi) = circ_std(circ_std(tempFreqData(zbaseline(1):xaxis(end),:)'))*RejThreshERSP;
%             end
            
            for j=1:NumConds
%                 C_bad_trials = [];
%                 for fi=1:num_frex
%                     power1 = squeeze(PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(j))))))(fi,:,:));
%                     stdDAT = std(power1(zbaseline(1):xaxis(end),:));
%                     for i=1:length(stdDAT)
%                         if stdDAT(i)>STDthresh(fi)
%                             C_bad_trials(fi,i) = 1;
%                         else
%                             C_bad_trials(fi,i) = 0;
%                         end
%                     end
%                 end
%                 [I,J] = find(C_bad_trials);
%                 C_good_trials = 1:size(C_bad_trials,2);
%                 C_good_trials(J)=[];
                
                tempITPC = [];
                for fi=1:num_frex
                    angle1 = squeeze(PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(j))))))(fi,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))));
                    
                    if sum(strcmp(PHASEcond(PHASEcLength),'ITPC')>0)
                        tempdat = abs(mean(exp(1i*angle1),2))';
                    elseif sum(strcmp(PHASEcond(PHASEcLength),'R2')>0)
                        tempdat = abs(mean(exp(1i*(2*angle1)),2))';
                    end
            
                    % mean and norm of ITPC
                    tempITPC(fi,:) = tempdat-mean(tempdat(plotbaseline));
                end
                AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(j)),'_avg'))))(chanx,:) = mean(tempITPC,1); %avg ITPC across freqs
                if sum(strcmp(PLOT.types,PHASEcond(PHASEcLength))>0),AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(j)),'_spect'))))(chanx,:,:) = tempITPC;end
                AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(j))))))(chanx,:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))) = PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(j))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j)))))));
                
                if sum(strcmp(PLOT.types,strcat(PHASEcond(PHASEcLength),'-line')))>0
                    % run perms for stdev for the ERSP-line analysis
                    if ismember(NumConds,[1 3:4]) && PThresh~=-1
                        tempPerm = [];
                        for k=1:NPerm
                            permAvg = [];
                            permX = randi(length(GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))),[1 length(GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j)))))))]);
                            permX = GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(permX);
                            tempITPC = [];
                            for fi=1:num_frex % calculate power for each freq
                                % pull the rand freqs many times then repeat the avging above
                                angle1 = squeeze(PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(j))))))(fi,:,permX));
                                if sum(strcmp(PHASEcond(PHASEcLength),'ITPC')>0)
                                    tempdat = abs(mean(exp(1i*angle1),2))';
                                elseif sum(strcmp(PHASEcond(PHASEcLength),'R2')>0)
                                    tempdat = abs(mean(exp(1i*(2*angle1)),2))';
                                end
                                % mean and norm of ITPC
                                tempITPC(fi,:) = tempdat-mean(tempdat(plotbaseline));
                            end
                            tempPerm(k,:) = mean(tempITPC,1);
                            %             disp(['Perm: ' int2str(k) ' of ' int2str(NPerm)])
                        end
                        AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(j)),'_err'))))(chanx,:) = std(tempPerm);
                    end
                end
            end
            
            if sum(strcmp(PLOT.types,strcat(PHASEcond(PHASEcLength),'-line'))>0)
                % Run stats on Power
                if PThresh~=-1
                if NumConds==1 % 95% confidence interval
                    % One-Sample Permuation using Sign-Switching
                    % for now still use the 95% bootstrapped CIs above
                    % Bootstrapped gives similar values and easier to compute CIs
                    
                    plotSize = size(PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1))))))));
                    realMean = AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(1)),'_avg'))))(chanx,:);
                    Analysis.PermTest = zeros(1,plotSize(2));
                    Analysis.PermValues = zeros(NPerm,plotSize(2));
                    for k=1:NPerm
                        
%                         signx = randi(2,1,plotSize(3));
%                         signx = reshape(signx,[1 1 plotSize(3)]);
%                         signx = repmat(signx,plotSize(1),plotSize(2));
%                         signx( signx==2 )=-1;
%                         tempPerm = ((rand(1,plotSize(3))*2*pi)-pi).*signx;
                     
%                         tempPerm = (rand(plotSize(2),plotSize(3))*2*pi)-pi;
%                         tempPerm = reshape(tempPerm,[1 plotSize(2) plotSize(3)]);
%                         tempPerm = repmat(tempPerm,plotSize(1));

                        signx = (rand(1,plotSize(3))*2*pi)-pi;
                        signx = reshape(signx,[1 1 plotSize(3)]);
                        signx = repmat(signx,plotSize(1),plotSize(2));
                        tempPerm = wrapToPi(PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(1))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1)))))))+signx);

                        tempDat = [];
                        for fi=1:num_frex % calculate power for each freq
                            
                            angle1 = squeeze(tempPerm(fi,:,:));
                            
                            if sum(strcmp(PHASEcond(PHASEcLength),'ITPC')>0)
                                tempDat(fi,:) = abs(mean(exp(1i*angle1),2))';
                            elseif sum(strcmp(PHASEcond(PHASEcLength),'R2')>0)
                                tempDat(fi,:) = abs(mean(exp(1i*(2*angle1)),2))';
                            end
                            
                            % mean and norm of ITPC
                            tempDat(fi,:) = tempDat(fi,:)-mean(tempDat(fi,plotbaseline));
                            
                        end
                        Analysis.PermValues(k,:) = mean(tempDat,1);
                        for j = 1:plotSize(2)
                            if realMean(j) > 0
                                if realMean(j)-Analysis.PermValues(k,j)> 0
                                    Analysis.PermTest(j) = Analysis.PermTest(j)+1;
                                end
                            elseif realMean(j) < 0
                                if realMean(j)-Analysis.PermValues(k,j)< 0
                                    Analysis.PermTest(j) = Analysis.PermTest(j)-1;
                                end
                            end
                        end
                    end
                    Analysis.PermTest = (Analysis.PermTest'/NPerm)*PTail;
                    Analysis.(genvarname(char(strcat(PHASEcond(PHASEcLength),'_PermTest')))) = Analysis.PermTest;
                    Analysis.(genvarname(char(strcat(PHASEcond(PHASEcLength),'_PermValues')))) = Analysis.PermValues;
                    for j = 1:plotSize(2)
                        if realMean(j) < 0, Analysis.PermTest(j,1) = -1-Analysis.PermTest(j,1);end
                        if realMean(j) > 0, Analysis.PermTest(j,1) = 1-Analysis.PermTest(j,1);end
                    end
                    CutoffValmax = ceil(quantile(1:NPerm,1-(PThresh/PTail)));
                    if NPerm*(PThresh/PTail) < 1
                        error('Not enough permutations to calculate this p-value');
                    end
                    tempA = sort(Analysis.PermValues,1);
                    % AVG.ERSPstats(chanx,:) = tempA(CutoffValmax',:)/2;
                    AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),'stats'))))(chanx,:) = tempA(CutoffValmax',:);
                    AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),'P'))))(chanx,:) = Analysis.PermTest;
                    
                elseif NumConds==2 % ttest
                    
                    % randomly switch condition labels, resampling without replacement
                    % permutation analysis
                    
                    tempdata1 = PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(1))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1)))))));
                    tempdata2 = PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(2))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2)))))));
                    diff = AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(1)),'_avg'))))(chanx,:)-AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(2)),'_avg'))))(chanx,:);
                    plotSize = size(diff);
                    plotSize1 = size(tempdata1);
                    plotSize2 = size(tempdata2);
                    Analysis.PermTest = zeros(1,plotSize(2));
                    Analysis.PermValues = zeros(NPerm,plotSize(2));
                    
                    Analysis.test1n = 1:size(tempdata1,3);
                    Analysis.test2n = plotSize1(3)+1:plotSize1(3)+plotSize2(3);
                    Analysis.All = tempdata1;
                    Analysis.All(:,:,Analysis.test2n) = tempdata2;
                    for k=1:NPerm
                        Analysis.orderX = randperm(size(Analysis.All,3));
                        Analysis.perm1 = Analysis.All(:,:,Analysis.orderX(Analysis.test1n));
                        Analysis.perm2 = Analysis.All(:,:,Analysis.orderX(Analysis.test2n));
                        tempITPC1 = [];
                        tempITPC2 = [];
                        for fi=1:num_frex % calculate power for each freq
                            % pull the rand freqs many times then repeat the avging above
                            angle1 = squeeze(Analysis.perm1(fi,:,:));
                            angle2 = squeeze(Analysis.perm2(fi,:,:));
                            if sum(strcmp(PHASEcond(PHASEcLength),'ITPC')>0)
                                tempdat1 = abs(mean(exp(1i*angle1),2));
                                tempdat2 = abs(mean(exp(1i*angle2),2));
                            elseif sum(strcmp(PHASEcond(PHASEcLength),'R2')>0)
                                tempdat1 = abs(mean(exp(1i*(2*angle1)),2))';
                                tempdat2 = abs(mean(exp(1i*(2*angle2)),2))';
                            end
                            % mean and norm of ITPC
                            tempITPC1(fi,:) = tempdat1-mean(tempdat1(plotbaseline));
                            tempITPC2(fi,:) = tempdat2-mean(tempdat2(plotbaseline));
                        end
                        
                        Analysis.PermValues(k,:) = mean(tempITPC1,1)-mean(tempITPC2,1);
                        for j = 1:plotSize(2)
                            
                            if diff(j) > 0
                                if Analysis.PermValues(k,j) > diff(j)
                                    Analysis.PermTest(j) = Analysis.PermTest(j)+1;
                                end
                            elseif diff(j) < 0
                                if Analysis.PermValues(k,j) < diff(j)
                                    Analysis.PermTest(j) = Analysis.PermTest(j)-1;
                                end
                            end
                        end
                        %             if ismember(k,round(linspace(1,NPerm,20)))
                        %                 (k/NPerm)
                        %             end
                    end
                    Analysis.PermTest = (Analysis.PermTest'/NPerm)*PTail;
                    Analysis.(genvarname(char(strcat(PHASEcond(PHASEcLength),'_PermTest')))) = Analysis.PermTest;
                    Analysis.(genvarname(char(strcat(PHASEcond(PHASEcLength),'_PermValues')))) = Analysis.PermValues;
                    for j = 1:plotSize(2)
                        if diff(j) < 0, Analysis.PermTest(j,1) = -1-Analysis.PermTest(j,1);end
                        if diff(j) > 0, Analysis.PermTest(j,1) = 1-Analysis.PermTest(j,1);end
                    end
                    
                    CutoffValmax = ceil(quantile(1:NPerm,1-(PThresh/PTail)));
                    if NPerm*(PThresh/PTail) < 1
                        error('Not enough permutations to calculate this p-value');
                    end
                    tempA = sort(Analysis.PermValues,1);
                    % AVG.ERSPstats(chanx,:) = tempA(CutoffValmax',:)/2;
                    AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),'stats'))))(chanx,:) = tempA(CutoffValmax',:)/2;
                    AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),'P'))))(chanx,:) = Analysis.PermTest;
                elseif NumConds==3 % ANOVA
                    
                    % Calculate the sum of squared deviations of group means from the grand mean
                    % https://www.uvm.edu/~dhowell/StatPages/Resampling/Bootst1way/bootstrapping_oneway.html
                    % 1. find real: sum((mean of all Z - mean of each Z)^2)
                    
                    
                    tempdata1 = PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(1))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1)))))));
                    tempdata2 = PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(2))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2)))))));
                    tempdata3 = PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(3))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(3)))))));
                    tempMean1 = AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(1)),'_avg'))))(chanx,:);
                    tempMean2 = AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(2)),'_avg'))))(chanx,:);
                    tempMean3 = AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(3)),'_avg'))))(chanx,:);
                    tempMu = mean([tempMean1;tempMean2;tempMean3],1);
                    Analysis.SS = ((tempMean1-tempMu).^2)+((tempMean2-tempMu).^2)+((tempMean3-tempMu).^2);
                    plotSize = size(Analysis.SS);
                    plotSize1 = size(tempdata1);
                    plotSize2 = size(tempdata2);
                    plotSize3 = size(tempdata3);
                    Analysis.PermTest = zeros(1,plotSize(2));
                    Analysis.PermValues = zeros(NPerm,plotSize(2));
                    
                    Analysis.test1n = 1:size(tempdata1,3);
                    Analysis.test2n = plotSize1(3)+1:plotSize1(3)+plotSize2(3);
                    Analysis.test3n = plotSize1(3)+plotSize2(3)+1:plotSize1(3)+plotSize2(3)+plotSize3(3);
                    Analysis.All = tempdata1;
                    Analysis.All(:,:,Analysis.test2n) = tempdata2;
                    Analysis.All(:,:,Analysis.test3n) = tempdata3;
                    
                    for k=1:NPerm
                        Analysis.orderX = randperm(size(Analysis.All,3));
                        Analysis.perm1 = Analysis.All(:,:,Analysis.orderX(Analysis.test1n));
                        Analysis.perm2 = Analysis.All(:,:,Analysis.orderX(Analysis.test2n));
                        Analysis.perm3 = Analysis.All(:,:,Analysis.orderX(Analysis.test3n));
                        
                        tempITPC1 = [];
                        tempITPC2 = [];
                        tempITPC3 = [];
                        for fi=1:num_frex % calculate power for each freq
                            % pull the rand freqs many times then repeat the avging above
                            angle1 = squeeze(Analysis.perm1(fi,:,:));
                            angle2 = squeeze(Analysis.perm2(fi,:,:));
                            angle3 = squeeze(Analysis.perm3(fi,:,:));
                            if sum(strcmp(PHASEcond(PHASEcLength),'ITPC')>0)
                                tempdat1 = abs(mean(exp(1i*angle1),2));
                                tempdat2 = abs(mean(exp(1i*angle2),2));
                                tempdat3 = abs(mean(exp(1i*angle3),2));
                            elseif sum(strcmp(PHASEcond(PHASEcLength),'R2')>0)
                                tempdat1 = abs(mean(exp(1i*(2*angle1)),2))';
                                tempdat2 = abs(mean(exp(1i*(2*angle2)),2))';
                                tempdat3 = abs(mean(exp(1i*(2*angle3)),2))';
                            end
                            % mean and norm of ITPC
                            tempITPC1(fi,:) = tempdat1-mean(tempdat1(plotbaseline));
                            tempITPC2(fi,:) = tempdat2-mean(tempdat2(plotbaseline));
                            tempITPC3(fi,:) = tempdat3-mean(tempdat3(plotbaseline));
                        end
                        
                        
                        tempDat1 = mean(tempITPC1,1);
                        tempDat2 = mean(tempITPC2,1);
                        tempDat3 = mean(tempITPC3,1);
                        permMu = mean([tempDat1;tempDat2;tempDat3],1);
                        Analysis.PermValues(k,:) = ((tempDat1-permMu).^2)+((tempDat2-permMu).^2)+((tempDat3-permMu).^2);
                        
                        for j = 1:plotSize(2)
                            
                            if Analysis.SS(j) > 0
                                if Analysis.PermValues(k,j) > Analysis.SS(j)
                                    Analysis.PermTest(j) = Analysis.PermTest(j)+1;
                                end
                            elseif Analysis.SS(j) < 0
                                if Analysis.PermValues(k,j) < Analysis.SS(j)
                                    Analysis.PermTest(j) = Analysis.PermTest(j)+1;
                                end
                            end
                        end
                        %             if ismember(k,round(linspace(1,NPerm,20)))
                        %                 (k/NPerm)
                        %             end
                    end
                    Analysis.PermTest = (Analysis.PermTest'/NPerm);
                    Analysis.(genvarname(char(strcat(PHASEcond(PHASEcLength),'_PermTest')))) = Analysis.PermTest;
                    Analysis.(genvarname(char(strcat(PHASEcond(PHASEcLength),'_PermValues')))) = Analysis.PermValues;
                    for timex=1:plotSize(2)
                        if Analysis.PermTest(timex)< PThresh
                            AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),'stats'))))(chanx,timex) = 1;
                        else
                            AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),'stats'))))(chanx,timex) = 0;
                        end
                    end
                    AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),'P'))))(chanx,:) = Analysis.PermTest;
                    
                elseif NumConds==4 % ANOVA
                    
                    
                    
                    tempdata1 = PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(1))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1)))))));
                    tempdata2 = PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(2))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2)))))));
                    tempdata3 = PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(3))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(3)))))));
                    tempdata4 = PLOT.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(4))))))(:,:,GoodTrials.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(4)))))));
                    tempMean1 = AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(1)),'_avg'))))(chanx,:);
                    tempMean2 = AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(2)),'_avg'))))(chanx,:);
                    tempMean3 = AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(3)),'_avg'))))(chanx,:);
                    tempMean4 = AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),num2str(plot_conds(4)),'_avg'))))(chanx,:);
                    tempMu = mean([tempMean1;tempMean2;tempMean3,;tempMean4],1);
                    Analysis.SS = ((tempMean1-tempMu).^2)+((tempMean2-tempMu).^2)+((tempMean3-tempMu).^2)+((tempMean4-tempMu).^2);
                    plotSize = size(Analysis.SS);
                    plotSize1 = size(tempdata1);
                    plotSize2 = size(tempdata2);
                    plotSize3 = size(tempdata3);
                    plotSize4 = size(tempdata4);
                    Analysis.PermTest = zeros(1,plotSize(2));
                    Analysis.PermValues = zeros(NPerm,plotSize(2));
                    
                    Analysis.test1n = 1:size(tempdata1,3);
                    Analysis.test2n = plotSize1(3)+1:plotSize1(3)+plotSize2(3);
                    Analysis.test3n = plotSize1(3)+plotSize2(3)+1:plotSize1(3)+plotSize2(3)+plotSize3(3);
                    Analysis.test4n = plotSize1(3)+plotSize2(3)+plotSize3(3)+1:plotSize1(3)+plotSize2(3)+plotSize3(3)+plotSize4(3);
                    Analysis.All = tempdata1;
                    Analysis.All(:,:,Analysis.test2n) = tempdata2;
                    Analysis.All(:,:,Analysis.test3n) = tempdata3;
                    Analysis.All(:,:,Analysis.test4n) = tempdata4;
                    
                    for k=1:NPerm
                        Analysis.orderX = randperm(size(Analysis.All,3));
                        Analysis.perm1 = Analysis.All(:,:,Analysis.orderX(Analysis.test1n));
                        Analysis.perm2 = Analysis.All(:,:,Analysis.orderX(Analysis.test2n));
                        Analysis.perm3 = Analysis.All(:,:,Analysis.orderX(Analysis.test3n));
                        Analysis.perm4 = Analysis.All(:,:,Analysis.orderX(Analysis.test4n));
                        
                        tempITPC1 = [];
                        tempITPC2 = [];
                        tempITPC3 = [];
                        tempITPC4 = [];
                        for fi=1:num_frex % calculate power for each freq
                            % pull the rand freqs many times then repeat the avging above
                            angle1 = squeeze(Analysis.perm1(fi,:,:));
                            angle2 = squeeze(Analysis.perm2(fi,:,:));
                            angle3 = squeeze(Analysis.perm3(fi,:,:));
                            angle4 = squeeze(Analysis.perm4(fi,:,:));
                            if sum(strcmp(PHASEcond(PHASEcLength),'ITPC')>0)
                                tempdat1 = abs(mean(exp(1i*angle1),2));
                                tempdat2 = abs(mean(exp(1i*angle2),2));
                                tempdat3 = abs(mean(exp(1i*angle3),2));
                                tempdat4 = abs(mean(exp(1i*angle4),2));
                            elseif sum(strcmp(PHASEcond(PHASEcLength),'R2')>0)
                                tempdat1 = abs(mean(exp(1i*(2*angle1)),2))';
                                tempdat2 = abs(mean(exp(1i*(2*angle2)),2))';
                                tempdat3 = abs(mean(exp(1i*(2*angle3)),2))';
                                tempdat4 = abs(mean(exp(1i*(2*angle4)),2))';
                            end
                            % mean and norm of ITPC
                            tempITPC1(fi,:) = tempdat1-mean(tempdat1(plotbaseline));
                            tempITPC2(fi,:) = tempdat2-mean(tempdat2(plotbaseline));
                            tempITPC3(fi,:) = tempdat3-mean(tempdat3(plotbaseline));
                            tempITPC4(fi,:) = tempdat4-mean(tempdat4(plotbaseline));
                        end
                        
                        
                        tempDat1 = mean(tempITPC1,1);
                        tempDat2 = mean(tempITPC2,1);
                        tempDat3 = mean(tempITPC3,1);
                        tempDat4 = mean(tempITPC4,1);
                        permMu = mean([tempDat1;tempDat2;tempDat3;tempDat4],1);
                        Analysis.PermValues(k,:) = ((tempDat1-permMu).^2)+((tempDat2-permMu).^2)+((tempDat3-permMu).^2)+((tempDat4-permMu).^2);
                        
                        for j = 1:plotSize(2)
                            
                            if Analysis.SS(j) > 0
                                if Analysis.PermValues(k,j) > Analysis.SS(j)
                                    Analysis.PermTest(j) = Analysis.PermTest(j)+1;
                                end
                            elseif Analysis.SS(j) < 0
                                if Analysis.PermValues(k,j) < Analysis.SS(j)
                                    Analysis.PermTest(j) = Analysis.PermTest(j)+1;
                                end
                            end
                        end
                        %             if ismember(k,round(linspace(1,NPerm,20)))
                        %                 (k/NPerm)
                        %             end
                    end
                    Analysis.PermTest = (Analysis.PermTest'/NPerm);
                    Analysis.(genvarname(char(strcat(PHASEcond(PHASEcLength),'_PermTest')))) = Analysis.PermTest;
                    Analysis.(genvarname(char(strcat(PHASEcond(PHASEcLength),'_PermValues')))) = Analysis.PermValues;
                    for timex=1:plotSize(2)
                        if Analysis.PermTest(timex)< PThresh
                            AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),'stats'))))(chanx,timex) = 1;
                        else
                            AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),'stats'))))(chanx,timex) = 0;
                        end
                    end
                    AVG.(genvarname(char(strcat(PHASEcond(PHASEcLength),'P'))))(chanx,:) = Analysis.PermTest;
                end
                end
            end
        end
    end
    counter = counter+1;
end
% if SPECS.Save==1
%     disp('Saving Data')
%     save(strcat(subject,'_',task,'_PLOT'),'-v7.3', '-regexp', '^(?!(Evoked|evoked_sig|tempdata)$).');
% end
disp('Plotting Channels')
if SPECS.FigSize ~=-1 && PThresh~=-1
%     if SPECS.Yaxis==1
%     for chanx = 1:length(plot_chans)
%         for condx = 1:length(plot_conds)
%             if sum(strcmp(PLOT.types,'ERSP-line')>0) || sum(strcmp(PLOT.types,'ERSP-individ')>0)
%                 Analysis.Yaxis.(genvarname('ERSP-line'))(chanx,condx,1) = min(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%                 Analysis.Yaxis.(genvarname('ERSP-line'))(chanx,condx,2) = max(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%             end
%             if sum(strcmp(PLOT.types,'ITPC-line')>0)
%                 Analysis.Yaxis.(genvarname('ITPC-line'))(chanx,condx,1) = min(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%                 Analysis.Yaxis.(genvarname('ITPC-line'))(chanx,condx,2) = max(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%             end
%             if sum(strcmp(PLOT.types,'R2-line')>0)
%                 Analysis.Yaxis.(genvarname('R2-line'))(chanx,condx,1) = min(AVG.(genvarname(char(strcat('R2',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%                 Analysis.Yaxis.(genvarname('R2-line'))(chanx,condx,2) = max(AVG.(genvarname(char(strcat('R2',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%             end
%             if sum(strcmp(PLOT.types,'Induced-line')>0) || sum(strcmp(PLOT.types,'Induced-individ')>0)
%                 Analysis.Yaxis.(genvarname('Induced-line'))(chanx,condx,1) = min(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%                 Analysis.Yaxis.(genvarname('Induced-line'))(chanx,condx,2) = max(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%             end
%             if sum(strcmp(PLOT.types,'ERP')>0)
%                 Analysis.Yaxis.(genvarname('ERP'))(chanx,condx,1) = min(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%                 Analysis.Yaxis.(genvarname('ERP'))(chanx,condx,2) = max(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%             end
%         end
%     end
%     
%     % check min across channels
%     %     ceilfix = @(x)ceil(abs(x)).*sign(x);
%     ylim_array = [0 .0001 .0005 .001 .005 .01 .05 .1 .5 1 5:5:100 200 500 1000 2000 5000 1000 5000 10000 25000 50000 100000];
% %     prod_array = repmat([5 2],1,100);
% %     ylim_array = zeros(1,length(prod_array)+1);
% %     start_val = .0001;
% %     for i=1:length(prod_array)
% %     ylim_array(i+1) = start_val*(prod(prod_array(1:i)));
% %     end
%     
%     if sum(strcmp(PLOT.types,'ERSP-line')>0) || sum(strcmp(PLOT.types,'ERSP-individ')>0)
%         Analysis.Yaxis.(genvarname('ERSP-line'))(1) = min(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%         Analysis.Yaxis.(genvarname('ERSP-line'))(2) = max(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%     end
%     if sum(strcmp(PLOT.types,'ITPC-line')>0)
%         Analysis.Yaxis.(genvarname('ITPC-line'))(chanx,condx,1) = min(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%         Analysis.Yaxis.(genvarname('ITPC-line'))(chanx,condx,2) = max(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%     end
%     if sum(strcmp(PLOT.types,'R2-line')>0)
%         Analysis.Yaxis.(genvarname('R2-line'))(chanx,condx,1) = min(AVG.(genvarname(char(strcat('R2',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%         Analysis.Yaxis.(genvarname('R2-line'))(chanx,condx,2) = max(AVG.(genvarname(char(strcat('R2',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%     end
%     if sum(strcmp(PLOT.types,'Induced-line')>0) || sum(strcmp(PLOT.types,'Induced-individ')>0)
%         Analysis.Yaxis.(genvarname('Induced-line'))(chanx,condx,1) = min(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%         Analysis.Yaxis.(genvarname('Induced-line'))(chanx,condx,2) = max(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(condx)),'_avg'))))(plot_chans(chanx),xaxis));
%     end
%     if sum(strcmp(PLOT.types,'ERP')>0)
%         Analysis.Ylim.(genvarname('ERP'))(1) = min(Analysis.Yaxis.(genvarname('ERP'))(:,:,1));
%         Analysis.Ylim.(genvarname('ERP'))(2) = max(Analysis.Yaxis.(genvarname('ERP'))(:,:,2));
%         [c MinIndex] = min(abs(ylim_array-(abs(Analysis.Ylim.(genvarname('ERP'))(1)))));
%         if(ylim_array(MinIndex)<abs(Analysis.Ylim.(genvarname('ERP'))(1))),MinIndex=MinIndex+1;end
%         if(Analysis.Ylim.(genvarname('ERP'))(1)<0),
%             Analysis.Ylim.(genvarname('ERP'))(1)=-ylim_array(MinIndex);
%         else
%             Analysis.Ylim.(genvarname('ERP'))(1)=ylim_array(MinIndex);
%         end
%         [c MaxIndex] = min(abs(ylim_array-(abs(Analysis.Ylim.(genvarname('ERP'))(2)))));
%         if(ylim_array(MaxIndex)<abs(Analysis.Ylim.(genvarname('ERP'))(2))),MaxIndex=MaxIndex+1;end
%         if(Analysis.Ylim.(genvarname('ERP'))(2)<0),
%             Analysis.Ylim.(genvarname('ERP'))(2)=-ylim_array(MaxIndex);
%         else
%             Analysis.Ylim.(genvarname('ERP'))(2)=ylim_array(MaxIndex);
%         end
%     end
%     
%     
%     end
for chanx = fliplr(plot_chans) %[29 30 37 38 43:46 86:91]
    if sum(strcmp(PLOT.types,'ERSP-line')>0) || sum(strcmp(PLOT.types,'ERSP-individ')>0)
        zerobarERSP = zeros(1,trial_length);
        zerobarERSP(zerotime) = min(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis));
        zerobarERSP = zerobarERSP(xaxis);
    end
    if sum(strcmp(PLOT.types,'ITPC-line')>0)
        zerobarITPC = zeros(1,trial_length);
        zerobarITPC(zerotime) = min(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis));
        zerobarITPC = zerobarITPC(xaxis);
    end
    if sum(strcmp(PLOT.types,'R2-line')>0)
        zerobarR2 = zeros(1,trial_length);
        zerobarR2(zerotime) = min(AVG.(genvarname(char(strcat('R2',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis));
        zerobarR2 = zerobarR2(xaxis);
    end
    if sum(strcmp(PLOT.types,'Induced-line')>0) || sum(strcmp(PLOT.types,'Induced-individ')>0)
        zerobarInduced = zeros(1,trial_length);
        zerobarInduced(zerotime) = min(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis));
        zerobarInduced = zerobarInduced(xaxis);
    end
    if sum(strcmp(PLOT.types,'ERP')>0)
        zerobarERP = zeros(1,trial_length);
        zerobarERP(zerotime) = min(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis));
        zerobarERP = zerobarERP(xaxis);
    end
    
    f1=figure(100+chanx);
    if FullSrnFig==1
        set(0,'Units','pixels')
        screen_size = get(0,'ScreenSize');
        set(f1,'Position',[0 0 screen_size(3) screen_size(4)]);
    end
    
    analysescounter=1;
    while analysescounter <= length(PLOT.types)
        clear h
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'ERP'))>0
            subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:))
            if ismember(NumConds,1:2)
                h(1) = plot(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),'Color',colorarray(1,:),'LineWidth',LineThickness);
                hold on
                errbar = AVG.ERPstats(chanx,xaxis);
                shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),errbar,char(colorarrayShade(1)),.5);
                if NumConds>1
                    j=2;
                    h(j) = plot(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',colorarray(j,:),'LineWidth',LineThickness);
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),errbar,char(colorarrayShade(j)),.5);
                end
            end
            if ismember(NumConds,[3:4])
                tempmin = [];
                for j=1:NumConds
                    tempmin = [tempmin (min(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis)))];
                end
                tempmin = min(tempmin);
                bar([AVG.ERPstats(chanx,xaxis)*tempmin;AVG.ERPstats(chanx,xaxis)*tempmin]','stacked'), colormap([1 1 1;0 0 0])
                hold on
                for j=1:NumConds
                    h(j) = plot(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',colorarray(j,:),'LineWidth',LineThickness);
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err'))))(chanx,xaxis),char(colorarrayShade(j)),.5);
                end
            end
            if ismember(NumConds,[5:14])
                tempmin = [];
                for j=1:NumConds
                    tempmin = [tempmin (min(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis)))];
                end
                tempmin = min(tempmin);
                bar([AVG.ERPstats(chanx,xaxis)*tempmin;AVG.ERPstats(chanx,xaxis)*tempmin]','stacked'), colormap([1 1 1;0 0 0])
                hold on
                for j=1:NumConds
                    h(j) = plot(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',colorarray(j,:),'LineWidth',LineThickness);
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err'))))(chanx,xaxis),'k',.5);
                end
            end
            plot(zerobarERP)
            hold off
            title(strcat('Evoked'),'FontWeight','bold');
            set(gca, 'XTick', xtickrange)
            set(gca, 'XTickLabel', xticklabel)
            xlim([1 length(xaxis)])
            if isfield(SPECS,'Yaxis'),ylim(SPECS.Yaxis);end
            legend(h,legend_array,interpreter_string);
            analysescounter=analysescounter+1;
        end
        end
        
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'ERSP-line'))>0
           ha= subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:));
            
            if ismember(NumConds,1:2)
                h(1) = plot(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),'Color',colorarray(1,:),'LineWidth',LineThickness);
                hold on
                errbar = AVG.ERSPstats(chanx,xaxis);
                shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),errbar,char(colorarrayShade(1)),.5);
                if NumConds>1 
                    j=2;
                    h(j) = plot(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',colorarray(j,:),'LineWidth',LineThickness);
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),errbar,char(colorarrayShade(j)),.5);
                end
            end
            if ismember(NumConds,[3:4])
                tempmin = [];
                for j=1:NumConds
                    tempmin = [tempmin (min(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis)))];
                end
                tempmin = min(tempmin);
                bar([AVG.ERSPstats(chanx,xaxis)*tempmin;AVG.ERSPstats(chanx,xaxis)*tempmin]','stacked'), colormap([1 1 1;0 0 0])
                hold on
                for j=1:NumConds
                    h(j) = plot(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',colorarray(j,:),'LineWidth',LineThickness);
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_err'))))(chanx,xaxis),char(colorarrayShade(j)),.5);
                end
            end
            if ismember(NumConds,[5:14])
                tempmin = [];
                for j=1:NumConds
                    tempmin = [tempmin (min(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis)))];
                end
                tempmin = min(tempmin);
                hold on
                for j=1:NumConds
                    h(j) = plot(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',colorarray(j,:),'LineWidth',LineThickness);
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_err'))))(chanx,xaxis),'k',.5);
                end
            end
            
            
            
            plot(zerobarERSP)
            hold off
            title(strcat('Average ERSP:',{' '}, num2str(SPECS.Freq),{' '},'hz'),'FontWeight','bold');
            set(gca, 'XTick', xtickrange)
            set(gca, 'XTickLabel', xticklabel)
            xlim([1 length(xaxis)])
            xlabel('Time (sec)');
            ylabel('Normalized Amplitude');
            legend(h,legend_array,interpreter_string);
            analysescounter=analysescounter+1;
        end
        end
        
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'Distribution'))>0
           ha= subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:));   
           tempDist = []; for j=1:NumConds,tempDist(j) = (mean(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis)));end;
           h=bar([tempDist;zeros(length(tempDist),1)']);xlim([.5 1.5]);
           title(strcat('Average ERSP Dist:',{' '}, num2str(SPECS.Freq),{' '},'hz, between',{' '},num2str(SPECS.PLOTaxis(1)),{' '},'and',{' '},num2str(SPECS.PLOTaxis(2)),{' '},'seconds'),'FontWeight','bold');
           xlabel('Conditions');
           ylabel('Normalized Amplitude');
           legend(h,legend_array,interpreter_string);
           analysescounter=analysescounter+1;
        end
        end
        
        
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'ERSP-individ'))>0
            subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:))
            hold on
            for j=1:NumConds
                h(j) = plot(squeeze(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_indiv'))))(chanx,xaxis,1)),'Color',colorarray(j,:),'LineWidth',LineThickness);
            end
            legend(h,legend_array,interpreter_string);
            for j=1:NumConds
                plot(squeeze(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_indiv'))))(chanx,xaxis,:)),'Color',colorarray(j,:));
            end
            
            plot(zerobarERSP)
            hold off
            title(strcat('Average ERSP:',{' '}, num2str(SPECS.Freq),{' '},' hz'),'FontWeight','bold');
            set(gca, 'XTick', xtickrange)
            set(gca, 'XTickLabel', xticklabel)
            xlim([1 length(xaxis)])
            analysescounter=analysescounter+1;
        end
        end
        
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'ERSP'))>0
            clear tempmax tempmin
            for j=1:NumConds
                tempmax(j) = max(max(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_spect'))))(chanx,:,xaxis)));
                tempmin(j) = min(min(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_spect'))))(chanx,:,xaxis)));
            end
            templimERSP = round(10*max([tempmax abs(tempmin)]))/10;
            
            for j=1:NumConds
                subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:))
                contourf(time(xaxis),frex,squeeze(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_spect'))))(chanx,:,xaxis)),40,'linecolor','none')
                set(gca,'clim',[-templimERSP templimERSP])
                title(strcat('ERSP Cond: ',legend_array(j)),'FontWeight','bold');
                colorbar;
                analysescounter=analysescounter+1;
            end
        end    
        end
        
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'Induced-line'))>0
            subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:))
            if ismember(NumConds,1:2)
                h(1) = plot(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),'Color',colorarray(1,:),'LineWidth',LineThickness);
                hold on
                errbar = AVG.Inducedstats(chanx,xaxis);
                shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),errbar,char(colorarrayShade(1)),.5);
                if NumConds>1 
                    j=2;
                    h(j) = plot(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',colorarray(j,:),'LineWidth',LineThickness);
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),errbar,char(colorarrayShade(j)),.5);
                end
            end
            if ismember(NumConds,[3:4])
                tempmin = [];
                for j=1:NumConds
                    tempmin = [tempmin (min(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis)))];
                end
                tempmin = min(tempmin);
                bar([AVG.Inducedstats(chanx,xaxis)*tempmin;AVG.Inducedstats(chanx,xaxis)*tempmin]','stacked'), colormap([1 1 1;0 0 0])
                hold on
                for j=1:NumConds
                    h(j) = plot(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',colorarray(j,:),'LineWidth',LineThickness);
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_err'))))(chanx,xaxis),char(colorarrayShade(j)),.5);
                end
            end
            if ismember(NumConds,[5:14])
                tempmin = [];
                for j=1:NumConds
                    tempmin = [tempmin (min(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis)))];
                end
                tempmin = min(tempmin);
                bar([AVG.Inducedstats(chanx,xaxis)*tempmin;AVG.Inducedstats(chanx,xaxis)*tempmin]','stacked'), colormap([1 1 1;0 0 0])
                hold on
                for j=1:NumConds
                    h(j) = plot(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',colorarray(j,:),'LineWidth',LineThickness);
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_err'))))(chanx,xaxis),'k',.5);
                end
            end
            plot(zerobarInduced)
            hold off
            title(strcat('Average Induced:',{' '}, num2str(SPECS.Freq),{' '},' hz'),'FontWeight','bold');
            set(gca, 'XTick', xtickrange)
            set(gca, 'XTickLabel', xticklabel)
            xlim([1 length(xaxis)])
            legend(h,legend_array,interpreter_string);
            analysescounter=analysescounter+1;
        end
        end
        
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'Induced-individ'))>0
            subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:))
            hold on
            for j=1:NumConds
                h(j) = plot(squeeze(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_indiv'))))(chanx,xaxis,1)),'Color',colorarray(j,:));
            end
            legend(h,legend_array,interpreter_string);
            for j=1:NumConds
                plot(squeeze(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_indiv'))))(chanx,xaxis,:)),'Color',colorarray(j,:));
            end
            
            plot(zerobarInduced)
            hold off
            title(strcat('Average Induced:',{' '}, num2str(SPECS.Freq),{' '},' hz'),'FontWeight','bold');
            set(gca, 'XTick', xtickrange)
            set(gca, 'XTickLabel', xticklabel)
            xlim([1 length(xaxis)])
            analysescounter=analysescounter+1;
        end
        end
        
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'Induced'))>0
            clear tempmax tempmin
            for j=1:NumConds
                tempmax(j) = max(max(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_spect'))))(chanx,:,xaxis)));
                tempmin(j) = min(min(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_spect'))))(chanx,:,xaxis)));
            end
            templimInduced = round(10*max([tempmax abs(tempmin)]))/10;
            
            for j=1:NumConds
                subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:))
                contourf(time(xaxis),frex,squeeze(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_spect'))))(chanx,:,xaxis)),40,'linecolor','none')
                set(gca,'clim',[-templimInduced templimInduced])
                title(strcat('Induced Cond: ',legend_array(j)),'FontWeight','bold');
                colorbar;
                analysescounter=analysescounter+1;
            end
        end    
        end
        
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'ITPC-line'))>0
            subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:))
            
            
            if ismember(NumConds,1:2)
                h(1) = plot(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),'Color',colorarray(1,:),'LineWidth',LineThickness);
                hold on
                errbar = AVG.ITPCstats(chanx,xaxis);
                shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),errbar,char(colorarrayShade(1)),.5);
                if NumConds>1 
                    j=2;
                    h(j) = plot(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',colorarray(j,:),'LineWidth',LineThickness);
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),errbar,char(colorarrayShade(j)),.5);
                end
            end
            
            if ismember(NumConds,3:4)
                tempmin = [];
                for j=1:NumConds
                    tempmin = [tempmin (min(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis)))];
                end
                tempmin = min(tempmin);
                bar([AVG.ITPCstats(chanx,xaxis)*tempmin;AVG.ITPCstats(chanx,xaxis)*tempmin]','stacked'), colormap([1 1 1;0 0 0])
                hold on
                for j=1:NumConds
                    h(j) = plot(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',colorarray(j,:),'LineWidth',LineThickness);
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_err'))))(chanx,xaxis),char(colorarrayShade(j)),.5);
                end
            end
            if ismember(NumConds,5:14)
                tempmin = [];
                for j=1:NumConds
                    tempmin = [tempmin (min(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis)))];
                end
                tempmin = min(tempmin);
                bar([AVG.ITPCstats(chanx,xaxis)*tempmin;AVG.ITPCstats(chanx,xaxis)*tempmin]','stacked'), colormap([1 1 1;0 0 0])
                hold on
                for j=1:NumConds
                    h(j) = plot(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',colorarray(j,:),'LineWidth',LineThickness);
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_err'))))(chanx,xaxis),'k',.5);
                end
            end
            
            plot(zerobarITPC)
            hold off
            title(strcat('Average ITPC: ',{' '}, num2str(SPECS.Freq),{' '},' hz'),'FontWeight','bold');
            set(gca, 'XTick', xtickrange)
            set(gca, 'XTickLabel', xticklabel)
            xlim([1 length(xaxis)])
            legend(h,legend_array,interpreter_string);
            analysescounter=analysescounter+1;
        end
        end
        
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'ITPC'))>0
            clear tempmax tempmin
            for j=1:NumConds
                tempmax(j) = max(max(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_spect'))))(chanx,:,xaxis)));
                tempmin(j) = min(min(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_spect'))))(chanx,:,xaxis)));
            end
            templimERSP = round(10*max([tempmax abs(tempmin)]))/10;
            
            for j=1:NumConds
                subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:))
                contourf(time(xaxis),frex,squeeze(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)),'_spect'))))(chanx,:,xaxis)),40,'linecolor','none')
                set(gca,'clim',[-templimERSP templimERSP])
                title(strcat('ITPC Cond: ',legend_array(j)),'FontWeight','bold');
                colorbar;
                analysescounter=analysescounter+1;
            end
        end    
        end
        
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'R2-line'))>0
            subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:))
            if ismember(NumConds,1:2)
                h(1) = plot(AVG.(genvarname(char(strcat('R2',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),'Color',colorarray(1,:),'LineWidth',LineThickness);
                hold on
                errbar = AVG.R2stats(chanx,xaxis);
                shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('R2',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),errbar,char(colorarrayShade(1)),.5);
                if NumConds>1 
                    j=2;
                    h(j) = plot(AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',colorarray(j,:),'LineWidth',LineThickness);
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),errbar,char(colorarrayShade(j)),.5);
                end
            end
            if ismember(NumConds,[3:4])
                tempmin = [];
                for j=1:NumConds
                    tempmin = [tempmin (min(AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis)))];
                end
                tempmin = min(tempmin);
                bar([AVG.R2stats(chanx,xaxis)*tempmin;AVG.R2stats(chanx,xaxis)*tempmin]','stacked'), colormap([1 1 1;0 0 0])
                hold on
                for j=1:NumConds
                    h(j) = plot(AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',colorarray(j,:),'LineWidth',LineThickness);
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_err'))))(chanx,xaxis),char(colorarrayShade(j)),.5);
                end
            end
            if ismember(NumConds,[5:14])
                tempmin = [];
                for j=1:NumConds
                    tempmin = [tempmin (min(AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis)))];
                end
                tempmin = min(tempmin);
                bar([AVG.R2stats(chanx,xaxis)*tempmin;AVG.R2stats(chanx,xaxis)*tempmin]','stacked'), colormap([1 1 1;0 0 0])
                hold on
                for j=1:NumConds
                    h(j) = plot(AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',colorarray(j,:),'LineWidth',LineThickness);
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_err'))))(chanx,xaxis),char(colorarrayShade(j)),.5);
                end
            end
            
            plot(zerobarR2)
            hold off
            title(strcat('Average R2: ', num2str(SPECS.Freq),' hz'),'FontWeight','bold');
            set(gca, 'XTick', xtickrange)
            set(gca, 'XTickLabel', xticklabel)
            xlim([1 length(xaxis)])
            legend(h,legend_array,interpreter_string);
            analysescounter=analysescounter+1;
        end
        end
        
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'R2'))>0
            clear tempmax tempmin
            for j=1:NumConds
                tempmax(j) = max(max(AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_spect'))))(chanx,:,xaxis)));
                tempmin(j) = min(min(AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_spect'))))(chanx,:,xaxis)));
            end
            templimERSP = round(10*max([tempmax abs(tempmin)]))/10;
            
            for j=1:NumConds
                subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:))
                contourf(time(xaxis),frex,squeeze(AVG.(genvarname(char(strcat('R2',num2str(plot_conds(j)),'_spect'))))(chanx,:,xaxis)),40,'linecolor','none')
                set(gca,'clim',[-templimERSP templimERSP])
                title(strcat('R2 Cond: ',legend_array(j)),'FontWeight','bold');
                colorbar;
                analysescounter=analysescounter+1;
            end
        end    
        end
        
    end
    suptitle(strcat('Channel:',{' '}, chan_names(chanx)))
    set(gca, 'FontUnits', 'points', 'FontWeight','normal', 'FontSize', 14, 'FontName', 'Times')
%     xlabel('label','FontUnits','points', 'interpreter','latex', 'FontSize', 25, 'FontName', 'Times');
%     ylabel('Label','FontUnits','points', 'interpreter','latex', 'FontSize', 25, 'FontName', 'Times');


end
end



%Add statistics to subvariable "Statistics.(chans)" and mark on the plots
%where these intervals lie
if SPECS.DoStats ==1
    Sig_marks=zeros(length(chan_names),length(xaxis)); clear Significance
    sigbox = zeros(1,length(time));
    fq=EEGpnts/(epoch(2)-epoch(1));
    for j= SPECS.plot_chans
        rectcoords=[]; tempsigs=[];temptimepoints=[];sig_counter=1;
        ActiveFig = figure(j+100);
        yl = ylim;
        sig_intervals=1; Significance.(genvarname(['chan' num2str(j)])).timepointsMerged =[];
        for i = xaxis
            if NumConds ==2 %values for ttest
                if 1-abs(AVG.ERSPP(j,i)) <=.05
                    tempsigs(sig_counter)=1-abs(AVG.ERSPP(j,i));
                    temptimepoints(sig_counter) = i;
                    sig_counter=sig_counter+1;
                else
                    if sig_counter >= 10
                        Significance.(genvarname(['chan' num2str(j)])).(genvarname(['pvals' num2str(sig_intervals)]))  = tempsigs;
                        Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)])) = temptimepoints;
                        
                        Significance.(genvarname(['chan' num2str(j)])).(genvarname(['pvals' num2str(sig_intervals) '_Mean'])) = mean(Significance.(genvarname(['chan' num2str(j)])).(genvarname(['pvals' num2str(sig_intervals)])));
                        Significance.(genvarname(['chan' num2str(j)])).(genvarname(['pvals' num2str(sig_intervals) '_Min'])) = min(Significance.(genvarname(['chan' num2str(j)])).(genvarname(['pvals' num2str(sig_intervals)])));
                        
                        Significance.(genvarname(['chan' num2str(j)])).(genvarname(['time_ms' num2str(sig_intervals)])) = ...
                            [(Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)]))...
                            (1,1)/fq+epoch(1)), (Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)]))(1,end)/fq+epoch(1))];
                        
                        %                 rectcoords(sig_intervals,:) = [Significance.(genvarname(['chan' num2str(j)])).(genvarname(['time_ms' num2str(sig_intervals)]))(1), yl(2) , Significance.(genvarname(['chan' num2str(j)])).(genvarname(['time_ms' num2str(sig_intervals)]))(2)-Significance.(genvarname(['chan' num2str(j)])).(genvarname(['time_ms' num2str(sig_intervals)]))(1), yl(1)];
                        %                 rectcoords(sig_intervals,:) = [Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)]))(1), yl(2) , ((Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)]))(end))-(Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)]))(1))), yl(1)];
                        rectcoords(sig_intervals,:) = [Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)]))(1), yl(2) , Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)]))(end), yl(1)];
                        
                        sig_intervals = sig_intervals+1;
                    end
                    
                    tempsigs = [];
                    temptimepoints=[];
                    sig_counter=1;
                end
            end
            if NumConds>2 %values for ANOVA
                if AVG.ERSPP(j,i) <= .05
                    tempsigs(sig_counter)=1-abs(AVG.ERSPP(j,i));
                    temptimepoints(sig_counter) = i;
                    sig_counter=sig_counter+1;
                else
                    if sig_counter >= 10
                        Significance.(genvarname(['chan' num2str(j)])).(genvarname(['pvals' num2str(sig_intervals)]))  = tempsigs;
                        Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)])) = temptimepoints;
                        
                        Significance.(genvarname(['chan' num2str(j)])).(genvarname(['pvals' num2str(sig_intervals) '_Mean'])) = mean(Significance.(genvarname(['chan' num2str(j)])).(genvarname(['pvals' num2str(sig_intervals)])));
                        Significance.(genvarname(['chan' num2str(j)])).(genvarname(['pvals' num2str(sig_intervals) '_Min'])) = min(Significance.(genvarname(['chan' num2str(j)])).(genvarname(['pvals' num2str(sig_intervals)])));
                        
                        Significance.(genvarname(['chan' num2str(j)])).(genvarname(['time_ms' num2str(sig_intervals)])) = ...
                            [(Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)]))...
                            (1,1)/fq+epoch(1)), (Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)]))(1,end)/fq+epoch(1))];
                        
                        %                 rectcoords(sig_intervals,:) = [Significance.(genvarname(['chan' num2str(j)])).(genvarname(['time_ms' num2str(sig_intervals)]))(1), yl(2) , Significance.(genvarname(['chan' num2str(j)])).(genvarname(['time_ms' num2str(sig_intervals)]))(2)-Significance.(genvarname(['chan' num2str(j)])).(genvarname(['time_ms' num2str(sig_intervals)]))(1), yl(1)];
                        %                 rectcoords(sig_intervals,:) = [Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)]))(1), yl(2) , ((Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)]))(end))-(Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)]))(1))), yl(1)];
                        rectcoords(sig_intervals,:) = [Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)]))(1), yl(2) , Significance.(genvarname(['chan' num2str(j)])).(genvarname(['timepoints' num2str(sig_intervals)]))(end), yl(1)];
                        
                        sig_intervals = sig_intervals+1;
                    end
                    
                    tempsigs = [];
                    temptimepoints=[];
                    sig_counter=1;
                end
            end 
            
            
        end
        
        %Draw significant rectangles on figures with annotations
        %These work for the figures that are output for me (JZ) in R2014b, but
        %I'm not positive that they'll work for everyone. We might need a
        %better solution down the line...
        hold on
        MinX = .12957; MaxX=.909928571428571; MinY=.1; HeightY = .73333;
        MinXT=min(xaxis); MaxXT = max(xaxis);
        for i = 1:size(rectcoords,1)
            
            xCoord = rectcoords(i,1);
            xMaxCoord = rectcoords(i,3);
            scaledpointX = (xCoord-MinXT)*(MaxX-MinX)/(MaxXT-MinXT) + MinX;
            scaledpointWidth = ((xMaxCoord-MinXT)*(MaxX-MinX)/(MaxXT-MinXT) + MinX)-scaledpointX;
            h = annotation('rectangle','position',[scaledpointX .1 scaledpointWidth HeightY], 'FaceColor','k','FaceAlpha',.15,'LineStyle','none');
            hold on
        end
        hold off
        
        
    end
end


%Print Images to file
if SPECS.PrintFigs ==1 
    cd(SPECS.imgdir);
    for j=SPECS.plot_chans
        ActiveFig = figure(j+100);
        cd(SPECS.imgdir);
        printvariable = strcat(SPECS.imglabels, '_', chan_names(j), '.eps');
        print(printvariable{1}, '-depsc2', '-loose')
    end
end

 

