% ECoG_plotter
% This function takes the output of the function ECOG_filter_analysis7
% along with the following inputs to plot the processed ECoG data:
%
% SPECS.PLOTaxis = [-1.2 1];
% SPECS.plot_conds = [304:305]; %[301 304:305]; %11 13 % Plot 1-4 conditions, with ordered colors of red, green, blue, orange
% RejThreshEvoked = 3;
% ECOG_epoch_cleaner3
% SPECS.ERPbaseline = [-1.2 -.9]; %[-.05 0]
% SPECS.Powerbaseline = [-1.2 -.9];
% SPECS.xticksize = .1;
% SPECS.db_on = 0; % Z if 0, dB if 1
% SPECS.plot_chans = [33];% 34 46:47];%1:length(chan_names);
% SPECS.NPerm=100; % Num permutations when needed
% SPECS.Freq = [70 110];%[5 40];
% SPECS.FreqWidth = [5];%2.5;
% SPECS.FreqCycles = [20];%[3];
% SPECS.bad_trials = []; % Bad trials by trial number (e.g. for blinks)
% SPECS.types = {'ERSP-line'};% SPECS.types = {'ERP','ERSP-line','ERSP','ERSP-individ','Induced-line','Induced','ITC-line','ITC','R2-line','R2'};
% ECOG_plotter_11

clear AVG PLOT
plot_conds=SPECS.plot_conds;
plot_chans=SPECS.plot_chans;
SPECS.Options.types = {'ERP','ERSP-line','ERSP','ERSP-individ','Induced-line','Induced','Induced-individ','ITC-line','ITC','R2-line','R2'};
SPECS.Options.plots = [1 1 .5 1 1 .5 1 1 .5 1 .5]; % how many of two rows will each plot take in the subplot
tempPLOTtypes = SPECS.types;
numanalyses = length(tempPLOTtypes);

PLOT.plotnum = 0;
PLOT.rowcount = length(plot_conds);
PLOT.colcount = 0;
PLOT.halfwhole = [];
for counter=1:numanalyses
    PLOT.colcount=PLOT.colcount+1;
    temp = SPECS.Options.plots(find(strcmp(tempPLOTtypes(counter),SPECS.Options.types)==1));
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
        for j=1:length(plot_conds)
            PLOT.types(counter2)=tempPLOTtypes(i);
            PLOT.plotlocs(counter2,:) = ones(1,PLOT.rowcount)*PLOT.plotarray(counter1);
            counter1=counter1+1;
            counter2=counter2+1;
        end    
    else
        temp = counter1:1:counter1+length(plot_conds)-1;
        PLOT.types(counter2)=tempPLOTtypes(i);
        PLOT.plotlocs(counter2,:) =PLOT.plotarray(temp);
        counter1=counter1+length(plot_conds);
        counter2=counter2+1;
    end
end

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

db_on=SPECS.db_on;
NPerm=SPECS.NPerm;
NumConds = length(plot_conds);
xticklabel = SPECS.PLOTaxis(1):xticksize:SPECS.PLOTaxis(2);
xtickrange = 1:srate*xticksize:trial_length;
colorarray = {'r','g','b','y'};
if exist('SPECS.bad_trials') == 0, bad_trials = [];else bad_trials = SPECS.bad_trials; end

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
    temp = find(CondNums==plot_conds(j));
    legend_array(length(legend_array)+1) = Conditions(temp);
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
            if Evoked_chan_bad_array(chanx,i)==0 && ~ismember(i,bad_trials)
                if stim(analyzed_trials(i),2)==plot_conds(j)
                    tempsizeERP=tempsizeERP+1;
                end
            end
            if ERSP_chan_bad_array(chanx,i)==0 && ~ismember(i,bad_trials)
                if stim(analyzed_trials(i),2)==plot_conds(j)
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
    end
    
end
if sum(strcmp(PLOT.types,'ERP')>0),AVG.ERPstats= nan(MaxChans,trial_length);end
if sum(strcmp(PLOT.types,'ERSP-line')>0),AVG.ERSPstats= nan(MaxChans,trial_length);end
if sum(strcmp(PLOT.types,'Induced-line')>0),AVG.Inducedstats= nan(MaxChans,trial_length);end

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
        end
    end
    
    NTrials_ERP = [];
    NTrials_ERSP = [];
    for j=1:NumConds
        
        tempsizeERP=1;
        tempsizeERSP=1;
        for i=1:length(analyzed_trials);
            if Evoked_chan_bad_array(chanx,i)==0 && ~ismember(i,bad_trials)
                if stim(analyzed_trials(i),2)==plot_conds(j)
                    temp = (squeeze(evoked_sig(chanx,:,analyzed_trials(i))));
                    PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,tempsizeERP) = temp;
                    tempsizeERP=tempsizeERP+1;
                end
            end
            if ERSP_chan_bad_array(chanx,i)==0 && ~ismember(i,bad_trials)
                if stim(analyzed_trials(i),2)==plot_conds(j)
                    temp = (squeeze(evoked_sig(chanx,:,analyzed_trials(i))));
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
        %subtract ERP from individual trials
        if sum(strcmp(PLOT.types,'Induced-line')>0) || sum(strcmp(PLOT.types,'Induced')>0) ||  sum(strcmp(PLOT.types,'Induced-individ')>0)
            for i=1:NTrials_ERSP(j)
                PLOT.(genvarname(char(strcat('Induced_ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j))))))(:,i)-AVG.(genvarname(char(strcat('ERSP_ERP',num2str(plot_conds(j)),'_avg'))))(chanx,:)';
            end
        end
    end
    
    % conduct stats on the ERP
    if NumConds==1 % 95% confidence interval
        [h,p,ci,stats] = ttest(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))))');
        tempMean = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2)';
        AVG.ERPstats(chanx,:) = abs(ci(1,:)-tempMean);
    elseif NumConds==2 % ttest
        tempdata1 = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))))';
        tempdata2 = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(2))))))';
        diff = mean(tempdata1,1)-mean(tempdata2,1);
        [h,p,ci,stats] = ttest2(tempdata1, tempdata2);
        AVG.ERPstats(chanx,:) = abs(ci(2,:)-diff)/2;
    elseif NumConds==3 % ANOVA
        tempgroup = [ones(1,NTrials_ERP(1)),ones(1,NTrials_ERP(2))+1,ones(1,NTrials_ERP(3))+2];
        for timex = 1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),1)
            tempdata1 = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))))(timex,:);
            tempdata2 = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(2))))))(timex,:);
            tempdata3 = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(3))))))(timex,:);
            tempstat=anova1([tempdata1,tempdata2,tempdata3],tempgroup,'off');
            if tempstat< .05
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
            if tempstat< .05
                AVG.ERPstats(chanx,timex)=1;
            else
                AVG.ERPstats(chanx,timex)=0;
            end
        end
    end
    
    
    
    if SPECS.Options.ERSP==1
        ERSPcond = {};
        if sum(strcmp(PLOT.types,'ERSP-line')>0) || sum(strcmp(PLOT.types,'ERSP')>0) ||  sum(strcmp(PLOT.types,'ERSP-individ')>0)
            ERSPcond = [ERSPcond 'ERSP'];
        end
        if sum(strcmp(PLOT.types,'Induced-line')>0) || sum(strcmp(PLOT.types,'Induced')>0) ||  sum(strcmp(PLOT.types,'Induced-individ')>0)
            ERSPcond = [ERSPcond 'Induced'];
        end
        
        
        for ERSPcLength=1:length(ERSPcond)
            % Power calculated on freq averages for ERSP
            for j=1:NumConds
                temp = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),'_ERP',num2str(plot_conds(j)))))); % time x trials format
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
                    power1 = sig1.*conj(sig1); % same as mean(abs(sig1).^2,2);
                    
                    % baseline norm the power to enable sign switching perm analyses
                    %                 power1 = power1 - repmat(mean(power1(zbaseline,:),1),size(power1,1),1);
                    PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(fi,:,:) = power1;
                end
            end
            
            tempFreqData = [];
            STDthresh = []; % in case we wanted to do a second level of trial rejs just on the induced. currently off (effectively)
            for fi=1:num_frex
                for j=1:NumConds
                    tempFreqData = [tempFreqData squeeze(PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(fi,:,:))];
                end
                STDthresh(fi) = std(std(tempFreqData(zbaseline(1):xaxis(end),:)))*20;
            end
            
            for j=1:NumConds
                C_bad_trials = [];
                for fi=1:num_frex
                    power1 = squeeze(PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(fi,:,:));
                    stdDAT = std(power1(zbaseline(1):xaxis(end),:));
                    for i=1:length(stdDAT)
                        if stdDAT(i)>STDthresh(fi)
                            C_bad_trials(fi,i) = 1;
                        else
                            C_bad_trials(fi,i) = 0;
                        end
                    end
                end
                [I,J] = find(C_bad_trials);
                C_good_trials = 1:size(C_bad_trials,2);
                C_good_trials(J)=[];
                
                tempDatIndiv = [];
                tempDat = [];
                for fi=1:num_frex
                    power1 = squeeze(PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(fi,:,C_good_trials));
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
                AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))))(chanx,:,:,1:NTrials_ERSP(j)) = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j))))));
                if sum(strcmp(PLOT.types,strcat(ERSPcond(ERSPcLength),'-individ'))>0),AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j)),'_indiv'))))(chanx,:,1:NTrials_ERSP(j)) = mean(tempDatIndiv,3)';end
                if sum(strcmp(PLOT.types,ERSPcond(ERSPcLength))>0),AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(j)),'_spect'))))(chanx,:,:) = tempDat;end
                
                if sum(strcmp(PLOT.types,strcat(ERSPcond(ERSPcLength),'-line')))>0
                    % run perms for stdev for the ERSP-line analysis
                    if ismember(NumConds,[1 3:4])
                        tempPerm = [];
                        for k=1:NPerm
                            permAvg = [];
                            permX = randi(NTrials_ERSP(j),[1 NTrials_ERSP(j)]);
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
                if NumConds==1 % 95% confidence interval
                    % One-Sample Permuation using Sign-Switching
                    % for now still use the 95% bootstrapped CIs above
                    % Bootstrapped gives similar values and easier to compute CIs
                    % % % % %
                    % % % % %         plotSize = size(PLOT.(genvarname(char(strcat('ERSP',num2str(plot_conds(1)))))));
                    % % % % %         realMean = AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(1)),'_avg'))))(chanx,:);
                    % % % % %         tempPermVals = zeros(NPerm,plotSize(2));
                    % % % % %         for k=1:NPerm
                    % % % % %             signx = randi(2,1,plotSize(3));
                    % % % % %             signx = reshape(signx,[1 1 60]);
                    % % % % %             signx = repmat(signx,plotSize(1),plotSize(2));
                    % % % % %             signx( signx==2 )=-1;
                    % % % % %             tempPerm = PLOT.(genvarname(char(strcat('ERSP',num2str(plot_conds(1)))))).*signx;
                    % % % % %             tempDat = [];
                    % % % % %             for fi=1:num_frex % calculate power for each freq
                    % % % % %                 meanDAT = mean(tempPerm(fi,:,:),3);
                    % % % % %                 meanC = mean(meanDAT(zbaseline));
                    % % % % %                 stdx = std(meanDAT(zbaseline));
                    % % % % %                 if db_on==1
                    % % % % %                     tempDB = [10*log10(meanDAT/meanC)]';
                    % % % % %                     tempDat(fi,:) = tempDB-mean(tempDB(plotbaseline));
                    % % % % %                 else
                    % % % % %                     tempDat(fi,:) = (meanDAT-meanC)/stdx;
                    % % % % %                     tempDat(fi,:) = tempDat(fi,:)-mean(tempDat(fi,plotbaseline));
                    % % % % %                 end
                    % % % % %             end
                    % % % % %             tempVal = mean(tempDat,1);
                    % % % % %             for j = 1:plotSize(2)
                    % % % % %                 if realMean(j) > 0
                    % % % % %                     if realMean(j)-tempVal(j)> 0
                    % % % % %                         tempPermVals(k,j) = tempPermVals(k,j)+1;
                    % % % % %                     end
                    % % % % %                 elseif realMean(j) < 0
                    % % % % %                     if realMean(j)-tempVal(j)< 0
                    % % % % %                         tempPermVals(k,j) = tempPermVals(k,j)+1;
                    % % % % %                     end
                    % % % % %                 end
                    % % % % %             end
                    % % % % %         end
                    % % % % %         tempMean = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2)';
                    % % % % %         AVG.ERPstats(chanx,:) = abs(ci(1,:)-tempMean);
                elseif NumConds==2 % ttest
                    
                    % randomly switch condition labels, resampling without replacement
                    % permutation analysis
                    
                    tempdata1 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1))))));
                    tempdata2 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2))))));
                    diff = AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1)),'_avg'))))(chanx,:)-AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2)),'_avg'))))(chanx,:);
                    plotSize = size(diff);
                    plotSize1 = size(tempdata1);
                    plotSize2 = size(tempdata2);
                    Analysis.PermTest = zeros(1,plotSize(2));
                    Analysis.PermVals = zeros(NPerm,plotSize(2));
                    
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
                                    Analysis.PermTest(j) = Analysis.PermTest(j)+1;
                                end
                            end
                        end
                        %             if ismember(k,round(linspace(1,NPerm,20)))
                        %                 (k/NPerm)
                        %             end
                    end
                    Analysis.PermTest = Analysis.PermTest/NPerm;
                    Analysis.PermTest = Analysis.PermTest';
                    
                    CutoffValmax = [quantile(1:NPerm,1-.05)+.5];
                    tempA = sort(Analysis.PermValues,1);
                    % AVG.ERSPstats(chanx,:) = tempA(CutoffValmax',:)/2;
                    AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'stats'))))(chanx,:) = tempA(CutoffValmax',:)/2;
                    
                elseif NumConds==3 % ANOVA
                    
                    % Calculate the sum of squared deviations of group means from the grand mean
                    % https://www.uvm.edu/~dhowell/StatPages/Resampling/Bootst1way/bootstrapping_oneway.html
                    % 1. find real: sum((mean of all Z - mean of each Z)^2)
                    
                    
                    tempdata1 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1))))));
                    tempdata2 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2))))));
                    tempdata3 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(3))))));
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
                    Analysis.PermVals = zeros(NPerm,plotSize(2));
                    
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
                    Analysis.PermTest = Analysis.PermTest/NPerm;
                    Analysis.PermTest = Analysis.PermTest;
                    for timex=1:plotSize(2)
                        if Analysis.PermTest(timex)< .05
                            AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'stats'))))(chanx,timex) = 1;
                        else
                            AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'stats'))))(chanx,timex) = 0;
                        end
                    end
                    
                    
                elseif NumConds==4 % ANOVA
                    
                    
                    tempdata1 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(1))))));
                    tempdata2 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(2))))));
                    tempdata3 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(3))))));
                    tempdata4 = PLOT.(genvarname(char(strcat(ERSPcond(ERSPcLength),num2str(plot_conds(4))))));
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
                    Analysis.PermVals = zeros(NPerm,plotSize(2));
                    
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
                    Analysis.PermTest = Analysis.PermTest/NPerm;
                    Analysis.PermTest = Analysis.PermTest;
                    for timex=1:plotSize(2)
                        if Analysis.PermTest(timex)< .05
                            AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'stats'))))(chanx,timex) =1;
                        else
                            AVG.(genvarname(char(strcat(ERSPcond(ERSPcLength),'stats'))))(chanx,timex) =0;
                        end
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

for chanx = plot_chans %[29 30 37 38 43:46 86:91]
    if sum(strcmp(PLOT.types,'ERSP-line')>0) || sum(strcmp(PLOT.types,'ERSP-individ')>0)
        zerobarERSP = zeros(1,trial_length);
        zerobarERSP(zerotime) = min(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis));
        zerobarERSP = zerobarERSP(xaxis);
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
    set(0,'Units','pixels')
    screen_size = get(0,'ScreenSize');
    set(f1,'Position',[0 0 screen_size(3) screen_size(4)]);
    
    analysescounter=1;
    while analysescounter <= length(PLOT.types)
        clear h
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'ERP'))>0
            subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:))
            if ismember(NumConds,1:2)
                h(1) = plot(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),'Color',char(colorarray(1)));
                hold on
                errbar = AVG.ERPstats(chanx,xaxis);
                shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),errbar,colorarray(1),.5);
                if NumConds>1
                    j=2;
                    h(j) = plot(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',char(colorarray(j)));
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),errbar,colorarray(j),.5);
                end
            end
            if NumConds>2
                bar((AVG.ERPstats(chanx,xaxis)*min(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis))),'FaceColor',[.8 .8 .8]);
                hold on
                for j=1:NumConds
                    h(j) = plot(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',char(colorarray(j)));
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err'))))(chanx,xaxis),colorarray(j),.5);
                end
                
            end
            plot(zerobarERP)
            hold off
            title(strcat('Evoked'),'FontWeight','bold');
            set(gca, 'XTick', xtickrange)
            set(gca, 'XTickLabel', xticklabel)
            xlim([1 length(xaxis)])
            legend(h,legend_array,'Location','NorthWest');
            analysescounter=analysescounter+1;
        end
        end
        
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'ERSP-line'))>0
            subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:))
            if ismember(NumConds,2)
                h(1) = plot(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),'Color',char(colorarray(1)));
                hold on
                errbar = AVG.ERSPstats(chanx,xaxis);
                shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),errbar,colorarray(1),.5);
                j=2;
                h(j) = plot(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',char(colorarray(j)));
                shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),errbar,colorarray(j),.5);
            end
            if ismember(NumConds,[1,3:4])
                bar((AVG.ERSPstats(chanx,xaxis)*min(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis))),'FaceColor',[.8 .8 .8]);
                hold on
                for j=1:NumConds
                    h(j) = plot(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',char(colorarray(j)));
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_err'))))(chanx,xaxis),colorarray(j),.5);
                end
            end
            
            plot(zerobarERSP)
            hold off
            title(strcat('Average ERSP: ', num2str(SPECS.Freq),' hz'),'FontWeight','bold');
            set(gca, 'XTick', xtickrange)
            set(gca, 'XTickLabel', xticklabel)
            xlim([1 length(xaxis)])
            legend(h,legend_array,'Location','NorthWest');
            analysescounter=analysescounter+1;
        end
        end
        
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'ERSP-individ'))>0
            subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:))
            hold on
            for j=1:NumConds
                h(j) = plot(squeeze(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_indiv'))))(chanx,xaxis,1)),'Color',char(colorarray(j)));
            end
            legend(h,legend_array,'Location','NorthWest');
            for j=1:NumConds
                plot(squeeze(AVG.(genvarname(char(strcat('ERSP',num2str(plot_conds(j)),'_indiv'))))(chanx,xaxis,:)),'Color',char(colorarray(j)));
            end
            
            plot(zerobarERSP)
            hold off
            title(strcat('Average ERSP: ', num2str(SPECS.Freq),' hz'),'FontWeight','bold');
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
            templimERSP = round(max([tempmax abs(tempmin)]));
            
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
            if ismember(NumConds,2)
                h(1) = plot(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),'Color',char(colorarray(1)));
                hold on
                errbar = AVG.Inducedstats(chanx,xaxis);
                shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),errbar,colorarray(1),.5);
                j=2;
                h(j) = plot(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',char(colorarray(j)));
                shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),errbar,colorarray(j),.5);
            end
            if ismember(NumConds,[1,3:4])
                bar((AVG.Inducedstats(chanx,xaxis)*min(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis))),'FaceColor',[.8 .8 .8]);
                hold on
                for j=1:NumConds
                    h(j) = plot(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',char(colorarray(j)));
                    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_err'))))(chanx,xaxis),colorarray(j),.5);
                end
            end
            
            plot(zerobarInduced)
            hold off
            title(strcat('Average Induced: ', num2str(SPECS.Freq),' hz'),'FontWeight','bold');
            set(gca, 'XTick', xtickrange)
            set(gca, 'XTickLabel', xticklabel)
            xlim([1 length(xaxis)])
            legend(h,legend_array,'Location','NorthWest');
            analysescounter=analysescounter+1;
        end
        end
        
        if analysescounter <= length(PLOT.types)
        if sum(strcmp(PLOT.types(analysescounter),'Induced-individ'))>0
            subplot(PLOT.rowcount,PLOT.colcount,PLOT.plotlocs(analysescounter,:))
            hold on
            for j=1:NumConds
                h(j) = plot(squeeze(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_indiv'))))(chanx,xaxis,1)),'Color',char(colorarray(j)));
            end
            legend(h,legend_array,'Location','NorthWest');
            for j=1:NumConds
                plot(squeeze(AVG.(genvarname(char(strcat('Induced',num2str(plot_conds(j)),'_indiv'))))(chanx,xaxis,:)),'Color',char(colorarray(j)));
            end
            
            plot(zerobarInduced)
            hold off
            title(strcat('Average Induced: ', num2str(SPECS.Freq),' hz'),'FontWeight','bold');
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
            templimInduced = round(max([tempmax abs(tempmin)]));
            
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
    end
    suptitle(strcat('CHANNEL:   ', chan_names(chanx)))
end
