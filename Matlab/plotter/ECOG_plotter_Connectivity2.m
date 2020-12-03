% ECOG_plotter

% Enter condition numbers interested in plotting
time = epoch(1):1/srate:epoch(2);
xaxis = dsearchn(time',SPECS.PLOTaxis');
xaxis = xaxis(1):xaxis(2);
zerotime=dsearchn(time',0);
plotbaseline = dsearchn(time',SPECS.ERPbaseline');
plotbaseline = plotbaseline(1):plotbaseline(2);
zbaseline = dsearchn(time',SPECS.Powerbaseline');
zbaseline = zbaseline(1):zbaseline(2);

xticksize=SPECS.xticksize;
% % % % % change_yaxis=SPECS.change_yaxis;
% % % % % yaxis_ERP=SPECS.yaxis_ERP;
% % % % % yaxis_C=SPECS.yaxis_C;

db_on=SPECS.db_on;
plot_conds=SPECS.plot_conds;
plot_chansA=SPECS.plot_chansA;
plot_chansB=SPECS.plot_chansB;

NPerm=SPECS.NPerm;

NumConds = length(plot_conds);

xticklabel = SPECS.PLOTaxis(1):xticksize:SPECS.PLOTaxis(2);
xtickrange = 1:srate*xticksize:trial_length;
colorarray = {'r','g','b','y'};
if exist('bad_trials') == 0, bad_trials = [];end

% define wavelet parameters
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

ISPCwindow = (cyclex)./frex;
WindowPts = round(srate*ISPCwindow);

% create legend
legend_array = {};

for j=1:NumConds
    temp = find(CondNums==plot_conds(j));
    legend_array(length(legend_array)+1) = Conditions(temp);
end
clear AVG
counter = 1;

%Preallocate up to the max size of any var, then make sure to pull only
%relevant info
NTrials=nan(length(plot_chansA),length(plot_chansB),NumConds);
for chanxA = plot_chansA %[29 30 37 38 43:46 86:91]
    for chanxB = plot_chansB %[29 30 37 38 43:46 86:91]
    for j=1:NumConds
        tempsizeERP=0;
        for i=1:numevents
            if Evoked_chan_bad_array(chanxA,i)==0 && ~ismember(i,bad_trials) && Evoked_chan_bad_array(chanxB,i)==0
                if stim(i,2)==plot_conds(j)
                    tempsizeERP=tempsizeERP+1;
                end
            end
        end
        NTrials(chanxA,chanxB,j) = tempsizeERP;
    end
    end
end


MaxChansA = length(plot_chansA);
MaxChansB = length(plot_chansB);
MaxTrials = max(max(NTrials));

for j=1:NumConds
%     AVG.(genvarname(char(strcat('ISPCtrials',num2str(plot_conds(j)))))) = nan(MaxChans,MaxChansB,num_frex,trial_length,MaxTrials);
%     AVG.(genvarname(char(strcat('ISPCtime',num2str(plot_conds(j)))))) = nan(MaxChans,MaxChansB,num_frex,trial_length,MaxTrials);
    AVG.(genvarname(char(strcat('ISPCtrials',num2str(plot_conds(j)),'_avg')))) = nan(MaxChansA,MaxChansB,num_frex,trial_length);
    AVG.(genvarname(char(strcat('ISPCtime',num2str(plot_conds(j)),'_avg')))) = nan(MaxChansA,MaxChansB,num_frex,trial_length);
    AVG.(genvarname(char(strcat('POWERtrials',num2str(plot_conds(j)),'_avg')))) = nan(MaxChansA,MaxChansB,num_frex,trial_length);
    AVG.(genvarname(char(strcat('POWERtime',num2str(plot_conds(j)),'_avg')))) = nan(MaxChansA,MaxChansB,num_frex,trial_length);
end
if SPECS.matlabpool>0
    matlabpool('open','local',SPECS.matlabpool)
end
for chanxA = plot_chansA %[29 30 37 38 43:46 86:91]
    for chanxB = plot_chansB %[29 30 37 38 43:46 86:91]
        disp(['Channel ' int2str(counter) ' of ' int2str(length(plot_chansA)*length(plot_chansB))])
        clear PLOT
        for j=1:NumConds
            tempsizeERP=1;
            for i=1:numevents
                if Evoked_chan_bad_array(chanxA,i)==0 && ~ismember(i,bad_trials) && Evoked_chan_bad_array(chanxB,i)==0
                    if stim(i,2)==plot_conds(j)
                        temp = (squeeze(evoked_sig(chanxA,:,i)));
                        PLOT.(genvarname(char(strcat('ERP_A',num2str(plot_conds(j))))))(:,tempsizeERP) = temp;
                        temp = (squeeze(evoked_sig(chanxB,:,i)));
                        PLOT.(genvarname(char(strcat('ERP_B',num2str(plot_conds(j))))))(:,tempsizeERP) = temp;
                        tempsizeERP=tempsizeERP+1;
                    end
                end
            end
            for i=1:size(PLOT.(genvarname(char(strcat('ERP_A',num2str(plot_conds(j)))))),2)
                PLOT.(genvarname(char(strcat('ERP_A',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP_A',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP_A',num2str(plot_conds(j))))))(:,i),1));
                PLOT.(genvarname(char(strcat('ERP_B',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP_B',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP_B',num2str(plot_conds(j))))))(:,i),1));
            end
        end
        
        % freq analysis
        for j=1:NumConds
            disp(char(strcat('Condition',{' '},int2str(j),{' '},'of',{' '},int2str(NumConds))))
            
            tempA = PLOT.(genvarname(char(strcat('ERP_A',num2str(plot_conds(j))))));
            tempB = PLOT.(genvarname(char(strcat('ERP_B',num2str(plot_conds(j))))));
            
            EEGpnts = size(tempA,1); % same for both A and B
            EEGtrials = size(tempA,2); % same for both A and B
            
            % definte convolution parameters
            n_wavelet            = length(time);
            n_data               = EEGpnts*EEGtrials;
            n_convolution        = n_wavelet+n_data-1;
            half_of_wavelet_size = (n_wavelet-1)/2;
            
            % get FFT of data
            eegfftA = fft(reshape(tempA,1,EEGpnts*EEGtrials),n_convolution);
            eegfftB = fft(reshape(tempB,1,EEGpnts*EEGtrials),n_convolution);
            
            % initialize
            powerA = zeros(EEGpnts,EEGtrials); % frequencies X time X trials
            powerB = zeros(EEGpnts,EEGtrials); % frequencies X time X trials
            angleA = zeros(EEGpnts,EEGtrials); % frequencies X time X trials
            angleB = zeros(EEGpnts,EEGtrials); % frequencies X time X trials
            
            % loop through frequencies and compute synchronization
            s_val = cyclex./(2*pi*frex);
            
            tempDat = [];
            for fi=1:num_frex
                disp(strcat('Frequency',{' '},int2str(fi),{' '},'of',{' '},int2str(num_frex)))
                
                wavelet_fft = fft( exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s_val(fi)^2))) ,n_convolution);
                
                convolution_result_fftA = ifft(wavelet_fft.*eegfftA,n_convolution);
                convolution_result_fftA = convolution_result_fftA(half_of_wavelet_size+1:end-half_of_wavelet_size);
                sigA = reshape(convolution_result_fftA,EEGpnts,EEGtrials);
                
                convolution_result_fftB = ifft(wavelet_fft.*eegfftB,n_convolution);
                convolution_result_fftB = convolution_result_fftB(half_of_wavelet_size+1:end-half_of_wavelet_size);
                sigB = reshape(convolution_result_fftB,EEGpnts,EEGtrials);
                
                powerA = sigA.*conj(sigA); % same as abs(sig1).^2;
                powerB = sigB.*conj(sigB); % same as abs(sig1).^2;
                angleA = angle(sigA);
                angleB = angle(sigB);
                
                % calculate Inter Signals Phase Coh - TRIALS
                temp = abs(mean(exp(1i*(angleA-angleB)),2))';
                meanx = mean(temp(plotbaseline));
                AVG.(genvarname(char(strcat('ISPCtrials',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,fi,:) = (temp-meanx);
                
                % calculate Inter Signals Phase Coh - TIME
                tempISPCtime = zeros(1,EEGpnts);
                for k=xaxis % time range
                    timex = [k-WindowPts(fi):k+WindowPts(fi)];
                    tempA = angleA(timex,:);
                    tempB = angleB(timex,:);
                    temp = abs(mean(exp(1i*(tempA-tempB)'),2));
                    tempISPCtime(k) = mean(temp);
                end
                meanx = mean(tempISPCtime(plotbaseline));
                AVG.(genvarname(char(strcat('ISPCtime',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,fi,:) = (tempISPCtime-meanx);
            
                % calculate Inter Signals Power Time - Correlation between electrodes across time window
                tempPOWERtime = zeros(EEGtrials,EEGpnts);
                for trialx = 1:EEGtrials
                    parfor k=xaxis % time range
                        tempPOWERtime(trialx,k) = corr(powerA(k-WindowPts(fi):k+WindowPts(fi),trialx),powerB(k-WindowPts(fi):k+WindowPts(fi),trialx),'type','s');
                    end
                end
                meanPow = mean(tempPOWERtime,1);
                meanx = mean(meanPow(plotbaseline));
                AVG.(genvarname(char(strcat('POWERtime',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,fi,:) = (meanPow-meanx);
                
                % calculate Inter Signals Power Trials - Correlation across electrodes at each time point
                tempPOWERtrials = zeros(1,EEGpnts);
                for k=xaxis % time range
                    tempA = powerA(k,:);
                    tempB = powerB(k,:);
                    tempPOWERtrials(k)=corr(tempA',tempB','type','s');
                end
                meanx = mean(tempPOWERtrials(plotbaseline));
                AVG.(genvarname(char(strcat('POWERtrials',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,fi,:) = (tempPOWERtrials-meanx);
                
                
                
% % % % %                 % calulate Spectral Coherence
% % % % %                 tempAphase = (squeeze(mphaseA(:,i,:)));
% % % % %                 tempBphase = (squeeze(mphaseB(:,i,:)));
% % % % %                 tempApow = (squeeze(mpowA(:,i,:)));
% % % % %                 tempBpow = (squeeze(mpowB(:,i,:)));
% % % % %                 
% % % % %                 %spec1 = mean(abs(sig1).^2,2);
% % % % %                 %spec2 = mean(abs(sig2).^2,2);
% % % % %                 %specX = abs(mean( abs(sig1).*abs(sig2) .* exp(1i*(angle(sig1)-angle(sig2))) ,2)).^2;
% % % % %                 
% % % % %                 % compute spectral coherence, using only requested time points
% % % % %                 spectcoher(fi,:) = specX(times2saveidx)./(spec1(times2saveidx).*spec2(times2saveidx));
% % % % %                 
            end
            
            
        end
        counter = counter+1;
    end
end
if SPECS.matlabpool>0
    matlabpool('close')
end
disp('Plotting Channels')
for chanxA = plot_chansA %[29 30 37 38 43:46 86:91]
    for chanxB = plot_chansB %[29 30 37 38 43:46 86:91]
        clear tempmax tempmin
        for j=1:NumConds
            tempmax(j) = max(max(AVG.(genvarname(char(strcat('ISPCtrials',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,:,xaxis)));
            tempmin(j) = min(min(AVG.(genvarname(char(strcat('ISPCtrials',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,:,xaxis)));
        end
        templimISPCtrials = ceil(max([tempmax abs(tempmin)])*10)/10;
        
        clear tempmax tempmin
        for j=1:NumConds
            tempmax(j) = max(max(AVG.(genvarname(char(strcat('ISPCtime',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,:,xaxis)));
            tempmin(j) = min(min(AVG.(genvarname(char(strcat('ISPCtime',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,:,xaxis)));
        end
        templimISPCtime = ceil(max([tempmax abs(tempmin)])*10)/10;
        
        clear tempmax tempmin
        for j=1:NumConds
            tempmax(j) = max(max(AVG.(genvarname(char(strcat('POWERtime',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,:,xaxis)));
            tempmin(j) = min(min(AVG.(genvarname(char(strcat('POWERtime',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,:,xaxis)));
        end
        templimPOWERtime = ceil(max([tempmax abs(tempmin)])*10)/10;
        
        clear tempmax tempmin
        for j=1:NumConds
            tempmax(j) = max(max(AVG.(genvarname(char(strcat('POWERtrials',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,:,xaxis)));
            tempmin(j) = min(min(AVG.(genvarname(char(strcat('POWERtrials',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,:,xaxis)));
        end
        templimPOWERtrials = ceil(max([tempmax abs(tempmin)])*10)/10;
        
        f1=figure;
        set(0,'Units','pixels')
        screen_size = get(0,'ScreenSize');
        set(f1,'Position',[0 0 screen_size(3) screen_size(4)]);
        for j=1:NumConds
            subplot(4,NumConds,j)
            contourf(time(xaxis),frex,squeeze(AVG.(genvarname(char(strcat('ISPCtrials',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,:,xaxis)),40,'linecolor','none')
            set(gca,'clim',[-templimISPCtrials templimISPCtrials])
            title(strcat('ISPC Trials Cond: ',legend_array(j)),'FontWeight','bold');
            
            subplot(4,NumConds,j+NumConds)
            contourf(time(xaxis),frex,squeeze(AVG.(genvarname(char(strcat('ISPCtime',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,:,xaxis)),40,'linecolor','none')
            set(gca,'clim',[-templimISPCtime templimISPCtime])
            title(strcat('ISPC Time Cond: ',legend_array(j)),'FontWeight','bold');
            
            subplot(4,NumConds,j+(2*NumConds))
            contourf(time(xaxis),frex,squeeze(AVG.(genvarname(char(strcat('POWERtrials',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,:,xaxis)),40,'linecolor','none')
            set(gca,'clim',[-templimPOWERtrials templimPOWERtrials])
            title(strcat('Power Trials Cond: ',legend_array(j)),'FontWeight','bold');
            
            subplot(4,NumConds,j+(3*NumConds))
            contourf(time(xaxis),frex,squeeze(AVG.(genvarname(char(strcat('POWERtime',num2str(plot_conds(j)),'_avg'))))(chanxA,chanxB,:,xaxis)),40,'linecolor','none')
            set(gca,'clim',[-templimPOWERtime templimPOWERtime])
            title(strcat('Power Time Cond: ',legend_array(j)),'FontWeight','bold');
        end
        suptitle(strcat('Connectivity between',{' '}, chan_names(chanxA),{' '},'and',{' '},chan_names(chanxB)))
    end
end

