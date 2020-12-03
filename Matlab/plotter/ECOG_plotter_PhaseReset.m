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
plot_chans=SPECS.plot_chans;

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


% create legend
legend_array = {};

for j=1:NumConds
    temp = find(CondNums==plot_conds(j));
    legend_array(length(legend_array)+1) = Conditions(temp);
end
clear AVG PLOT
counter = 1;

%Preallocate up to the max size of any var, then make sure to pull only
%relevant info
NTrials=[];
for chanx = plot_chans
    for j=1:NumConds
        tempsizeERP=0;
        for i=1:numevents
            if Evoked_chan_bad_array(chanx,i)==0 && ~ismember(i,bad_trials)
                if stim(i,2)==plot_conds(j)
                    tempsizeERP=tempsizeERP+1;
                end
            end
        end
        NTrials(chanx,j) = tempsizeERP;
    end
end
MaxChans = size(NTrials,1);
MaxTrials = max(max(NTrials));
for j=1:NumConds
    AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))) = nan(MaxChans,trial_length,MaxTrials);
    AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))) = nan(MaxChans,trial_length);
    AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg')))) = nan(MaxChans,trial_length);
    AVG.(genvarname(char(strcat('C',num2str(plot_conds(j)))))) = nan(MaxChans,num_frex,trial_length,MaxTrials);
    AVG.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_avg')))) = nan(MaxChans,num_frex,trial_length);
    AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)))))) = nan(MaxChans,num_frex,trial_length);
end

for chanx = plot_chans %[29 30 37 38 43:46 86:91]
    disp(['Channel ' int2str(counter) ' of ' int2str(length(plot_chans))])
    for j=1:NumConds
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))) = [];
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))) = [];
        PLOT.(genvarname(char(strcat('ITPC',num2str(plot_conds(j)))))) = [];
    end
    
    NTrials = [];
    for j=1:NumConds
         
        tempsizeERP=1;
        for i=1:numevents
            if Evoked_chan_bad_array(chanx,i)==0 && ~ismember(i,bad_trials)
                if stim(i,2)==plot_conds(j)
                    temp = (squeeze(evoked_sig(chanx,:,i)));
                    PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,tempsizeERP) = temp;
                    tempsizeERP=tempsizeERP+1;
                end
            end 
        end
        NTrials(j) = tempsizeERP-1;
        
        for i=1:NTrials(j)
            PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(plotbaseline,i),1));
        end
        AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(chanx,:,1:NTrials(j)) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))));
        AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err'))))(chanx,:) = std(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))')/sqrt(NTrials(j));
        AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,:) = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)';
    end
    
    % Power calculated on freq averages
    for j=1:NumConds
        temp = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))); % time x trials format
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
        
        tempDat = [];
        for fi=1:num_frex
            wavelet_fft = fft( exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s_val(fi)^2))) ,n_convolution);
            
            convolution_result_fft = ifft(wavelet_fft.*eegfft,n_convolution);
            convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
            sig1 = reshape(convolution_result_fft,EEGpnts,EEGtrials);
            
            power1 = sig1.*conj(sig1); % same as mean(abs(sig1).^2,2);
            PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(fi,:,:) = power1;
            
            angle1 = angle(sig1);
            PLOT.(genvarname(char(strcat('ITPC',num2str(plot_conds(j))))))(fi,:,:) = angle1;
            % tempITPC = abs(mean(exp(1i*angle1),2))'; %ITPC
            tempITPC = squeeze(abs(mean(exp(1i*(2*angle1'))))); % R2
            
            % mean and norm of ITPC
            tempITPC = tempITPC-mean(tempITPC(plotbaseline));
            AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j))))))(chanx,fi,:) = tempITPC;
            
            % mean and norm of ERSP
            meanDAT = mean(power1,2)';
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
            
        AVG.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_avg'))))(chanx,:,:) = tempDat;
        AVG.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(chanx,:,:,1:NTrials(j)) = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))));

    end
    counter = counter+1;
end
% % % % % disp('Saving Data')
% % % % % save(strcat(subject,'_',task,'_ERSP'),'-v7.3', '-regexp', '^(?!(Evoked|evoked_sig|tempdata)$).');
disp('Plotting Channels')
for chanx = plot_chans %[29 30 37 38 43:46 86:91]
    clear tempmax tempmin
    for j=1:NumConds
        tempmax(j) = max(max(AVG.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_avg'))))(chanx,:,xaxis)));
        tempmin(j) = min(min(AVG.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_avg'))))(chanx,:,xaxis)));
    end
    templimERSP = round(max([tempmax abs(tempmin)]));
    
    clear tempmax tempmin
    for j=1:NumConds
        tempmax(j) = max(max(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j))))))(chanx,:,xaxis)));
        tempmin(j) = min(min(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j))))))(chanx,:,xaxis)));
    end
    templimITPC = ceil(max([tempmax abs(tempmin)])*10)/10;
    
    zerobarERP = zeros(1,trial_length);
    zerobarERP(zerotime) = min(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis));
    zerobarERP = zerobarERP(xaxis);
    
    
    f1=figure(100+chanx);
    subplot(NumConds,4,[1:4:NumConds*4 2:4:NumConds*4])
%     subplot(NumConds,3,[1:3:NumConds*3])
    h(1) = plot(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),'Color',char(colorarray(1)));
    hold on
    shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_avg'))))(chanx,xaxis),AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_err'))))(chanx,xaxis),colorarray(1),.5);
    if NumConds>1
        for j=2:NumConds
            h(j) = plot(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),'Color',char(colorarray(j)));
            shadedErrorBar(1:length(xaxis),AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_avg'))))(chanx,xaxis),AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err'))))(chanx,xaxis),colorarray(j),.5);
        end
    end
    plot(zerobarERP)
    hold off
    title(strcat('Evoked'),'FontWeight','bold');
    set(gca, 'XTick', xtickrange)
    set(gca, 'XTickLabel', xticklabel)
    %     xlim([xaxis(1) xaxis(2)])
    xlim([1 length(xaxis)])
% % % % %     if change_yaxis ==1
% % % % %         ylim([yaxis_ERP(1) yaxis_ERP(2)])
% % % % %     end
    legend(h,legend_array,'Location','NorthWest');
    set(0,'Units','pixels')
    screen_size = get(0,'ScreenSize');
    set(f1,'Position',[0 0 screen_size(3) screen_size(4)]);
    
    
    for j=1:NumConds
        subplot(NumConds,4, ((j-1)*4)+3)
        contourf(time(xaxis),frex,squeeze(AVG.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_avg'))))(chanx,:,xaxis)),40,'linecolor','none')
        set(gca,'clim',[-templimERSP templimERSP])
        title(strcat('ERSP Cond: ',legend_array(j)),'FontWeight','bold');
        
        subplot(NumConds,4,((j-1)*4)+4)
        contourf(time(xaxis),frex,squeeze(AVG.(genvarname(char(strcat('ITPC',num2str(plot_conds(j))))))(chanx,:,xaxis)),40,'linecolor','none')
        set(gca,'clim',[-templimITPC templimITPC])
        title(strcat('ITPC Cond: ',legend_array(j)),'FontWeight','bold');
    end 
    suptitle(strcat('CHANNEL:   ', chan_names(chanx)))
end
    
   
