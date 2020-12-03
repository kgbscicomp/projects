% ECOG_plotter

% Enter condition numbers interested in plotting

plotaxis=SPECS.plotaxis;
xaxis=SPECS.xaxis;
zerotime=SPECS.zerotime;
zbaseline=SPECS.zbaseline;
xticksize=SPECS.xticksize;
change_yaxis=SPECS.change_yaxis;
yaxis_ERP=SPECS.yaxis_ERP;
yaxis_C=SPECS.yaxis_C;
db_on=SPECS.db_on;
plot_conds=SPECS.plot_conds;
plot_chans=SPECS.plot_chans;

NumConds = length(plot_conds);
zerotime = zerotime-xaxis(1);
zbaseline = zbaseline;

zbaseline = [round((zbaseline(1)*srate)+1):round((zbaseline(2)*srate)+1)];
zerotime = round((zerotime*srate)+1);
xaxis = [round((xaxis(1)*srate)+1) round((xaxis(2)*srate)+1)];
xticklabel = epoch(1):xticksize:epoch(2);
xtickrange = 1:srate*xticksize:size(evoked_sig,2);
colorarray = {'r','g','b','y'};
% baseidx = dsearchn(time',[-.5 -.2]')

% define wavelet parameters
Freq = SPECS.Freq;
FreqWidth = SPECS.FreqWidth;
FreqCycles = SPECS.FreqCycles;
frex = Freq(1):FreqWidth:Freq(2);
time = (epoch(1)*srate:epoch(2)*srate)./srate;
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

if exist('bad_trials') == 0, bad_trials = [];end

% create legend
legend_array = {};

for j=1:NumConds
    temp = find(CondNums==plot_conds(j));
    legend_array(length(legend_array)+1) = Conditions(temp);
end

clear PLOT AVG
for j=1:NumConds
    PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))=[];
    AVG.(genvarname(char(strcat('C',num2str(plot_conds(j)))))) = [];
%     PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(j))))))=[];
end
for chanx = plot_chans %[29 30 37 38 43:46 86:91]
     
    for j=1:NumConds
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))) = [];
        PLOT.(genvarname(char(strcat('avgERP',num2str(plot_conds(j)))))) = [];
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
        
        for i=1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)
            PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i),1));
        end        
        
    end
    
    % freq analysis
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
        
        for fi=1:num_frex
            wavelet_fft = fft( exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s_val(fi)^2))) ,n_convolution);
            
            convolution_result_fft = ifft(wavelet_fft.*eegfft,n_convolution);
            convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
            sig1 = reshape(convolution_result_fft,EEGpnts,EEGtrials);
            power1 = mean(sig1.*conj(sig1),2); % same as mean(abs(sig1).^2,2);
            
            meanx = mean(power1(zbaseline,:));
            stdx = std(power1(zbaseline,:));
            if db_on==1
                PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(fi,:) = [10*log10(power1/meanx)]';
            else
                PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(fi,:) = [(power1-meanx)/stdx]';
            end
        end
        
        figure(10+j)
        contourf(time(xaxis(1):xaxis(2)),frex,PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,xaxis(1):xaxis(2)),40,'linecolor','none')
%         tempmax = max(max(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,xaxis(1):xaxis(2))));
%         tempmin = min(min(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,xaxis(1):xaxis(2))));
%         templim = max([tempmax abs(tempmin)]);
        set(gca,'xlim',plotaxis)
%         set(gca,'xlim',plotaxis,'clim',[-templim templim])
        title(strcat('CH:   ', chan_names(chanx),      ' Cond: ',legend_array(j)),'FontWeight','bold');
        
        AVG.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(chanx,:,:) = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))));
    end
    pause();
end

