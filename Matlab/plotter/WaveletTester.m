function WaveletTester(EEG,filter_range,filter_width,MorletWidth,variablefreq)

% Test Wavelet construction

% VARY wavelet width by frequency
% Set the wavelet width for the very first cycle then it will match the
% FWHM at all other frequencies -- 

srate = EEG.srate;

linefq = srate/2; % changed from 60 hz line activity to include all freqs
fqstep = filter_width; 
fqrng = filter_range; 
line_param = 6; 
numcycles = MorletWidth;

allfqs = (fqrng(1):fqstep:fqrng(end))';
badfqs = and(allfqs > line_param, abs(mod(allfqs,linefq) - linefq/2) > (linefq/2 - line_param));
allfqs = allfqs(~badfqs);

frequency = allfqs;
time      = -3:1/srate:3;

wavecolors = 'br';

% calculate FWHM for first wavelet for each cycle
numcycles_arry = [];
FWHM_arry = [];
for i=1:length(numcycles)
    j=1;
    wavelet = exp(2*1i*pi*frequency(j).*time) .* exp(-time.^2./(2*(numcycles(i)/(2*pi*frequency(j)))^2));
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
%     plot(hz_wav,test)
%     set(gca,'xlim',[40 100])
    numcycles_arry(i,1) = numcycles(i);
    FWHM_arry(i,1) = refFWHM;
end

for i=1:length(numcycles)
    refFWHM = numcycles_arry(i,1);
    for j=2:length(frequency)
        test_cycle_arry = zeros(1,50);
        for k=1:length(test_cycle_arry) % test cycles
            wavelet = exp(2*1i*pi*frequency(j).*time) .* exp(-time.^2./(2*((numcycles(i)+k-1)/(2*pi*frequency(j)))^2));
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
% 
%             plot(hz_wav,test)
%             set(gca,'xlim',[40 100])
%             pause()
        end
        [a b] =min(abs(test_cycle_arry-FWHM_arry(i,1)));
        numcycles_arry(i,j) = numcycles_arry(i,1)+b-1;
    end
end

if variablefreq==false
    for i=1:length(numcycles)
        refFWHM = numcycles_arry(i,1);
        for j=2:length(frequency)
            numcycles_arry(i,j) = refFWHM;
        end
    end
end

test = zeros([length(numcycles),length(frequency),length(linspace(0,srate/2,round(length(time)/2)+1))]);
wavlettarry = zeros([length(numcycles),length(frequency),length(time)]);
for i=1:length(numcycles)
    for j=1:length(frequency)
        % make wavelet
        wavelet = exp(2*1i*pi*frequency(j).*time) .* exp(-time.^2./(2*(numcycles_arry(i,j)/(2*pi*frequency(j)))^2));

%         subplot(2,2,i)
%         plot(time,real(wavelet),wavecolors(i))
%         xlabel('Time')
%         title([ 'Wavelet at ' num2str(frequency(j)) ' Hz with ' num2str(numcycles(i)) ' cycles' ])

        fft_wav = 2*abs(fft(wavelet));
        hz_wav  = linspace(0,srate/2,round(length(wavelet)/2)+1);
        test(i,j,:) = fft_wav(1:length(hz_wav));
        wavlettarry(i,j,:) = wavelet;
    end
end
figure
clf
hold on
for i=1:length(numcycles)
    for j=1:length(frequency)
        plot(hz_wav,squeeze(test(i,j,:)),wavecolors(1))
    end
end
   
plot(hz_wav,squeeze(mean(test(i,:,:),2)*length(frequency)/2),wavecolors(2), 'LineWidth',5)

hold off
xlabel('Frequency (Hz)')
set(gca,'xlim',[1 250])
title('Wavelet Frequency Response')

% legend({[ num2str(numcycles(1)) ' cycles' ];[ num2str(numcycles(2)) ' cycles' ]})

figure
clf
hold on
for i=1:length(numcycles)
    for j=1:length(frequency)
        tempwavlet = squeeze(abs(hilbert(real(wavlettarry(i,j,:)))));
%         plot(time,abs(real(squeeze(mean(wavlettarry(i,j,:),2)))),wavecolors(i))
        plot(time,tempwavlet,wavecolors(i))
    end
end
hold off
xlabel('Time')
set(gca,'xlim',[-.15 .15])
set(gca,'ylim',[0 1])
title('Average Wavelet Temporal Distribution')


end
