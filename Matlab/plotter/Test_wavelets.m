
srate=1024;


%% Figure 13.13
% Test Wavelet construction

srate=1024;

frequency = 200;
time      = -1:1/srate:1;
numcycles = [10 15 ];

wavecolors = 'br';

figure
for i=1:length(numcycles)
    % make wavelet
    wavelet = exp(2*1i*pi*frequency.*time) .* exp(-time.^2./(2*(numcycles(i)/(2*pi*frequency))^2));
    
    subplot(2,2,i)
    plot(time,real(wavelet),wavecolors(i))
    xlabel('Time')
    title([ 'Wavelet at ' num2str(frequency) ' Hz with ' num2str(numcycles(i)) ' cycles' ])
    set(gca,'xlim',[-.2 .2])
    
    subplot(2,1,2)
    hold on
    fft_wav = 2*abs(fft(wavelet));
    hz_wav  = linspace(0,srate/2,round(length(wavelet)/2)+1);
    plot(hz_wav,fft_wav(1:length(hz_wav)),wavecolors(i))
end

subplot(212)
xlabel('Frequency (Hz)')
set(gca,'xlim',[0 250])
legend({[ num2str(numcycles(1)) ' cycles' ];[ num2str(numcycles(2)) ' cycles' ]})

%% Figure 13.13
% Test Wavelet construction

% VARY FREQ
% frequency = [71:5:100] ;


linefq = 60; 
fqstep = 5; 
fqrng = [70 200]; 
line_param = 6; 
allfqs = (fqrng(1):fqstep:fqrng(end))';
badfqs = and(allfqs > line_param, abs(mod(allfqs,linefq) - linefq/2) > (linefq/2 - line_param));
allfqs = allfqs(~badfqs);

frequency = allfqs;
time      = -1:1/srate:1;
numcycles = [15 20];
srate=1024;

wavecolors = 'br';

test = zeros([length(numcycles),length(frequency),length(linspace(0,srate/2,round(length(time)/2)+1))]);
wavlettarry = zeros([length(numcycles),length(frequency),length(time)]);
for i=1:length(numcycles)
    for j=1:length(frequency)
        % make wavelet
        wavelet = exp(2*1i*pi*frequency(j).*time) .* exp(-time.^2./(2*(numcycles(i)/(2*pi*frequency(j)))^2));

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
figure(2)
clf
hold on
for i=1:length(numcycles)
    for j=1:length(frequency)
        plot(hz_wav,squeeze(test(i,j,:)),wavecolors(i))
    end
end
hold off
xlabel('Frequency (Hz)')
set(gca,'xlim',[60 250])
% legend({[ num2str(numcycles(1)) ' cycles' ];[ num2str(numcycles(2)) ' cycles' ]})

figure(3)
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
set(gca,'xlim',[-.2 .2])
% legend({[ num2str(numcycles(1)) ' cycles' ];[ num2str(numcycles(2)) ' cycles' ]})


% plot(time,tempwavlet)
% xlabel('Time')
% set(gca,'xlim',[-.2 .2])



figure(4)
clf
hold on
for i=1:length(numcycles)
    plot(hz_wav,squeeze(mean(test(i,:,:),2)),wavecolors(i))
end
hold off
xlabel('Frequency (Hz)')
set(gca,'xlim',[60 250])



%%
frequency = [100:10:200] ;
time      = -1:1/srate:1;
numcycles = [5 10];
srate=1024;
figure(2)
hold on
for i=1:length(numcycles)
    for j=1:length(frequency)
        % make wavelet
        
        
%         wavelet = exp(2*1i*pi*frequency(j).*time) .* exp(-time.^2./(2*(numcycles(i)/(2*pi*frequency(j)))^2));

%         subplot(2,2,i)
%         plot(time,real(wavelet),wavecolors(i))
%         xlabel('Time')
%         title([ 'Wavelet at ' num2str(frequency(j)) ' Hz with ' num2str(numcycles(i)) ' cycles' ])

        fft_wav = 2*abs(fft(wavelet));
        hz_wav  = linspace(0,srate/2,round(length(wavelet)/2)+1);
        plot(hz_wav,fft_wav(1:length(hz_wav)))
    end
end
hold off
xlabel('Frequency (Hz)')
set(gca,'xlim',[0 300])


%%


%% Figure 13.13
% Test Wavelet construction

% VARY wavelet width by frequency
% Set the wavelet width for the very first cycle then it will match the
% FWHM at all other frequencies -- 

linefq = srate/2; % changed from 60 hz line activity to include all freqs
fqstep = 5; 
fqrng = [70 200]; 
line_param = 5; 
numcycles = [15 20];

allfqs = (fqrng(1):fqstep:fqrng(end))';
badfqs = and(allfqs > line_param, abs(mod(allfqs,linefq) - linefq/2) > (linefq/2 - line_param));
allfqs = allfqs(~badfqs);

frequency = allfqs;
time      = -1:1/srate:1;
srate=1024;

wavecolors = 'brg';

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
figure(2)
clf
hold on
for i=1:length(numcycles)
    for j=1:length(frequency)
        plot(hz_wav,squeeze(test(i,j,:)),wavecolors(i))
    end
end
hold off
xlabel('Frequency (Hz)')
set(gca,'xlim',[50 250])
% legend({[ num2str(numcycles(1)) ' cycles' ];[ num2str(numcycles(2)) ' cycles' ]})

figure(3)
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
% legend({[ num2str(numcycles(1)) ' cycles' ];[ num2str(numcycles(2)) ' cycles' ]})

% plot(time,tempwavlet)
% xlabel('Time')
% set(gca,'xlim',[-.2 .2])



figure(4)
clf
hold on
for i=1:length(numcycles)
    plot(hz_wav,squeeze(mean(test(i,:,:),2)),wavecolors(i))
end
hold off
xlabel('Frequency (Hz)')
set(gca,'xlim',[50 250])

%%

% test on data

% data = squeeze(EEG.data(39,:,:));

linefq = srate/2; 
fqstep = 5; 
fqrng = [70 200]; 
line_param = 10; 
wav_width = 20;
fs=1024;
srate=1024;
do_log = false;
variablefreq = true;

% generate sine waves of known properties
% list some frequencies
frex = [80 70 150];

% list some random amplitudes... make sure there are 
% the same number of amplitudes as there are frequencies!
amplit = [5 10 2];

% phases... list some random numbers between -pi and pi
phases = [  pi/3  pi/8 pi];%  pi  pi/2  -pi/4 ];

% define time...
time1=-2:1/fs:0-1/fs;
time2=0:1/fs:2-1/fs;
times = [time1 time2];

% now we loop through frequencies and create sine waves
sine_waves = zeros(length(frex),length(time1)+length(time2)); % remember: always initialize!
for fi=1:2
    sine_waves(fi,:) = [amplit(fi) * sin(2*pi*frex(fi).*time1 + phases(fi)), 5*amplit(fi) * sin(2*pi*frex(fi).*time1 + phases(fi))];
end
sine_waves(3,:) = amplit(fi) * sin(2*pi*frex(fi).*times + phases(fi));

% bb = bbpower(sum(sine_waves)', fs, fqrng, fqstep, wav_width, do_log, linefq, variablefreq, baseline1,baseline2);
[bb freqmap] = bbpower(data, fs, fqrng, fqstep, wav_width, do_log, linefq, variablefreq, baseline1,baseline2);
% test = mean(bb(1024:3072,:),2);
test = mean(bb,2);
xmean = mean(test(1:923));
xstd = std(test(1:923));
test = (test-xmean)/xstd;

for i=1:size(freqmap,2)
    xmean = mean(freqmap(1:923,i));
    xstd = std(freqmap(1:923,i));
    freqmap(:,i) = (freqmap(:,i)-xmean)/xstd;
end

times = -1:1/srate:1;

figure(4)
plot(times,test)
ylim([-5 30])
set(gca,'xlim',[-1 1.000])
figure
frex = linspace(fqrng(1),fqrng(2),size(freqmap,2));
contourf(times,frex,freqmap',40,'linecolor','none')
% mesh(times,frex,freqmap')
set(gca,'xlim',[-1 1.000])
% rotate3d on






