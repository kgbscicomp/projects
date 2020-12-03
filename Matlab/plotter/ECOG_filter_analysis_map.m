% analyze the frquency spectrum plots

MorletWidth=Morlet_Cycles;
input_bad_chans = bad_ch;
srate=EEG.srate;
data = EEG.data';
chan_lbls_orig=1:size(data,2);
chan_lbls = chan_lbls_orig;
base1 = round((baseline(1)/srate)-(epoch(1)/srate)+1);
base2 = round((baseline(2)*srate)-(baseline(1)*srate));

data(:,[input_bad_chans])=[]; 
chan_lbls([input_bad_chans])=[];

RejThresh = 5;
while RejThresh>0
    a=var(data);
    b=find(a>(RejThresh*median(a)));
    c=find(a<(median(a)/RejThresh));
    bad_chans=[b c];
    disp([' bad channels additional ind: ' int2str(bad_chans)])
    disp([' bad channels additional orig lbls: ' int2str(chan_lbls(bad_chans))])
    for i=1:length(bad_chans)
        disp([' bad channels names: ' (EEG.chanlocs(chan_lbls(bad_chans(i))).labels)])
    end
    figure,plot(data(1:10000,bad_chans));gcf;
    disp([' Current Threshold: ' int2str(RejThresh)])
    RejThresh= input('Enter in 0 to exit or new threshold: ');
end

% delete bad channels
removed_chan = [(input_bad_chans) chan_lbls(bad_chans)];

data(:,bad_chans)=[];
chan_lbls(bad_chans)=[];
clear a b c
disp('bad channels deleted');

chan_names = {};
for i=1:length(chan_lbls)
    temp = char(EEG.chanlocs(1,chan_lbls(i)).labels);
    chan_names{i} = temp;
end

% re-ref to common-average-reference
data=car(double(data));
disp('data in CAR')

% Morlet Filters
disp('Morlet Filters')
linefq = srate/2;
[hband_morlet] = bbpower_continuous(data, srate, filter_range, filter_width, MorletWidth, false, linefq, variablefreq);
 
% remake event (trigger) timeseries
tempnumevents = size(EEG.event,2);
stim=[];
counter = 1;
for n=1:tempnumevents
    tempevent = EEG.event(1,n).type;
    if ischar(tempevent)
        tempevent  = str2double(tempevent);
    end
    if ismember(tempevent,CondNums)
    stim(counter,:) = [EEG.event(1,n).latency+1 tempevent];
    counter=counter+1;
   end
end
numevents = size(stim,1);

% create epochs

epoch_adj = [round(epoch(1)*srate)/srate round(epoch(2)*srate)/srate];
trial_length = round((epoch_adj(2)-epoch_adj(1))*srate)+1;

disp(['Adjusted epoch time ' num2str(epoch_adj(1)*1000) ' to ' num2str(epoch_adj(2)*1000) ' ms'])

num_el = size(hband_morlet,2);

% event 92583 is at 90.4131 secs
stim_resamp = stim;
stim_resamp(:,1) = round(stim(:,1).*(srate/srate));
hband_sig_morlet = [];

%%


conds  = find(stim_resamp(:,2)==CondNums(7))';
elecs = [17:22];
plotcount = 1;
for elec = elecs
    figure
    testA = [];
    counter =1;
for yyy=conds;
    currsample = stim_resamp(yyy,1);
    temp2 = [];
    temp3 = [];
    for freqx = 1:size(hband_morlet,3)
        temp=hband_morlet(currsample+(epoch_adj(1)*srate):currsample+(epoch_adj(2)*srate),elec,freqx)';
%         temp = detrend(temp);
        base_line=(squeeze(temp(base1:base2)));
        meanx = mean(base_line);
        stdx = std(base_line);
        temp2 = (temp-meanx);
        testA(counter,freqx,:)=squeeze(mean(temp2,1));
    end 
    counter = counter+1;
end

% figure(4)
times = epoch_adj(1):1/srate:epoch_adj(2);
% plot(times,squeeze(mean(testA,1)))
% set(gca,'xlim',[-1 1.000])
testB = squeeze(mean(testA,1));
for freqx = 1:size(hband_morlet,3)    
    base_line=(testB(freqx,base1:base2));
        meanx = mean(base_line);
        stdx = std(base_line);
        testB(freqx,:) = (testB(freqx,:)-meanx)/stdx;
end
figure(2)
frex = linspace(filter_range(1),filter_range(2),size(testA,2));

subplot(1,length(elecs),plotcount)
contourf(times,frex,testB,40,'linecolor','none')
title(char(chan_names(elec)))
% ylim([-5 20])
% mesh(times,frex,freqmap')
set(gca,'xlim',[-.5 .5])
caxis([-6 6]); % color limits
% rotate3d on
hold on
plotcount = plotcount+1;
end




