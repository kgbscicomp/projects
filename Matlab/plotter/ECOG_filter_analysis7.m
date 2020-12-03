function ECOG_filter_analysis7(subject,task,ALLEEG, EEG,CURRENTSET, Conditions,CondNums,input_bad_chans,epoch)

% Conditions = {'2beeps','1flash','2flashes','Illusion','Ctrl','Blank','2flashes2beeps'};[Induced_Morlett,Induced_ButtOrd_Narrow,Induced_evoked,Induced_FIRband_sig,Induced_Wide_morlet,Induced_Narrow_morlet] = ECOG_filter_analysis('AKelly','DF',EEG,Conditions,[71 200],5,[-.1 1],[-.1 0],100,[64 74:96 123:128 16 68 97:122]);
% load ECoG data with evnts into eeglab like normal
% input into function EEG.data';

% subject = 'AKelly';
% task = 'DF';
% 
% % Conditions correspond to the event numbers in eeglab (e.g., 1-7 for DF)
% Conditions = {'2beeps','1flash','2flashes','Illusion','Ctrl','Blank','2flashes2beeps'};
% CondNums = 1:7; % cond numbers to look for
% 

srate = EEG.srate;
data = EEG.data';

chan_lbls_orig=1:size(data,2);
chan_lbls = chan_lbls_orig;
if nargin < 9
    epoch = [-2 2];
end

% delete input_bad_chans and eeg_chans

data(:,[input_bad_chans])=[]; 
chan_lbls([input_bad_chans])=[];

% cut bad chans -- bad channel detector here (could also be done manually) - this just finds excessive outliers in terms of power
% check data-- checks 5 times variance and then visually inspect each marked channel. 

% save test.mat
% 
% beep
% pause();


RejThresh = 3;
RejThreshNew = RejThresh;
while RejThreshNew>0
    RejThresh = RejThreshNew;
    tempdata = data;
    plottempdata = [];
    bad_chans_temp = 1;
    bad_chans = [];
    b = [];c=[];
    while isempty(b) && isempty(c) && ~isempty(bad_chans_temp)
        tempdata = car(tempdata);
        a=var(tempdata);
        b=find(a>(RejThreshNew*median(a)));
        c=find(a<(median(a)/RejThreshNew));
        bad_chans_temp=[b c];
        plottempdata = [plottempdata tempdata(1:5000,bad_chans_temp)];
        tempdata(:,bad_chans_temp)=[];
        tempdata = car(tempdata);
        bad_chans = [bad_chans bad_chans_temp];
    end
    disp([' bad channels additional ind: ' int2str(bad_chans)])
    disp([' bad channels additional orig lbls: ' int2str(chan_lbls(bad_chans))])
    for i=1:length(bad_chans)
        disp([' bad channels names: ' (EEG.chanlocs(chan_lbls(bad_chans(i))).labels)])
    end
    
    figure(2),
    clf
    plot(plottempdata);gcf;
    disp([' Current Threshold: ' int2str(RejThreshNew)])
    RejThreshNew= input('Enter in 0 to exit or new threshold: ');
end

% delete bad channels
removed_chan = [(input_bad_chans) chan_lbls(bad_chans)];

evoked = double(EEG.data)';
evoked(:,[input_bad_chans])=[]; 

evoked(:,bad_chans)=[];
chan_lbls(bad_chans)=[];
clear a b c
disp('bad channels deleted');

chan_names = {};
for i=1:length(chan_lbls)
    temp = char(EEG.chanlocs(1,chan_lbls(i)).labels);
    chan_names{i} = temp;
end

% re-ref to common-average-reference

% evoked = laplacian(evoked);
evoked = car(evoked);
disp('data in CAR')

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
        stim(counter,:) = [EEG.event(1,n).latency+1 tempevent]; %EEG.event(1,x).latency is in samples
        counter=counter+1;
    end
end
numevents = size(stim,1);

% create epochs

trial_length = round((epoch(2)-epoch(1))*srate)+1;
num_el = size(evoked,2);
evoked_sig = [];
% save('test.mat');
for yyy=1:numevents
    disp(['trial ' int2str(yyy) ' of ' int2str(numevents)])
    currsample = round(stim(yyy,1));
    tempwindow = [(epoch(1)*srate) (epoch(2)*srate)];
    temprng = currsample+tempwindow(1):currsample+tempwindow(2);
    temp=double(evoked(temprng,:));
    temp2 = [];
    for i=1:size(temp,2)
        temp2(:,i) = temp(:,i)-mean(temp(1:1+(epoch(1)*-srate),i));
    end
    evoked_sig(:,:,yyy)=temp2';
end

% parse data into individual conditions
disp(['Parsing data into individual conditions'])

for i=1:length(Conditions)
    temp = find(stim(:,2)==CondNums(i));
    Evoked.(genvarname(char(Conditions(i)))) = squeeze(mean(evoked_sig(:,:,temp),3));
end

disp('Saving Data')

% save all data at the end
save(strcat(subject,'_',task),'-v7.3', '-regexp', '^(?!(ALLCOM|ALLEEG|ALLERP|ALLERPCOM|CURRENTERP|CURRENTSET|CURRENTSTUDY|EEG|ERP|LASTCOM|PLUGINLIST|STUDY|eeglabUpdater|plotset|data|evoked)$).');

end





