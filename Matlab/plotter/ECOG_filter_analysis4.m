function ECOG_filter_analysis4(subject,task,ALLEEG, EEG,CURRENTSET, Conditions,CondNums,filter_range,filter_width,epoch,baseline,input_bad_chans,MorletWidth,variablefreq)

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
% % Epoch length
% epoch = [-.1 1];
% baseline = [-.1 0];
% 

% ZscoreOn = 1 means we'll z-transform each individual epoch before averaging,
% off and we demean each epoch before averaging instead

% downsampled rate
srate=EEG.srate;

data = EEG.data';
ERPband = [.1 40];


chan_lbls_orig=1:size(data,2);
chan_lbls = chan_lbls_orig;

% filter_range = [71 200];
% filter_width = 5;

base1 = round((baseline(1)/srate)-(epoch(1)/srate)+1);
base2 = round((baseline(2)*srate)-(baseline(1)*srate));

% channels to remove from dataset
% input_bad_chans = [64 74:96 123:128 16 68 97:122];

% delete input_bad_chans and eeg_chans

data(:,[input_bad_chans])=[]; 
chan_lbls([input_bad_chans])=[];

% cut bad chans -- bad channel detector here (could also be done manually) - this just finds excessive outliers in terms of power
% check data-- checks 5 times variance and then visually inspect each marked channel. 

% save test.mat

RejThresh = 3;
while RejThresh>0
    tempdata = data;
    plottempdata = [];
    bad_chans_temp = 1;
    bad_chans = [];
    b = [];c=[];
    while isempty(b) && isempty(c) && ~isempty(bad_chans_temp)
        tempdata = car(tempdata);
        a=var(tempdata);
        b=find(a>(RejThresh*median(a)));
        c=find(a<(median(a)/RejThresh));
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
    disp([' Current Threshold: ' int2str(RejThresh)])
    RejThresh= input('Enter in 0 to exit or new threshold: ');
end

% delete bad channels
removed_chan = [(input_bad_chans) chan_lbls(bad_chans)];

disp('Low-Pass Filtering Evoked')
% % % % % EEG = pop_eegfilt( EEG, 0, 40, [], [0], 0, 0, 'fir1', 0);


% EEG = pop_ERPLAB_polydetrend( EEG, 5120, 5);


evoked = double(EEG.data)';
evoked(:,[input_bad_chans])=[]; 

data(:,bad_chans)=[];
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

data=car(double(data));
evoked = car(evoked);
disp('data in CAR')

% Create Filter ranges

save('test.mat');
pause();

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
evoked_sig = [];

for elec=1:num_el
    disp(['filter el ' int2str(elec) ' of ' int2str(num_el)])
    for yyy=1:numevents
        currsample = stim_resamp(yyy,1);
%         temp2 = [];
%         for freqx = 1:size(hband_morlet,3)
%             temp=hband_morlet(currsample+(epoch_adj(1)*srate):currsample+(epoch_adj(2)*srate),elec,freqx)';
        temp=hband_morlet(currsample+(epoch_adj(1)*srate):currsample+(epoch_adj(2)*srate),elec,:);
%             temp = temp - mean(temp);
%             base_line=squeeze(temp(base1:base2));
%             meanx = mean(base_line);
%             stdx = std(base_line);
%             temp2(freqx,:) = (temp-meanx)/5;
%         end
%         hband_sig_morlet(elec,yyy,:)=squeeze(mean(temp2,1));
        temp2 = squeeze(mean(temp,3))';
        base_line=squeeze(temp2(base1:base2));
        meanx = mean(base_line);
        stdx = std(base_line);
        temp2 = (temp2-meanx);
        hband_sig_morlet(elec,yyy,:)=temp2;
        evoked_sig(elec,yyy,:)=...
            evoked(currsample+(epoch_adj(1)*srate):currsample+(epoch_adj(2)*srate),elec)';
    end
end
% 
% test = squeeze(hband_sig_morlet(39,:,:));
% figure(3)
% plot(mean(test,1)');figure(gcf);

% parse data into individual conditions
disp(['Parsing data into individual conditions'])

% numevents_percond = zeros(length(Conditions),1);
% for i=1:length(Conditions)
%     numevents_percond(i) = size((find(stim_resamp==CondNums(i))),1);
%     Induced_Morlett.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
%     Induced_evoked.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
% end
% for i=1:numevents
%     Induced_Morlett.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums))))) = Induced_Morlett.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums)))))+squeeze(hband_sig_morlet(:,i,:));
%     Induced_evoked.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums))))) = Induced_evoked.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums)))))+squeeze(evoked_sig(:,i,:));
% end
for i=1:length(Conditions)
    temp = find(stim_resamp(:,2)==CondNums(i));
    Induced_Morlett.(genvarname(char(Conditions(i)))) = squeeze(mean(hband_sig_morlet(:,temp,:),2));
    Induced_evoked.(genvarname(char(Conditions(i)))) = squeeze(mean(evoked_sig(:,temp,:),2));
end


disp(['Baseline all channels to 0 and Z-norming to baseline period'])

for i=1:length(Conditions)
    Induced_Morlett.(genvarname(char(Conditions(i)))) = Induced_Morlett.(genvarname(char(Conditions(i))))';
    Induced_evoked.(genvarname(char(Conditions(i)))) = Induced_evoked.(genvarname(char(Conditions(i))))';
%     
    % Baseline all channels to 0
    
    for elec=1:num_el
        tempval = Induced_Morlett.(genvarname(char(Conditions(i))))(base1:base2,elec);
        tempmean = mean(tempval);
        tempstd = std(tempval);
        Induced_Morlett.(genvarname(char(Conditions(i))))(:,elec) = (Induced_Morlett.(genvarname(char(Conditions(i))))(:,elec)-tempmean)/tempstd;
        
        tempval = Induced_evoked.(genvarname(char(Conditions(i))))(base1:base2,elec);
        tempmean = mean(tempval);
        tempstd = std(tempval);
        Induced_evoked.(genvarname(char(Conditions(i))))(:,elec) = (Induced_evoked.(genvarname(char(Conditions(i))))(:,elec)-tempmean)/tempstd;
        
    end
    disp(['Norming Condition ' int2str(i) ' of ' int2str(length(Conditions))])

end


% save all data at the end
save(strcat(subject,'_',task),'-v7.3', '-regexp', '^(?!(ALLCOM|ALLEEG|ALLERP|ALLERPCOM|CURRENTERP|CURRENTSET|CURRENTSTUDY|EEG|ERP|LASTCOM|PLUGINLIST|STUDY|eeglabUpdater|plotset|data|tempEEG|hband_morlet|evoked)$).');

end





