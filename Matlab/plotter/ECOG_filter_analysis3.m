function ECOG_filter_analysis3(subject,task,ALLEEG, EEG,CURRENTSET, Conditions,CondNums,filter_range,filter_width,epoch,baseline,srate_downsample,input_bad_chans,ButterRp,ButterRs, ZscoreOn)

% Conditions = {'2beeps','1flash','2flashes','Illusion','Ctrl','Blank','2flashes2beeps'};[Induced_ButtOrd_Wide,Induced_ButtOrd_Narrow,Induced_evoked,Induced_FIRband_sig,Induced_Wide_morlet,Induced_Narrow_morlet] = ECOG_filter_analysis('AKelly','DF',EEG,Conditions,[71 200],5,[-.1 1],[-.1 0],100,[64 74:96 123:128 16 68 97:122]);
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
srate_new=EEG.srate; % had this in there before in case we wanted to downsample earlier in the program
% srate_downsample = 100; % rate that we'll downsample at after processing

data = EEG.data';
ERPband = [.1 40];

chan_lbls_orig=1:size(data,2);
chan_lbls = chan_lbls_orig;

% filter_range = [71 200];
% filter_width = 5;

base1 = round((baseline(1)/srate_new)-(epoch(1)/srate_new)+1);
base2 = round((baseline(2)*srate_new)-(baseline(1)*srate_new));

% channels to remove from dataset
% input_bad_chans = [64 74:96 123:128 16 68 97:122];

% delete input_bad_chans and eeg_chans

data(:,[input_bad_chans])=[]; 
chan_lbls([input_bad_chans])=[];

% cut bad chans -- bad channel detector here (could also be done manually) - this just finds excessive outliers in terms of power
% check data-- checks 5 times variance and then visually inspect each marked channel. 

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

% Create Filter ranges

hband_filters_temp = [filter_range(1):filter_width:filter_range(2);filter_range(1)+filter_width-1:filter_width:filter_range(2)+filter_width-1]';
hband_filters = [];
badfreqs = [];
for i=51:69
    badfreqs(1,size(badfreqs,2)+1:size(badfreqs,2)+size(i:60:2000,2)) = i:60:2000;
end
for i=1:size(hband_filters_temp,1)
    if sum(ismember(hband_filters_temp(i,1):hband_filters_temp(i,2),badfreqs))>0
    else
        hband_filters(size(hband_filters,1)+1,:) = hband_filters_temp(i,:);
    end
end

% Butterworth Filters

hband_buttord_wide=zeros(size(data));
hband_ERP_data=zeros(size(data));
hband_sig_data=zeros([size(data,1) size(hband_filters,1)]);
tempEEG = EEG;
for el=1:size(data,2)
    disp(['filter el ' int2str(el) ' of ' int2str(size(data,2))])
    for yy = 1:size(hband_filters,1)
        hband = hband_filters(yy,:);
        hband_sig_data(:,yy) = butterpass_eeglabdata(data(:,el),hband,srate,ButterRp,ButterRs);%3,60);
        hband_sig_data(:,yy) = abs(hilbert(hband_sig_data(:,yy)));
    end    
    hband_buttord_wide(:,el)=mean(hband_sig_data,2);
    hband_sig_data=zeros([size(data,1) size(hband_filters,1)]);

    % create evoked activity as well
    hband_ERP_data(:,el) = butterpass_eeglabdata(data(:,el),ERPband,srate,ButterRp,ButterRs);%3,60);
    evoked_downsmpl(:,el)=resample(data(:,el),srate_new,srate);
end
clear hband_sig_data
 
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

epoch_adj = [round(epoch(1)*srate_new)/srate_new round(epoch(2)*srate_new)/srate_new];
trial_length = round((epoch_adj(2)-epoch_adj(1))*srate_new)+1;

disp(['Adjusted epoch time ' num2str(epoch_adj(1)*1000) ' to ' num2str(epoch_adj(2)*1000) ' ms'])

num_el = size(hband_buttord_wide,2);

% event 92583 is at 90.4131 secs
stim_resamp = stim;
stim_resamp(:,1) = round(stim(:,1).*(srate_new/srate));

for elec=1:num_el
    %disp(['filter el ' int2str(elec) ' of ' int2str(num_el)])
    for yyy=1:numevents
        currsample = stim_resamp(yyy,1);
        
        hband_sig_buttord_wide(elec,yyy,:)=...
            hband_buttord_wide(currsample+(epoch_adj(1)*srate_new):currsample+(epoch_adj(2)*srate_new),elec)';
        base_line=squeeze(hband_sig_buttord_wide(elec,yyy,base1:base2));
        base_line=base_line(:);
        
        if ZscoreOn==1
            hband_sig_buttord_wide(elec,yyy,:)=(hband_sig_buttord_wide(elec,yyy,:)-mean(base_line))/...
                std(base_line);
        else
            hband_sig_buttord_wide(elec,yyy,:)=(hband_sig_buttord_wide(elec,yyy,:)-mean(base_line));
        end
        evoked_downsmpl_induced(elec,yyy,:)=...
            evoked_downsmpl(currsample+(epoch_adj(1)*srate_new):currsample+(epoch_adj(2)*srate_new),elec)';
    end
end

% parse data into individual conditions
disp(['Parsing data into individual conditions'])

numevents_percond = zeros(length(Conditions),1);
for i=1:length(Conditions)
    numevents_percond(i) = size((find(stim_resamp==CondNums(i))),1);
    Induced_ButtOrd_Wide.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
    Induced_evoked.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
end
for i=1:numevents
    Induced_ButtOrd_Wide.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums))))) = Induced_ButtOrd_Wide.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums)))))+squeeze(hband_sig_buttord_wide(:,i,:));
    Induced_evoked.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums))))) = Induced_evoked.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums)))))+squeeze(evoked_downsmpl_induced(:,i,:));
end

disp(['Baseline all channels to 0 and Z-norming to baseline period'])

for i=1:length(Conditions)
    Induced_ButtOrd_Wide.(genvarname(char(Conditions(i)))) = [Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))/numevents_percond(i)]';
    Induced_evoked.(genvarname(char(Conditions(i)))) = [Induced_evoked.(genvarname(char(Conditions(i))))/numevents_percond(i)]';
    
    % Baseline all channels to 0
    
    for elec=1:num_el
        tempval = Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(base1:base2,elec);
        tempmean = mean(tempval);
        Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(:,elec) = Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(:,elec)-tempmean;
        Induced_ButtOrd_WideRS.(genvarname(char(Conditions(i))))(:,elec) = resample(Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(:,elec),srate_downsample,srate);
      
        tempval = Induced_evoked.(genvarname(char(Conditions(i))))(base1:base2,elec);
        tempmean = mean(tempval);
        tempstd = std(tempval);
        if ZscoreOn==1
            Induced_evoked.(genvarname(char(Conditions(i))))(:,elec) = (Induced_evoked.(genvarname(char(Conditions(i))))(:,elec)-tempmean)/tempstd;
        else
            Induced_evoked.(genvarname(char(Conditions(i))))(:,elec) = (Induced_evoked.(genvarname(char(Conditions(i))))(:,elec)-tempmean);
        end
        Induced_evokedRS.(genvarname(char(Conditions(i))))(:,elec) = resample(Induced_evoked.(genvarname(char(Conditions(i))))(:,elec),srate_downsample,srate);
        
        
    end
    disp(['Norming Condition ' int2str(i) ' of ' int2str(length(Conditions))])

end


% save all data at the end
save(strcat(subject,'_',task),'-v7.3', '-regexp', '^(?!(ALLCOM|ALLEEG|ALLERP|ALLERPCOM|CURRENTERP|CURRENTSET|CURRENTSTUDY|EEG|ERP|LASTCOM|PLUGINLIST|STUDY|eeglabUpdater|plotset|data|tempEEG)$).');

end





