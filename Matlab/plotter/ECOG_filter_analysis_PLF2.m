function ECOG_filter_analysis_PLF2(subject,task,ALLEEG, EEG,CURRENTSET, Conditions,CondNums,filter_range,epoch,baseline,input_bad_chans,ButterRp,ButterRs)

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

% EEG  = pop_basicfilter( EEG,  1 , 'Cutoff',  60, 'Design', 'notch', 'Filter', 'PMnotch', 'Order',  180 );
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
% EEG = eeg_checkset( EEG );


% save (strcat(subject,'_',task));
% pause();

% Create Filter ranges
hband_filters = filter_range;

% High gamma band filter intervals
% hband_filters = [71 75;76 80;81 85;86 90;91 95;96 100;101 105;106 110;131 135;136 140;141 145;146 150;151 155;156 160;161 165;166 170;191 195;196 200];

% Butterworth Filters
hband_buttord_wide=zeros(size(data,1),size(hband_filters,1),size(data,2));
hband_sig_data=zeros([size(data,1) size(hband_filters,1)]);
tempEEG = EEG;
for el=1:size(data,2)
    disp(['filter el ' int2str(el) ' of ' int2str(size(data,2))])
    for yy = 1:size(hband_filters,1)
        hband = hband_filters(yy,:);
        hband_sig_data(:,yy) = butterpass_eeglabdata(data(:,el),hband,srate,ButterRp,ButterRs);%3,60);
        hband_sig_data(:,yy) = (hilbert(hband_sig_data(:,yy)));
        hband_buttord_wide(:,yy,el)=mean(hband_sig_data,2);
    end    
    
    hband_sig_data=zeros([size(data,1) size(hband_filters,1)]);
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

num_el = size(hband_buttord_wide,3);

% event 92583 is at 90.4131 secs
stim_resamp = stim;
stim_resamp(:,1) = round(stim(:,1).*(srate_new/srate));

% hband_buttord_wide = buttord 3 60
% hband_buttord_narrow = buttord 1 5
% evoked_downsmpl = evoked
% FIRband_sig = FIR Filter
% hband_morlet = Wavlets

clear hband_sig_buttord_wide
for elec=1:num_el
    disp(['filter el ' int2str(elec) ' of ' int2str(num_el)])
    for xxx=1:size(hband_filters,1)
        for yyy=1:numevents
            currsample = stim_resamp(yyy,1);
            hband_sig_buttord_wide(yyy,:,xxx,elec)=...
                squeeze(hband_buttord_wide(currsample+(epoch_adj(1)*srate_new):currsample+(epoch_adj(2)*srate_new),xxx,elec))';
        end
    end
end


% save (strcat(subject,'_',task));

% parse data into individual conditions
% % 
% % numevents_percond = zeros(length(Conditions),1);
% % for i=1:length(Conditions)
% %     numevents_percond(i) = size((find(stim_resamp==CondNums(i))),1);
% %     Induced_ButtOrd_Wide.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
% % end
% % % for i=[1:7 9:17]
% % for i=[11 14 15]
% %     CurrCondTrials = find(stim_resamp(:,2)==CondNums(i));
% %     for j=1:trial_length
% %         for k=1:num_el
% %             temp = squeeze(hband_sig_buttord_wide(k,CurrCondTrials',j));
% %             [pval,zval] = circ_rtest(temp);
% %             Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(k,j) = zval;
% %         end
% %     end
% %     i
% % end

% save all data at the end
clear ALLCOM ALLEEG ALLERP ALLERPCOM CURRENTERP CURRENTSET CURRENTSTUDY EEG ERP LASTCOM PLUGINLIST STUDY eeglabUpdater plotset
clear data tempEEG 
save(strcat(subject,'_',task),'-v7.3');

end





