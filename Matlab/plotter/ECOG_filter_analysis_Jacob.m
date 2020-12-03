% function [Induced_ButtOrd_Wide,Induced_ButtOrd_Narrow,Induced_evoked,Induced_FIRband_sig,Induced_Wide_morlet,Induced_Narrow_morlet] = ECOG_filter_analysis(subject,task,EEG,Conditions,filter_range,filter_width,epoch,baseline,srate_downsample,input_bad_chans)

% Conditions = {'2beeps','1flash','2flashes','Illusion','Ctrl','Blank','2flashes2beeps'};[Induced_ButtOrd_Wide,Induced_ButtOrd_Narrow,Induced_evoked,Induced_FIRband_sig,Induced_Wide_morlet,Induced_Narrow_morlet] = ECOG_filter_analysis('AKelly','DF',EEG,Conditions,[71 200],5,[-.1 1],[-.1 0],100,[64 74:96 123:128 16 68 97:122]);
% load ECoG data with evnts into eeglab like normal
% input into function EEG.data';

% data=eeglab2fieldtrip(EEG,'timelockanalysis','none');
cd('/Users/jacobzweig/Dropbox/ECOG_Jacob_AVCong/data');

subject = 'Rush1';
task = 'Jacob_AV';
% 
% % Conditions correspond to the event numbers in eeglab (e.g., 1-7 for DF)
Conditions = {'HighCongruent','HighIncongruent','HighAudAlone','HighVisAlone','MedCongruent','MedIncongruent','MedAudAlone','MedVisAlone','LowCongruent','LowIncongruent','LowAudAlone','LowVisAlone','DotVigilance'};
% 
% % Epoch length
epoch = [-.5 1];
baseline = [-.5 0];
% 
% downsampled rate
srate=EEG.srate;
srate_new=EEG.srate; % had this in there before in case we wanted to downsample earlier in the program
srate_downsample = 100; % rate that we'll downsample at after processing

data = EEG.data';

chan_lbls_orig=1:size(data,2);
chan_lbls = chan_lbls_orig;

filter_range = [70 200];
filter_width = 5;

% channels to remove from dataset
input_bad_chans = [41:50,81,81,105:128];

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



% High gamma band filter intervals
% hband_filters = [71 75;76 80;81 85;86 90;91 95;96 100;101 105;106 110;131 135;136 140;141 145;146 150;151 155;156 160;161 165;166 170;191 195;196 200];
hband_filtersFIR = hband_filters;% [70 200];
%hband_filters = [68 120];


% Morlet Wavelets filters

hband_morlet_wide = zeros(size(data));
for i=1:size(data,2)
    disp(['chan ' int2str(i) ' of ' int2str(size(data,2))])
    bb = bbpower(data(:,i),srate,filter_range,filter_width,30,false); %30
    hband_morlet_wide(:,i)=bb;
    clear bb
end
hband_morlet_narrow = zeros(size(data));
for i=1:size(data,2)
    disp(['chan ' int2str(i) ' of ' int2str(size(data,2))])
    bb = bbpower(data(:,i),srate,filter_range,filter_width,6,false); %30
    hband_morlet_narrow(:,i)=bb;
    clear bb
end

% Butterworth Filters

hband_buttord_wide=zeros(size(data));
hband_sig_data=zeros([size(data,1) size(hband_filters,1)]);
tempEEG = EEG;
for el=1:size(data,2)
    disp(['filter el ' int2str(el) ' of ' int2str(size(data,2))])
    for yy = 1:size(hband_filters,1)
        hband = hband_filters(yy,:);
        hband_sig_data(:,yy) = butterpass_eeglabdata(data(:,el),hband,srate,3,60);%3,60);
        hband_sig_data(:,yy) = abs(hilbert(hband_sig_data(:,yy)));
    end    
    hband_buttord_wide(:,el)=mean(hband_sig_data,2);
    hband_sig_data=zeros([size(data,1) size(hband_filters,1)]);
    
    for yy = 1:size(hband_filters,1)
        hband = hband_filters(yy,:);
        hband_sig_data(:,yy) = butterpass_eeglabdata(data(:,el),hband,srate,1,5); % winner
        hband_sig_data(:,yy) = abs(hilbert(hband_sig_data(:,yy)));
    end    
    hband_buttord_narrow(:,el)=mean(hband_sig_data,2);
    hband_sig_data=zeros([size(data,1) size(hband_filters,1)]);
    
    % fft filter as well
    
    tempEEG.data = EEG.data(chan_lbls(1,el),:);
    disp(['filter el ' int2str(el) ' of ' int2str(size(data,2))])
    for yy = 1:size(hband_filters,1)
        hband = hband_filters(yy,:);
        disp(['filter ' int2str(yy) ' of ' int2str(size(hband_filtersFIR,1))])
        firtemp = pop_eegfilt( tempEEG, hband(1), hband(2), [], [0], 0, 0, 'fir1', 0);
        hband_sig_data(:,yy) = firtemp.data';
        hband_sig_data(:,yy) = abs(hilbert(squeeze(hband_sig_data(:,yy))));
    end
    FIRband_sig(:,el) =mean(hband_sig_data,2);
    
end
clear hband_sig_data
 
% create evoked activity as well

for el=1:size(hband_buttord_wide,2)
    evoked_downsmpl(:,el)=resample(data(:,el),srate_new,srate);
end

% remake event (trigger) timeseries

numevents = size(EEG.event,2);
stim=zeros(numevents,2);
for n=1:numevents
    stim(n,:) = [EEG.event(1,n).latency+1 EEG.event(1,n).type];
end

% create epochs

epoch_adj = [round(epoch(1)*srate_new)/srate_new round(epoch(2)*srate_new)/srate_new];
trial_length = round((epoch_adj(2)-epoch_adj(1))*srate_new)+1;

disp(['Adjusted epoch time ' num2str(epoch_adj(1)*1000) ' to ' num2str(epoch_adj(2)*1000) ' ms'])

num_el = size(hband_buttord_wide,2);

% event 92583 is at 90.4131 secs
stim_resamp = stim;
stim_resamp(:,1) = round(stim(:,1).*(srate_new/srate));

% hband_buttord_wide = buttord 3 60
% hband_buttord_narrow = buttord 1 5
% evoked_downsmpl = evoked
% FIRband_sig = FIR Filter
% hband_morlet = Wavlets

for elec=1:num_el
    %disp(['filter el ' int2str(elec) ' of ' int2str(num_el)])
    for yyy=1:numevents
        currsample = stim_resamp(yyy,1);
        
        hband_sig_buttord_wide(elec,yyy,:)=...
            hband_buttord_wide(currsample+(epoch_adj(1)*srate_new):currsample+(epoch_adj(2)*srate_new),elec)';
        base_line=squeeze(hband_sig_buttord_wide(elec,yyy,1:epoch_adj(1)*srate_new*-1));
        base_line=base_line(:);
        hband_sig_buttord_wide(elec,yyy,:)=(hband_sig_buttord_wide(elec,yyy,:)-mean(base_line))/...
            std(base_line);
        
        hband_sig_buttord_narrow(elec,yyy,:)=...
            hband_buttord_narrow(currsample+(epoch_adj(1)*srate_new):currsample+(epoch_adj(2)*srate_new),elec)';
        base_line=squeeze(hband_sig_buttord_narrow(elec,yyy,1:epoch_adj(1)*srate_new*-1));
        base_line=base_line(:);
        hband_sig_buttord_narrow(elec,yyy,:)=(hband_sig_buttord_narrow(elec,yyy,:)-mean(base_line))/...
            std(base_line);
        
%         evoked_downsmpl_induced(elec,yyy,:)=...
%             evoked_downsmpl(currsample+(epoch_adj(1)*srate_new):currsample+(epoch_adj(2)*srate_new),elec)';
%         base_line=squeeze(evoked_downsmpl_induced(elec,yyy,1:epoch_adj(1)*srate_new*-1));
%         base_line=base_line(:);
%         evoked_downsmpl_induced(elec,yyy,:)=(evoked_downsmpl_induced(elec,yyy,:)-mean(base_line))/...
%             std(base_line);
%         
        evoked_downsmpl_induced(elec,yyy,:)=...
            evoked_downsmpl(currsample+(epoch_adj(1)*srate_new):currsample+(epoch_adj(2)*srate_new),elec)';
%         base_line=squeeze(evoked_downsmpl_induced(elec,yyy,1:epoch_adj(1)*srate_new*-1));
%         base_line=base_line(:);
%         evoked_downsmpl_induced(elec,yyy,:)=(evoked_downsmpl_induced(elec,yyy,:)-mean(base_line))/...
%             std(base_line);
        
        
        FIRband_sig_induced(elec,yyy,:)=...
            FIRband_sig(currsample+(epoch_adj(1)*srate_new):currsample+(epoch_adj(2)*srate_new),elec)';
        base_line=squeeze(FIRband_sig_induced(elec,yyy,1:epoch_adj(1)*srate_new*-1));
        base_line=base_line(:);
        FIRband_sig_induced(elec,yyy,:)=(FIRband_sig_induced(elec,yyy,:)-mean(base_line))/...
            std(base_line);
        
        hband_Wide_morlet_induced(elec,yyy,:)=...
            hband_morlet_wide(currsample+(epoch_adj(1)*srate_new):currsample+(epoch_adj(2)*srate_new),elec)';
        base_line=squeeze(hband_Wide_morlet_induced(elec,yyy,1:epoch_adj(1)*srate_new*-1));
        base_line=base_line(:);
        hband_Wide_morlet_induced(elec,yyy,:)=(hband_Wide_morlet_induced(elec,yyy,:)-mean(base_line))/...
            std(base_line);
        
        hband_Narrow_morlet_induced(elec,yyy,:)=...
            hband_morlet_narrow(currsample+(epoch_adj(1)*srate_new):currsample+(epoch_adj(2)*srate_new),elec)';
        base_line=squeeze(hband_Narrow_morlet_induced(elec,yyy,1:epoch_adj(1)*srate_new*-1));
        base_line=base_line(:);
        hband_Narrow_morlet_induced(elec,yyy,:)=(hband_Narrow_morlet_induced(elec,yyy,:)-mean(base_line))/...
            std(base_line);
        
    end
end


save (strcat(subject,'_2_',task));



% parse data into individual conditions

numevents_percond = zeros(length(Conditions),1);
for i=1:length(Conditions)
    numevents_percond(i) = size((find(stim_resamp==i)),1);
    Induced_ButtOrd_Wide.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
    Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
    Induced_evoked.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
    Induced_FIRband_sig.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
    Induced_Wide_morlet.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
	Induced_Narrow_morlet.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
end
for i=1:numevents
    Induced_ButtOrd_Wide.(genvarname(char(Conditions(stim_resamp(i,2))))) = Induced_ButtOrd_Wide.(genvarname(char(Conditions(stim_resamp(i,2)))))+squeeze(hband_sig_buttord_wide(:,i,:));
    Induced_ButtOrd_Narrow.(genvarname(char(Conditions(stim_resamp(i,2))))) = Induced_ButtOrd_Narrow.(genvarname(char(Conditions(stim_resamp(i,2)))))+squeeze(hband_sig_buttord_narrow(:,i,:));
    Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2))))) = Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2)))))+squeeze(evoked_downsmpl_induced(:,i,:));
    Induced_FIRband_sig.(genvarname(char(Conditions(stim_resamp(i,2))))) = Induced_FIRband_sig.(genvarname(char(Conditions(stim_resamp(i,2)))))+squeeze(FIRband_sig_induced(:,i,:));
    Induced_Wide_morlet.(genvarname(char(Conditions(stim_resamp(i,2))))) = Induced_Wide_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))+squeeze(hband_Wide_morlet_induced(:,i,:));
    Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2))))) = Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))+squeeze(hband_Narrow_morlet_induced(:,i,:));
end

%Parse to individual epochs
for i=1:length(Conditions)
    Count_Induced_ButtOrd_Wide.(genvarname(char(Conditions(i)))) = 1;
%     Count_Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i)))) = 1;
    Count_Induced_evoked.(genvarname(char(Conditions(i)))) = 1;
%     Count_Induced_FIRband_sig.(genvarname(char(Conditions(i)))) = 1;
%     Count_Induced_Wide_morlet.(genvarname(char(Conditions(i)))) = 1;
%     Count_Induced_Narrow_morlet.(genvarname(char(Conditions(i)))) = 1;
    Count_Induced_ButtOrd_WideRS.(genvarname(char(Conditions(i)))) = 1;
%     Count_Induced_ButtOrd_NarrowRS.(genvarname(char(Conditions(i)))) = 1;
%     Count_Induced_evokedRS.(genvarname(char(Conditions(i)))) = 1;
%     Count_Induced_FIRband_sigRS.(genvarname(char(Conditions(i)))) = 1;
%     Count_Induced_Wide_morletRS.(genvarname(char(Conditions(i)))) = 1;
%     Count_Induced_Narrow_morletRS.(genvarname(char(Conditions(i)))) = 1;
%     
    Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(i)))) = [];
%     Epoch_Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i)))) = [];
    Epoch_Induced_evoked.(genvarname(char(Conditions(i)))) = [];
%     Epoch_Induced_FIRband_sig.(genvarname(char(Conditions(i)))) = [];
%     Epoch_Induced_Wide_morlet.(genvarname(char(Conditions(i)))) = [];
%     Epoch_Induced_Narrow_morlet.(genvarname(char(Conditions(i)))) = [];
    Epoch_Induced_ButtOrd_WideRS.(genvarname(char(Conditions(i)))) = [];
%     Epoch_Induced_ButtOrd_NarrowRS.(genvarname(char(Conditions(i)))) = [];
%     Epoch_Induced_evokedRS.(genvarname(char(Conditions(i)))) = [];
%     Epoch_Induced_FIRband_sigRS.(genvarname(char(Conditions(i)))) = [];
%     Epoch_Induced_Wide_morletRS.(genvarname(char(Conditions(i)))) = [];
%     Epoch_Induced_Narrow_morletRS.(genvarname(char(Conditions(i)))) = [];
end
for i=1:numevents
    Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(stim_resamp(i,2)))))(:,:,Count_Induced_ButtOrd_Wide.(genvarname(char(Conditions(stim_resamp(i,2))))))...
        = squeeze(hband_sig_buttord_wide(:,i,:));
    
%     Epoch_Induced_ButtOrd_Narrow.(genvarname(char(Conditions(stim_resamp(i,2)))))(:,:,Count_Induced_ButtOrd_Narrow.(genvarname(char(Conditions(stim_resamp(i,2))))))...
%         =  squeeze(hband_sig_buttord_narrow(:,i,:));
%     
    Epoch_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2)))))(:,:,Count_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2))))))...
        = squeeze(evoked_downsmpl_induced(:,i,:));
%     
%     Epoch_Induced_FIRband_sig.(genvarname(char(Conditions(stim_resamp(i,2)))))(:,:,Count_Induced_FIRband_sig.(genvarname(char(Conditions(stim_resamp(i,2))))))...
%         = squeeze(FIRband_sig_induced(:,i,:));
%     
%     Epoch_Induced_Wide_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))(:,:,Count_Induced_Wide_morlet.(genvarname(char(Conditions(stim_resamp(i,2))))))...
%         = squeeze(hband_Wide_morlet_induced(:,i,:));
% 
%     Epoch_Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))(:,:,Count_Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2))))))...
%         = squeeze(hband_Narrow_morlet_induced(:,i,:));
%     
    %baseline and resample
    base1 = round((baseline(1)/srate_new)-(epoch(1)/srate_new)+1);
    base2 = round((baseline(2)*srate_new)-(baseline(1)*srate_new));
    for elec=1:num_el
        tempval = Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,base1:base2,Count_Induced_ButtOrd_Wide.(genvarname(char(Conditions(stim_resamp(i,2))))));
        tempmean = mean(tempval);
        Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_ButtOrd_Wide.(genvarname(char(Conditions(stim_resamp(i,2)))))) = Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_ButtOrd_Wide.(genvarname(char(Conditions(stim_resamp(i,2))))))-tempmean;
        Epoch_Induced_ButtOrd_WideRS.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_ButtOrd_WideRS.(genvarname(char(Conditions(stim_resamp(i,2)))))) = resample(Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_ButtOrd_Wide.(genvarname(char(Conditions(stim_resamp(i,2)))))),srate_downsample,srate);
        
%         tempval = Epoch_Induced_ButtOrd_Narrow.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,base1:base2,Count_Induced_ButtOrd_Narrow.(genvarname(char(Conditions(stim_resamp(i,2))))));
%         tempmean = mean(tempval);
%         Epoch_Induced_ButtOrd_Narrow.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_ButtOrd_Narrow.(genvarname(char(Conditions(stim_resamp(i,2)))))) = Epoch_Induced_ButtOrd_Narrow.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_ButtOrd_Narrow.(genvarname(char(Conditions(stim_resamp(i,2))))))-tempmean;
%         Epoch_Induced_ButtOrd_NarrowRS.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_ButtOrd_NarrowRS.(genvarname(char(Conditions(stim_resamp(i,2)))))) = resample(Epoch_Induced_ButtOrd_Narrow.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_ButtOrd_Narrow.(genvarname(char(Conditions(stim_resamp(i,2)))))),srate_downsample,srate);
%         

        tempval = Epoch_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,base1:base2,Count_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2))))));
        tempmean = mean(tempval);
        tempstd = std(tempval);
        if ZscoreOn==1
            Epoch_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2)))))) = (Epoch_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2))))))-tempmean)/tempstd;
        else
            Epoch_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2)))))) = (Epoch_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2))))))-tempmean);
        end
        Epoch_Induced_evokedRS.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2)))))) = resample(Epoch_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2)))))),srate_downsample,srate);
         
%         tempval = Epoch_Induced_FIRband_sig.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,base1:base2,Count_Induced_FIRband_sig.(genvarname(char(Conditions(stim_resamp(i,2))))));
%         tempmean = mean(tempval);
%         Epoch_Induced_FIRband_sig.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_FIRband_sig.(genvarname(char(Conditions(stim_resamp(i,2)))))) = Epoch_Induced_FIRband_sig.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_FIRband_sig.(genvarname(char(Conditions(stim_resamp(i,2))))))-tempmean;
%         Epoch_Induced_FIRband_sigRS.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_FIRband_sigRS.(genvarname(char(Conditions(stim_resamp(i,2)))))) = resample(Epoch_Induced_FIRband_sig.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_FIRband_sig.(genvarname(char(Conditions(stim_resamp(i,2)))))),srate_downsample,srate);
%           
%         tempval = Epoch_Induced_Wide_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,base1:base2,Count_Induced_Wide_morlet.(genvarname(char(Conditions(stim_resamp(i,2))))));
%         tempmean = mean(tempval);
%         Epoch_Induced_Wide_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_Wide_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))) = Epoch_Induced_Wide_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_Wide_morlet.(genvarname(char(Conditions(stim_resamp(i,2))))))-tempmean;
%         Epoch_Induced_Wide_morletRS.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_Wide_morletRS.(genvarname(char(Conditions(stim_resamp(i,2)))))) = resample(Epoch_Induced_Wide_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_Wide_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))),srate_downsample,srate);
%         
%         
%         tempval = Epoch_Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,base1:base2,Count_Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2))))));
%         tempmean = mean(tempval);
%         Epoch_Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))) = Epoch_Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2))))))-tempmean;
%         Epoch_Induced_Narrow_morletRS.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_Narrow_morletRS.(genvarname(char(Conditions(stim_resamp(i,2)))))) = resample(Epoch_Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))),srate_downsample,srate);
    end
    
    
    Count_Induced_ButtOrd_Wide.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_ButtOrd_Wide.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;
%     Count_Induced_ButtOrd_Narrow.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_ButtOrd_Narrow.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;
    Count_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_evoked.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;
%     Count_Induced_FIRband_sig.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_FIRband_sig.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;
%     Count_Induced_Wide_morlet.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_Wide_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;
%     Count_Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;
    Count_Induced_ButtOrd_WideRS.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_ButtOrd_WideRS.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;
%     Count_Induced_ButtOrd_NarrowRS.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_ButtOrd_NarrowRS.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;
%     Count_Induced_evokedRS.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_evokedRS.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;
%     Count_Induced_FIRband_sigRS.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_FIRband_sigRS.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;
%     Count_Induced_Wide_morletRS.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_Wide_morletRS.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;
%     Count_Induced_Narrow_morletRS.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_Narrow_morletRS.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;

end

%% Remove epochs with amplitudes above 4SD, currently only for buttord_wide,
%TODO: add for all conditions

% cleaning.Combined_Induced_ButtOrd_Wide =cat(3, Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(:)))))
% 
% 
% cleaning.CombinedEpochs = cat(3, Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(1)))), Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(2)))), Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(3)))), Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(4)))), Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(5)))), Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(6)))),Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(7)))),Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(8)))),Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(9)))),Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(10)))),Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(11)))),Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(12)))),Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(13)))));
% cleaning.smoosher= reshape(cleaning.CombinedEpochs,1,[]);
% cleaning.stdev = std(cleaning.smoosher);
% cleaning.removeabove = cleaning.stdev*4;
% 
% for i =1:length(Conditions)
%     cleaning.count = 1;
%     for j =1:size(Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(i)))),3)
%         cleaning.high = max(max(Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(:,:,j)));
%         cleaning.low = min(min(Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(:,:,j)));
%         if cleaning.high-cleaning.low < cleaning.removeabove
%             Cleaned_Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(:,:,cleaning.count) = Epoch_Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(:,:,j);
%             cleaning.count = cleaning.count+1;
%         end
%     end
% end
% 
% 
% RejThresh = 4;
% Induced_chan_var_array = zeros(size(hband_sig_buttord_wide,1),1);
% for currchan = 1:size(hband_sig_buttord_wide,1)
%     a = [];
%     for currevent = 1:size(hband_sig_buttord_wide,2)
%         a(currevent) = max(hband_sig_buttord_wide(currchan,currevent,:))-min(hband_sig_buttord_wide(currchan,currevent,:));
%     end
%     Induced_chan_var_array(currchan,1) = std(a);
% end
% Induced_chan_bad_array = zeros(size(hband_sig_buttord_wide,1),size(hband_sig_buttord_wide,2));
% for currchan = 1:size(hband_sig_buttord_wide,1)
%     mean_var = Induced_chan_var_array(currchan,1);
%     for currevent = 1:size(hband_sig_buttord_wide,2)
%         a = max(hband_sig_buttord_wide(currchan,currevent,:))-min(hband_sig_buttord_wide(currchan,currevent,:));
%         if a>mean_var*RejThresh
%             Induced_chan_bad_array(currchan,currevent) = 1;
%         end
%     end
% end
% temp = mean2(Induced_chan_bad_array);
% fprintf('%3.1f%% rejected \n',temp*100)
% 
% toremove = max(Induced_chan_bad_array,[],1);




%% AVERAGE and baseline
for i=1:length(Conditions)
    Induced_ButtOrd_Wide.(genvarname(char(Conditions(i)))) = [Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))/numevents_percond(i)]';
    Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i)))) = [Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i))))/numevents_percond(i)]';
    Induced_evoked.(genvarname(char(Conditions(i)))) = [Induced_evoked.(genvarname(char(Conditions(i))))/numevents_percond(i)]';
    Induced_FIRband_sig.(genvarname(char(Conditions(i)))) = [Induced_FIRband_sig.(genvarname(char(Conditions(i))))/numevents_percond(i)]';
    Induced_Wide_morlet.(genvarname(char(Conditions(i)))) = [Induced_Wide_morlet.(genvarname(char(Conditions(i))))/numevents_percond(i)]';
    Induced_Narrow_morlet.(genvarname(char(Conditions(i)))) = [Induced_Narrow_morlet.(genvarname(char(Conditions(i))))/numevents_percond(i)]';
    
    % Baseline all channels to 0
    base1 = round((baseline(1)/srate_new)-(epoch(1)/srate_new)+1);
    base2 = round((baseline(2)*srate_new)-(baseline(1)*srate_new));
    for elec=1:num_el
        tempval = Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(base1:base2,elec);
        tempmean = mean(tempval);
        Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(:,elec) = Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(:,elec)-tempmean;
        Induced_ButtOrd_WideRS.(genvarname(char(Conditions(i))))(:,elec) = resample(Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(:,elec),srate_downsample,srate);
        
        tempval = Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i))))(base1:base2,elec);
        tempmean = mean(tempval);
        Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i))))(:,elec) = Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i))))(:,elec)-tempmean;
        Induced_ButtOrd_NarrowRS.(genvarname(char(Conditions(i))))(:,elec) = resample(Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i))))(:,elec),srate_downsample,srate);
        
        tempval = Induced_evoked.(genvarname(char(Conditions(i))))(base1:base2,elec);
        tempmean = mean(tempval);
        tempstd = std(tempval);
        Induced_evoked.(genvarname(char(Conditions(i))))(:,elec) = (Induced_evoked.(genvarname(char(Conditions(i))))(:,elec)-tempmean)/tempstd;
        Induced_evokedRS.(genvarname(char(Conditions(i))))(:,elec) = resample(Induced_evoked.(genvarname(char(Conditions(i))))(:,elec),srate_downsample,srate);
        
        tempval = Induced_FIRband_sig.(genvarname(char(Conditions(i))))(base1:base2,elec);
        tempmean = mean(tempval);
        Induced_FIRband_sig.(genvarname(char(Conditions(i))))(:,elec) = Induced_FIRband_sig.(genvarname(char(Conditions(i))))(:,elec)-tempmean;
        Induced_FIRband_sigRS.(genvarname(char(Conditions(i))))(:,elec) = resample(Induced_FIRband_sig.(genvarname(char(Conditions(i))))(:,elec),srate_downsample,srate);
        
        tempval = Induced_Wide_morlet.(genvarname(char(Conditions(i))))(base1:base2,elec);
        tempmean = mean(tempval);
        Induced_Wide_morlet.(genvarname(char(Conditions(i))))(:,elec) = Induced_Wide_morlet.(genvarname(char(Conditions(i))))(:,elec)-tempmean;
        Induced_Wide_morletRS.(genvarname(char(Conditions(i))))(:,elec) = resample(Induced_Wide_morlet.(genvarname(char(Conditions(i))))(:,elec),srate_downsample,srate);
        
        tempval = Induced_Narrow_morlet.(genvarname(char(Conditions(i))))(base1:base2,elec);
        tempmean = mean(tempval);
        Induced_Narrow_morlet.(genvarname(char(Conditions(i))))(:,elec) = Induced_Narrow_morlet.(genvarname(char(Conditions(i))))(:,elec)-tempmean;
        Induced_Narrow_morletRS.(genvarname(char(Conditions(i))))(:,elec) = resample(Induced_Narrow_morlet.(genvarname(char(Conditions(i))))(:,elec),srate_downsample,srate);
    end
end


% channel quick stats by condition
for ii = 1:num_el
    chanx = [ii];
    C1 = [];
    C2 = [];
    
    tempcount=0;
    for i=1:numevents
        if stim_resamp(i,2)==1
            C1(:,size(C1,2)+1) = squeeze(hband_sig_buttord_wide(chanx,i,:));
        end
        if stim_resamp(i,2)==2
            C2(:,size(C2,2)+1) = squeeze(hband_sig_buttord_wide(chanx,i,:));
        end
        %         if stim_resamp(i,2)==3
        %             C3(:,size(C3,2)+1) = squeeze(hband_sig_buttord_wide(chanx,i,:));
        %         end
        %         if stim_resamp(i,2)==4
        %             C4(:,size(C4,2)+1) = squeeze(hband_sig_buttord_wide(chanx,i,:));
        %         end
    end
    % figure(10)
    % plot(beepindiv)
    % axis([1 308 -15 50])
    C1_err = std(C1')/sqrt(size(C1,2));
    C2_err = std(C2')/sqrt(size(C2,2));
    %     C3_err = std(C3')/sqrt(size(C3,2));
    %     C4_err = std(C4')/sqrt(size(C4,2));
    
   
    figure(3)
    hold on
    shadedErrorBar(1:551,mean(C1,2)',C1_err,'r',.5);
    shadedErrorBar(1:551,mean(C2,2)',C2_err,'b',.5);
    %     shadedErrorBar(1:564,mean(C3,2)',C3_err,'b',.5);
    %     shadedErrorBar(1:564,mean(C4,2)',C4_err,'g',.5);
    %     xlim([1 308])
    
%     figure(chanx+20)
%     temp=[Induced_ButtOrd_Wide.HighCongruent(:,chanx),Induced_ButtOrd_Wide.HighIncongruent(:,chanx)];
%     plot(temp)
    %     xlim([1 308])
    fprintf('electrode %d',chan_lbls(ii));
    fprintf('\n');
    pause;
    clf;
    hold off;
    
end






%% %Plot time frequency 
% tempHC = mean(Epoch_Induced_evoked.HighCongruent,1); %take the mean of all electrode activity
% tempHInc = mean(Epoch_Induced_evoked.HighIncongruent,1); %mean across all electrodes
% tempHC = mean(Induced_ButtOrd_Wide.HighCongruent,2);
% tempHInc = mean(Induced_ButtOrd_Wide.HighIncongruent,2);
% cycles =0;

% newtimef({tempHC tempHInc}, length(tempHC), [EEG.xmin EEG.xmax]*1000,EEG.srate,cycles, 'maxfreq',250,'plotitc','off','nfreqs',50); %,'baseline',0

% %Stats: identifying significant electrodes
% % DATA = {[Epoch_Induced_ButtOrd_Wide.HighCongruent] , [Epoch_Induced_ButtOrd_Wide.HighIncongruent] , [Epoch_Induced_ButtOrd_Wide.MedCongruent] , [Epoch_Induced_ButtOrd_Wide.MedIncongruent] , [Epoch_Induced_ButtOrd_Wide.LowCongruent] , [Epoch_Induced_ButtOrd_Wide.LowIncongruent]};
% DATA = {[Induced_ButtOrd_Wide.HighCongruent'] , [Induced_ButtOrd_Wide.HighAudAlone']}; 
% [t df pvals surog] = statcond(DATA, 'mode', 'perm', 'naccu', 1000);
% % [F df pvals] = statcond(a);

%% Treebagging
matlabpool
Treebagger = [];
paroptions = statset('UseParallel',true);

for i = 1:num_el 
    fprintf('bagging electrode %d of %d',i,num_el);fprintf('\n');
    for j = 1:size(Epoch_Induced_ButtOrd_WideRS.HighCongruent,3)
        Treebagger(i).HC(j,:) = Epoch_Induced_ButtOrd_WideRS.HighCongruent(i,:,j);
    end
    for j = 1:size(Epoch_Induced_ButtOrd_WideRS.HighIncongruent,3)
        Treebagger(i).HI(j,:) = Epoch_Induced_ButtOrd_WideRS.HighIncongruent(i,:,j);
    end
    

    Treebagger(i).HCandHI = [Treebagger(i).HI ;Treebagger(i).HC];
    Treebagger(i).HCandHI(1:20,112) = 0;
    Treebagger(i).HCandHI(21:40,112) = 1;
    Treebagger(i).HCandHI = Treebagger(i).HCandHI(randperm(size(Treebagger(i).HCandHI,1)),:); %random shuffle rows, don't think necessary
    
    model1 = TreeBagger(300,Treebagger(i).HCandHI(1:40,1:111),Treebagger(i).HCandHI(1:40,112),...
        'oobpred','on','oobvarimp','on', 'Options',paroptions );

    [typeOOB, typeOOBScores] = oobPredict(model1);
%     [conf,classorder] = confusionmat(Treebagger(i).HCandHI(1:40,112),str2num(cell2mat(typeOOB)));
    
    TotalOOBerror(i)=oobError(model1,'mode','ensemble');%oobMeanMargin
end

[sortedTotalOOBerror sortedOOBerrorindex] = sort(TotalOOBerror);
topelectrodes = chan_lbls(sortedOOBerrorindex(1:10));


%put together all electrodes for treebagging

% TopTree.HC = zeros(size(Epoch_Induced_ButtOrd_NarrowRS.HighCongruent,3),size(Epoch_Induced_ButtOrd_NarrowRS.HighCongruent,2)*size(topelectrodes,2));
% TopTree.HI = zeros(size(Epoch_Induced_ButtOrd_NarrowRS.HighIncongruent,3),size(Epoch_Induced_ButtOrd_NarrowRS.HighIncongruent,2)*size(topelectrodes,2));
TopTree.HC = []; TopTree.HI = [];
for j = 1:size(Epoch_Induced_ButtOrd_Wide.HighCongruent,3)
    temphc= [];
    temphi = [];
    for i = 1:num_el %[sortedOOBerrorindex(1:10)]
        temphc = [temphc Epoch_Induced_ButtOrd_Wide.HighCongruent(i,:,j)];
    end
    for  i =1:num_el %[sortedOOBerrorindex(1:10)]
        temphi = [temphi Epoch_Induced_ButtOrd_Wide.HighIncongruent(i,:,j)];
    end
    TopTree.HC(j,:) = temphc;
    TopTree.HI(j,:) = temphi;
end
TopTree.HCandHI = [TopTree.HC; TopTree.HI];
TopTree.size = length(TopTree.HCandHI);
TopTree.HCandHI(1:20,TopTree.size+1) = 0;
TopTree.HCandHI(21:40,TopTree.size+1) = 1;

TopTree.model = TreeBagger(300,TopTree.HCandHI(1:40,1:TopTree.size),TopTree.HCandHI(1:40,TopTree.size+1),...
    'oobpred','on','oobvarimp','on', 'Options',paroptions );

[TopTree.sortedDeltaCritDecisionSplit, TopTree.sortedVars] = sort(TopTree.model.DeltaCritDecisionSplit, 'descend');


[TopTree.typeOOB, TopTree.typeOOBScores] = oobPredict(TopTree.model);
%     [conf,classorder] = confusionmat(Treebagger(i).HCandHI(1:40,112),str2num(cell2mat(typeOOB)));
TopTree.TotalOOBerror=oobError(TopTree.model,'mode','ensemble');%oobMeanMargin


matlabpool close


%% FDR window selection and ttests
%Shuffle data to use for FDR detection of window length
%make sure comparison has equal sample rates
ttests.condstocompare = [1 3];
ttests.FDRalpha = .05;

ttests.SRate = floor(1.0/(sum(abs(epoch))/ size(Epoch_Induced_evoked.(genvarname(char(Conditions(ttests.condstocompare(1))))),2)));
ttests.CombinedEpochs = cat(3, Epoch_Induced_evoked.(genvarname(char(Conditions(ttests.condstocompare(1))))), Epoch_Induced_evoked.(genvarname(char(Conditions(ttests.condstocompare(2))))));
ttests.ShuffledEpochs = ttests.CombinedEpochs(:,:,randperm(size(ttests.CombinedEpochs,3))); %shuffle along trial(epoch) dimension
ttests.ShuffA = ttests.ShuffledEpochs(:,(abs(epoch(1))*ttests.SRate)+1:end,1:(size(ttests.ShuffledEpochs,3)/2)); %removes baseline period
ttests.ShuffB = ttests.ShuffledEpochs(:,(abs(epoch(1))*ttests.SRate)+1:end,(size(ttests.ShuffledEpochs,3)/2+1):end); %removes baseline period
% ttests.ShuffA = ttests.ShuffledEpochs(:,(abs(epoch(1))*ttests.SRate)+1:end,1:20); %removes baseline period
% ttests.ShuffB = ttests.ShuffledEpochs(:,(abs(epoch(1))*ttests.SRate)+1:end,21:40); %removes baseline period
ttests.MeanShuffA= mean(ttests.ShuffA,3);
ttests.MeanShuffB= mean(ttests.ShuffB,3);

%Minimum number of consecutive ms required to yield alpha*numelec of electrodes as
%sig between two conditions
% ttests.WindowLength = size(ttests.MeanShuffA,2);
ttests.WindowLength=1;
ttests.WindowDivider=size(ttests.MeanShuffA,2)/ttests.WindowLength;
FDRbreak = 0;
while FDRbreak == 0
    ttests.FDRh = []; ttests.FDRp=[];
    for i =1:num_el
        for j = 1:floor(size(ttests.MeanShuffA,2)/ttests.WindowLength)
             temp(1,j) = mean(ttests.MeanShuffA(i,j*ttests.WindowLength-ttests.WindowLength+1:j*ttests.WindowLength));
             temp(2,j) = mean(ttests.MeanShuffB(i,j*ttests.WindowLength-ttests.WindowLcength+1:j*ttests.WindowLength));
        end
        [ttests.FDRh(i,1),ttests.FDRp(i,1)] = ttest(temp(1,:),temp(2,:),'alpha',.01);
    end
    
    %now test if less than 5% of elecs are significant (4 elecs?)
    if sum(ttests.FDRh) <= floor(ttests.FDRalpha*num_el)
        FDRbreak =1;
    elseif sum(ttests.FDRh) > floor(ttests.FDRalpha*num_el)
        dividebreak=0;
        while dividebreak==0
            ttests.WindowDivider = ttests.WindowDivider-1;
            ttests.WindowLengthGen = floor(size(ttests.MeanShuffA,2)/ttests.WindowDivider);
            if ttests.WindowLengthGen > ttests.WindowLength
                dividebreak=1;
                ttests.WindowLength = ttests.WindowLengthGen;
            end
        end
        fprintf('Window length %d \n', ttests.WindowLength);
    end
    if ttests.WindowLength ==size(ttests.MeanShuffA,2);
        fprintf('No Appropriate Window Length Found')
        FDRbreak=1;
    end
end


for i =1:num_el
    temp=[];
    for j = 1:floor((length(Induced_evoked.(genvarname(char(Conditions(ttests.condstocompare(1))))))-((abs(epoch(1))*ttests.SRate)+1))/ttests.WindowLength)
        temp(1,j) = mean(Epoch_Induced_evoked.(genvarname(char(Conditions(ttests.condstocompare(1)))))(i,((abs(epoch(1))*ttests.SRate)+1)+j*ttests.WindowLength-ttests.WindowLength+1:((abs(epoch(1))*ttests.SRate)+1)+j*ttests.WindowLength));
        temp(2,j) = mean(Epoch_Induced_evoked.(genvarname(char(Conditions(ttests.condstocompare(2)))))(i,((abs(epoch(1))*ttests.SRate)+1)+j*ttests.WindowLength-ttests.WindowLength+1:((abs(epoch(1))*ttests.SRate)+1)+j*ttests.WindowLength));
    end
    [ttests.h(i,1),ttests.p(i,1)] = ttest(temp(1,:),temp(2,:),'alpha',.001);

end


        

%%  save all data at the end
clear ALLCOM ALLEEG ALLERP ALLERPCOM CURRENTERP CURRENTSET CURRENTSTUDY EEG ERP LASTCOM PLUGINLIST STUDY eeglabUpdater plotset
save (strcat(subject,'_2_',task));

% end





