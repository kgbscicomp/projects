function ECOG_filter_analysis_Jacob2(subject,task,ALLEEG, EEG,CURRENTSET, Conditions,CondNums,filter_range,filter_width,epoch,baseline,srate_downsample,input_bad_chans,ButterRp,ButterRs, ZscoreOn)

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


% reload data into EEGlab to use their notch filter for ERP analysis
% put the current data array back into the dim of EEG.data

for i=1:length(chan_lbls)
    EEG.data(chan_lbls(i),:) = data(:,i)';
end
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG  = pop_basicfilter( EEG,  1:size(EEG.data,1) , 'Cutoff',  60, 'Design', 'notch', 'Filter', 'PMnotch', 'Order',  180 );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

data_erp = double(EEG.data(chan_lbls,:))';

% EEG  = pop_basicfilter( EEG,  1 , 'Cutoff',  60, 'Design', 'notch', 'Filter', 'PMnotch', 'Order',  180 );
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
% EEG = eeg_checkset( EEG );


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

% % % % % 
% % % % % % Morlet Wavelets filters
% % % % % 
% % % % % hband_morlet_wide = zeros(size(data));
% % % % % for i=1:size(data,2)
% % % % %     disp(['chan ' int2str(i) ' of ' int2str(size(data,2))])
% % % % %     bb = bbpower(data(:,i),srate,filter_range,filter_width,30,false); %30
% % % % %     hband_morlet_wide(:,i)=bb;
% % % % %     clear bb
% % % % % end
% % % % % hband_morlet_narrow = zeros(size(data));
% % % % % for i=1:size(data,2)
% % % % %     disp(['chan ' int2str(i) ' of ' int2str(size(data,2))])
% % % % %     bb = bbpower(data(:,i),srate,filter_range,filter_width,6,false); %30
% % % % %     hband_morlet_narrow(:,i)=bb;
% % % % %     clear bb
% % % % % end

% Butterworth Filters

hband_buttord_wide=zeros(size(data));
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
    
% % % % %     for yy = 1:size(hband_filters,1)
% % % % %         hband = hband_filters(yy,:);
% % % % %         hband_sig_data(:,yy) = butterpass_eeglabdata(data(:,el),hband,srate,1,5); % winner
% % % % %         hband_sig_data(:,yy) = abs(hilbert(hband_sig_data(:,yy)));
% % % % %     end    
% % % % %     hband_buttord_narrow(:,el)=mean(hband_sig_data,2);
% % % % %     hband_sig_data=zeros([size(data,1) size(hband_filters,1)]);
% % % % %     
% % % % %     % fft filter as well
% % % % %     
% % % % %     tempEEG.data = EEG.data(chan_lbls(1,el),:);
% % % % %     disp(['filter el ' int2str(el) ' of ' int2str(size(data,2))])
% % % % %     for yy = 1:size(hband_filters,1)
% % % % %         hband = hband_filters(yy,:);
% % % % %         disp(['filter ' int2str(yy) ' of ' int2str(size(hband_filtersFIR,1))])
% % % % %         firtemp = pop_eegfilt( tempEEG, hband(1), hband(2), [], [0], 0, 0, 'fir1', 0);
% % % % %         hband_sig_data(:,yy) = firtemp.data';
% % % % %         hband_sig_data(:,yy) = abs(hilbert(squeeze(hband_sig_data(:,yy))));
% % % % %     end
% % % % %     FIRband_sig(:,el) =mean(hband_sig_data,2);
    
end
clear hband_sig_data
 
% create evoked activity as well

for el=1:size(hband_buttord_wide,2)
    evoked_downsmpl(:,el)=resample(data_erp(:,el),srate_new,srate);
end

% remake event (trigger) timeseries

% numevents = size(EEG.event,2);
% stim=zeros(numevents,2);
% for n=1:numevents
%     stim(n,:) = [EEG.event(1,n).latency+1 EEG.event(1,n).type];
% end

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
        base_line=squeeze(hband_sig_buttord_wide(elec,yyy,base1:base2));
        base_line=base_line(:);
        
        if ZscoreOn==1
            hband_sig_buttord_wide(elec,yyy,:)=(hband_sig_buttord_wide(elec,yyy,:)-mean(base_line))/...
                std(base_line);
        else
            hband_sig_buttord_wide(elec,yyy,:)=(hband_sig_buttord_wide(elec,yyy,:)-mean(base_line));
        end
% % % % %         hband_sig_buttord_narrow(elec,yyy,:)=...
% % % % %             hband_buttord_narrow(currsample+(epoch_adj(1)*srate_new):currsample+(epoch_adj(2)*srate_new),elec)';
% % % % %         base_line=squeeze(hband_sig_buttord_narrow(elec,yyy,1:epoch_adj(1)*srate_new*-1));
% % % % %         base_line=base_line(:);
% % % % %         hband_sig_buttord_narrow(elec,yyy,:)=(hband_sig_buttord_narrow(elec,yyy,:)-mean(base_line))/...
% % % % %             std(base_line);
        
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
        
        
% % % % %         FIRband_sig_induced(elec,yyy,:)=...
% % % % %             FIRband_sig(currsample+(epoch_adj(1)*srate_new):currsample+(epoch_adj(2)*srate_new),elec)';
% % % % %         base_line=squeeze(FIRband_sig_induced(elec,yyy,1:epoch_adj(1)*srate_new*-1));
% % % % %         base_line=base_line(:);
% % % % %         FIRband_sig_induced(elec,yyy,:)=(FIRband_sig_induced(elec,yyy,:)-mean(base_line))/...
% % % % %             std(base_line);
% % % % %         
% % % % %         hband_Wide_morlet_induced(elec,yyy,:)=...
% % % % %             hband_morlet_wide(currsample+(epoch_adj(1)*srate_new):currsample+(epoch_adj(2)*srate_new),elec)';
% % % % %         base_line=squeeze(hband_Wide_morlet_induced(elec,yyy,1:epoch_adj(1)*srate_new*-1));
% % % % %         base_line=base_line(:);
% % % % %         hband_Wide_morlet_induced(elec,yyy,:)=(hband_Wide_morlet_induced(elec,yyy,:)-mean(base_line))/...
% % % % %             std(base_line);
% % % % %         
% % % % %         hband_Narrow_morlet_induced(elec,yyy,:)=...
% % % % %             hband_morlet_narrow(currsample+(epoch_adj(1)*srate_new):currsample+(epoch_adj(2)*srate_new),elec)';
% % % % %         base_line=squeeze(hband_Narrow_morlet_induced(elec,yyy,1:epoch_adj(1)*srate_new*-1));
% % % % %         base_line=base_line(:);
% % % % %         hband_Narrow_morlet_induced(elec,yyy,:)=(hband_Narrow_morlet_induced(elec,yyy,:)-mean(base_line))/...
% % % % %             std(base_line);
        
    end
end


save (strcat(subject,'_',task));



% parse data into individual conditions

numevents_percond = zeros(length(Conditions),1);
for i=1:length(Conditions)
    numevents_percond(i) = size((find(stim_resamp==CondNums(i))),1);
    Induced_ButtOrd_Wide.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
% % % % %     Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
    Induced_evoked.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
% % % % %     Induced_FIRband_sig.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
% % % % %     Induced_Wide_morlet.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
% % % % % 	Induced_Narrow_morlet.(genvarname(char(Conditions(i)))) = zeros(num_el,trial_length); 
end
for i=1:numevents
    Induced_ButtOrd_Wide.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums))))) = Induced_ButtOrd_Wide.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums)))))+squeeze(hband_sig_buttord_wide(:,i,:));
% % % % %     Induced_ButtOrd_Narrow.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums))))) = Induced_ButtOrd_Narrow.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums)))))+squeeze(hband_sig_buttord_narrow(:,i,:));
    Induced_evoked.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums))))) = Induced_evoked.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums)))))+squeeze(evoked_downsmpl_induced(:,i,:));
% % % % %     Induced_FIRband_sig.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums))))) = Induced_FIRband_sig.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums)))))+squeeze(FIRband_sig_induced(:,i,:));
% % % % %     Induced_Wide_morlet.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums))))) = Induced_Wide_morlet.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums)))))+squeeze(hband_Wide_morlet_induced(:,i,:));
% % % % %     Induced_Narrow_morlet.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums))))) = Induced_Narrow_morlet.(genvarname(char(Conditions(find(stim_resamp(i,2)==CondNums)))))+squeeze(hband_Narrow_morlet_induced(:,i,:));
end
for i=1:length(Conditions)
    Induced_ButtOrd_Wide.(genvarname(char(Conditions(i)))) = [Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))/numevents_percond(i)]';
% % % % %     Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i)))) = [Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i))))/numevents_percond(i)]';
    Induced_evoked.(genvarname(char(Conditions(i)))) = [Induced_evoked.(genvarname(char(Conditions(i))))/numevents_percond(i)]';
% % % % %     Induced_FIRband_sig.(genvarname(char(Conditions(i)))) = [Induced_FIRband_sig.(genvarname(char(Conditions(i))))/numevents_percond(i)]';
% % % % %     Induced_Wide_morlet.(genvarname(char(Conditions(i)))) = [Induced_Wide_morlet.(genvarname(char(Conditions(i))))/numevents_percond(i)]';
% % % % %     Induced_Narrow_morlet.(genvarname(char(Conditions(i)))) = [Induced_Narrow_morlet.(genvarname(char(Conditions(i))))/numevents_percond(i)]';
    
    % Baseline all channels to 0
    
    for elec=1:num_el
        tempval = Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(base1:base2,elec);
        tempmean = mean(tempval);
        Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(:,elec) = Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(:,elec)-tempmean;
        Induced_ButtOrd_WideRS.(genvarname(char(Conditions(i))))(:,elec) = resample(Induced_ButtOrd_Wide.(genvarname(char(Conditions(i))))(:,elec),srate_downsample,srate);
        
% % % % %         tempval = Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i))))(base1:base2,elec);
% % % % %         tempmean = mean(tempval);
% % % % %         Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i))))(:,elec) = Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i))))(:,elec)-tempmean;
% % % % %         Induced_ButtOrd_NarrowRS.(genvarname(char(Conditions(i))))(:,elec) = resample(Induced_ButtOrd_Narrow.(genvarname(char(Conditions(i))))(:,elec),srate_downsample,srate);
        
        tempval = Induced_evoked.(genvarname(char(Conditions(i))))(base1:base2,elec);
        tempmean = mean(tempval);
        tempstd = std(tempval);
        if ZscoreOn==1
            Induced_evoked.(genvarname(char(Conditions(i))))(:,elec) = (Induced_evoked.(genvarname(char(Conditions(i))))(:,elec)-tempmean)/tempstd;
        else
            Induced_evoked.(genvarname(char(Conditions(i))))(:,elec) = (Induced_evoked.(genvarname(char(Conditions(i))))(:,elec)-tempmean);
        end
        Induced_evokedRS.(genvarname(char(Conditions(i))))(:,elec) = resample(Induced_evoked.(genvarname(char(Conditions(i))))(:,elec),srate_downsample,srate);
        
% % % % %         tempval = Induced_FIRband_sig.(genvarname(char(Conditions(i))))(base1:base2,elec);
% % % % %         tempmean = mean(tempval);
% % % % %         Induced_FIRband_sig.(genvarname(char(Conditions(i))))(:,elec) = Induced_FIRband_sig.(genvarname(char(Conditions(i))))(:,elec)-tempmean;
% % % % %         Induced_FIRband_sigRS.(genvarname(char(Conditions(i))))(:,elec) = resample(Induced_FIRband_sig.(genvarname(char(Conditions(i))))(:,elec),srate_downsample,srate);
% % % % %         
% % % % %         tempval = Induced_Wide_morlet.(genvarname(char(Conditions(i))))(base1:base2,elec);
% % % % %         tempmean = mean(tempval);
% % % % %         Induced_Wide_morlet.(genvarname(char(Conditions(i))))(:,elec) = Induced_Wide_morlet.(genvarname(char(Conditions(i))))(:,elec)-tempmean;
% % % % %         Induced_Wide_morletRS.(genvarname(char(Conditions(i))))(:,elec) = resample(Induced_Wide_morlet.(genvarname(char(Conditions(i))))(:,elec),srate_downsample,srate);
% % % % %         
% % % % %         tempval = Induced_Narrow_morlet.(genvarname(char(Conditions(i))))(base1:base2,elec);
% % % % %         tempmean = mean(tempval);
% % % % %         Induced_Narrow_morlet.(genvarname(char(Conditions(i))))(:,elec) = Induced_Narrow_morlet.(genvarname(char(Conditions(i))))(:,elec)-tempmean;
% % % % %         Induced_Narrow_morletRS.(genvarname(char(Conditions(i))))(:,elec) = resample(Induced_Narrow_morlet.(genvarname(char(Conditions(i))))(:,elec),srate_downsample,srate);
    end
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
    Count_Induced_evokedRS.(genvarname(char(Conditions(i)))) = 1;
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
    Epoch_Induced_evokedRS.(genvarname(char(Conditions(i)))) = [];
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
%         Epoch_Induced_Narrow_morlet.(genvarname(char(Conditions(sti m_resamp(i,2)))))(elec,:,Count_Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))) = Epoch_Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2)))))(elec,:,Count_Induced_Narrow_morlet.(genvarname(char(Conditions(stim_resamp(i,2))))))-tempmean;
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
    Count_Induced_evokedRS.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_evokedRS.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;
%     Count_Induced_FIRband_sigRS.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_FIRband_sigRS.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;
%     Count_Induced_Wide_morletRS.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_Wide_morletRS.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;
%     Count_Induced_Narrow_morletRS.(genvarname(char(Conditions(stim_resamp(i,2))))) = Count_Induced_Narrow_morletRS.(genvarname(char(Conditions(stim_resamp(i,2)))))+1;

end


% Rebaseline evoked epochsbecause i'm not sure it's working right up above
base1RS = round((baseline(1)/srate_downsample)-(epoch(1)/srate_downsample)+1);
base2RS = round((baseline(2)*srate_downsample)-(baseline(1)*srate_downsample));

for i = 1:length(Conditions)
    for epochloop =1:size(Epoch_Induced_evokedRS.(genvarname(char(Conditions(i)))),3)
        for elec = 1:size(Epoch_Induced_evokedRS.(genvarname(char(Conditions(i)))),1)
            tempval = Epoch_Induced_evokedRS.(genvarname(char(Conditions(i))))(elec,base1RS:base2RS,epochloop);
            tempmean = mean(tempval);
            tempstd = std(tempval);
            Epoch_Induced_evokedRS.(genvarname(char(Conditions(i))))(elec,:,epochloop) = ( Epoch_Induced_evokedRS.(genvarname(char(Conditions(i))))(elec,:,epochloop)-tempmean)/tempstd;
        end
    end
end



% save all data at the end
clear ALLCOM ALLEEG ALLERP ALLERPCOM CURRENTERP CURRENTSET CURRENTSTUDY EEG ERP LASTCOM PLUGINLIST STUDY eeglabUpdater plotset
clear data tempEEG 
save (strcat(subject,'_',task));

end





