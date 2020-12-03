
function ECOG_filter_analysis(subject,task,ALLEEG, EEG,CURRENTSET, Conditions,CondNums,input_bad_chans,epoch,reference,Ref_electrode,CarOption)
%Filter analysis for ECOG project

% % Conditions correspond to the event numbers in eeglab (e.g., 1-7 for DF)
% Conditions = {'2beeps','1flash','2flashes','Illusion','Ctrl','Blank','2flashes2beeps'};
% CondNums = 1:7; % cond numbers to look for
% Reference: 
% 0 for common average ref, 1 = lowest variance electrode, 2 = laplacian, 3 = common avg by grid, 4 = input electrode, 5 = no ref
% ref_electrode is optional
% CarOption: Defaults to 1 which does not remove electrode from itself to
% take common average ref, set to 0 to change this

%% Still need to implement

% For finding bad channels need to have the data selected be restricted to
% the epoch around events (find unique time in case overlapping epochs)
% We won't want to just concatenate the arrays since that could lead to
% artificial artifacts from jumpts in datapoints between samples. mean
% std across samples for each channel. effectively create artificial
% epochs.



%%

%Load data and remove input bad channels before going on to referencing
srate = round(EEG.srate);
data = EEG.data';
chan_lbls_orig=1:size(data,2);
chan_names_orig = {EEG.chanlocs(chan_lbls_orig).labels};

if nargin < 9
    epoch = [-2 2];
end

%Make ref_electrode optional input
if ~exist('Ref_electrode','var') || reference ~= 4
    Ref_electrode=0;
end
Ref_index = find(chan_lbls_orig == Ref_electrode);

if ~exist('CarOption','var')
    CarOption=1;
end

%% Initial step is to look for bad electrodes & reference
% For all techniques, data will be referenced and then bad electrodes
% selected 2 times. Referencing is not applied to the actual data until the
% final step after electrodes are removed

% 12/5/16
% Updated to be iterative process

RejThresh = 3; % Defaut SD threshold to find a bad channel
while RejThresh>0
    
    
    % second loop is to iteratively look for bad channels
    bad_chans_temp1 = [];
    loop_counter = 1;

    % Create figure to display bad channels
    figure(2);
    hold
    
    % Calculate appropriate spacing between plotted bad channels
    plot_diff = data; 
    plot_diff(:,input_bad_chans) = []; 
    chan_lbls = chan_lbls_orig;
    chan_lbls(:,input_bad_chans) = []; %remove input bad chans
    chan_names = {EEG.chanlocs(chan_lbls).labels};
    [plot_diff, ~] = ECOG_Reference(plot_diff, reference, Ref_electrode, CarOption, chan_names, chan_lbls);
    plot_diff = (mean(abs(min(plot_diff)-max(plot_diff))))/2;
    plot_time = 10; % secs
    
    while 1
        % Reset data and channel labels at start of each loop
        tempdata = data;
        chan_lbls = chan_lbls_orig;
        
        % Remove bad channels
        bad_chans_to_remove = [bad_chans_temp1 input_bad_chans];
        tempdata(:,bad_chans_to_remove) = []; %remove input bad chans
        chan_lbls(:,bad_chans_to_remove) = []; %remove input bad chans
        chan_names = {EEG.chanlocs(chan_lbls).labels};
        
        % find orig indices of bad channels
        bad_chan_index = 1:size(data,2);
        bad_chan_index(bad_chans_to_remove)=[];

        %Do initial referencing before looking for additional bad channels
        [tempdata, Ref_electrode_selected] = ECOG_Reference(tempdata, reference, Ref_electrode, CarOption, chan_names, chan_lbls);
        
        %First Pass
        %look for noisy channels
        
        a = std(tempdata);
        mean_var = mean(a);
        std_var = std(a);
        b = find(a>(mean_var+(RejThresh*std_var)));
        c = find(a<(mean_var-(RejThresh*std_var)));
        bad_chans = unique([b c]);
        
        
        bad_chans_temp1 = [bad_chans_temp1 bad_chan_index(bad_chans)];
        
        % Plot
        % Add 100 to separate plots
        % Should be in chronological order of reject and correspond to
        % order displayed below
        %figure;
        for i=1:length(bad_chans)
            plot(EEG.times(1:round(EEG.srate*plot_time)),tempdata(1:round(EEG.srate*plot_time),bad_chans(i))- mean(tempdata(1:5000,bad_chans(i)))-((length(bad_chans_temp1)-(length(bad_chans)))*plot_diff)-((i-1)*plot_diff))
            text(round(EEG.srate*plot_time),-((length(bad_chans_temp1)-(length(bad_chans)))*plot_diff)-((i-1)*plot_diff),strcat(chan_names_orig(bad_chan_index(bad_chans(i)))),'HorizontalAlignment','left');
        end

        if isempty(bad_chans)
            break;
        end
    end
    hold off
    
    %disp([' bad channels additional ind: ' int2str(bad_chans_temp1)])
    disp(char([' bad channels additional orig lbls: ' strcat(chan_names_orig(bad_chans_temp1))]))
    
    disp([' Current Threshold: ' num2str(RejThresh)])
    RejThresh= input('Enter in 0 to exit or new threshold: ');
    
    % if figure still open, close it
    if ishandle(2)==1
        close 2
    end
end
evoked = tempdata;
clear tempdata;

%% remake event (trigger) timeseries
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

%% create epochs
trial_length = round((epoch(2)-epoch(1))*srate)+1;
num_el = size(evoked,2);
evoked_sig = [];
for yyy=1:numevents
    disp(['trial ' int2str(yyy) ' of ' int2str(numevents)])
    currsample = round(stim(yyy,1));
    tempwindow = round([(epoch(1)*srate) (epoch(2)*srate)]);
    temprng = currsample+tempwindow(1):currsample+tempwindow(2);
    temp=double(evoked(temprng,:));
    temp2 = [];
    for i=1:size(temp,2)
        temp2(:,i) = temp(:,i)-mean(temp(1:1+(epoch(1)*-srate),i));
    end
    evoked_sig(:,:,yyy)=temp2';
end


%% parse data into individual conditions and save
disp(['Parsing data into individual conditions'])

for i=1:length(Conditions)
    temp = find(stim(:,2)==CondNums(i));
    Evoked.(genvarname(char(Conditions(i)))) = squeeze(mean(evoked_sig(:,:,temp),3));
end

%Parse to conditions with individual trials
for i=1:length(Conditions)
    temp = find(stim(:,2)==CondNums(i));
    EvokedWithTrials.(genvarname(char(Conditions(i)))) = squeeze((evoked_sig(:,:,temp)));
end

disp('Saving Data')

% save all data at the end
save(strcat(subject,'_',task),'-v7.3', '-regexp', '^(?!(ALLCOM|ALLEEG|ALLERP|ALLERPCOM|CURRENTERP|CURRENTSET|CURRENTSTUDY|EEG|ERP|LASTCOM|PLUGINLIST|STUDY|eeglabUpdater|plotset|data|evoked)$).');

end



