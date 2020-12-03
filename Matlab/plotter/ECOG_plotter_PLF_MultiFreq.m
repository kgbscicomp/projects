% ECOG_plotter

% Enter condition numbers interested in plotting

NumConds = length(plot_conds);
% % if NumConds>2
% %     beep
% %     pause
% % end

plotbaseline = [round((plotbaseline(1)*srate)+1):round((plotbaseline(2)*srate)+1)];
zerotime = round((zerotime*srate)+1);
xaxis = [round((xaxis(1)*srate)+1) round((xaxis(2)*srate)+1)];
xticklabel = epoch(1):xticksize:epoch(2);
xtickrange = 1:srate*xticksize:size(hband_sig_buttord_wide,2);
colorarray = {'b','r','g'};

% create legend
legend_array = {};

for j=1:NumConds
    temp = find(CondNums==plot_conds(j));
    legend_array(length(legend_array)+1) = Conditions(temp);
end

% first run the calculations then plot in 2nd for loop

CondA_BUTTER_ALL = [];
CondB_BUTTER_ALL = [];
CondA_MORLET_ALL = [];
CondB_MORLET_ALL = [];
PermBUTTER_ALL_Chan = [];
PermMORLET_ALL_Chan = [];
CondA_BUTTER_ALL_Chan_STD = [];
CondB_BUTTER_ALL_Chan_STD = [];
CondA_MORLET_ALL_Chan_STD = [];
CondB_MORLET_ALL_Chan_STD = [];

PermBUTTER_ALL_Chan_thresh = [];
PermMORLET_ALL_Chan_thresh = [];
CurrCondTrialsA = find(stim_resamp(:,2)==plot_conds(1)); %pred
CurrCondTrialsB = find(stim_resamp(:,2)==plot_conds(2)); %nonpred

for chanx = 1:length(plot_chans) %[29 30 37 38 43:46 86:91]
    disp(['Channel ' int2str(chanx) ' of ' int2str(size(plot_chans,2))])
    
    tempPredBUTTER = [];
    tempNonPredBUTTER= [];
    tempPredMORLET = [];
    tempNonPredMORLET=[];
    for yyy=1:size(hband_filters,1)
        for j=1:trial_length
            tempPredBUTTER(:,j,yyy) = squeeze(hband_sig_buttord_wide(CurrCondTrialsA',j,yyy,chanx));
            tempNonPredBUTTER(:,j,yyy) = squeeze(hband_sig_buttord_wide(CurrCondTrialsB',j,yyy,chanx));
            
            tempPredMORLET(:,j,yyy) = squeeze(mband_sig_morlet_wide(CurrCondTrialsA',j,yyy,chanx));
            tempNonPredMORLET(:,j,yyy) = squeeze(mband_sig_morlet_wide(CurrCondTrialsB',j,yyy,chanx));
        end
    end
    
    %     ITC, trials x time x freq
    BUTTERnonPred = squeeze(abs(mean(tempNonPredBUTTER./abs(tempNonPredBUTTER))));
    BUTTERPred = squeeze(abs(mean(tempPredBUTTER./abs(tempPredBUTTER))));
    
    MORLETnonPred = squeeze(abs(mean(tempNonPredMORLET./abs(tempNonPredMORLET))));
    MORLETPred = squeeze(abs(mean(tempPredMORLET./abs(tempPredMORLET))));
    
    
    [ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
                   newtimef({EEG.data(1,:,:) ALLEEG(2).data(1,:,:)}, ...
                        EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, cycles);
    
    
    
    
    
    
    
    
    CondA_BUTTER_ALL(chanx,:,:)=BUTTERPred;
    CondB_BUTTER_ALL(chanx,:,:)=BUTTERnonPred;
    CondA_MORLET_ALL(chanx,:,:)=MORLETPred;
    CondB_MORLET_ALL(chanx,:,:)=MORLETnonPred;
    
    permdataBUTTER = zeros(size(hband_filters,1),size(MORLETnonPred,1));
    permdataMORLET= zeros(size(hband_filters,1),size(MORLETnonPred,1));
    realdataMORLET = [tempNonPredMORLET;tempPredMORLET];
    realdataBUTTER = [tempNonPredBUTTER;tempPredBUTTER];
    realdataBUTTERPred = [tempPredBUTTER];
    realdataBUTTERNonPred = [tempNonPredBUTTER];
        
    real_PredBUTTER =BUTTERPred;
    real_NonPredBUTTER = BUTTERnonPred;
    real_DiffBUTTER = real_PredBUTTER - real_NonPredBUTTER;
    
    real_PredMORLET =MORLETPred;
    real_NonPredMORLET = MORLETnonPred;
    real_DiffMORLET = real_PredMORLET - real_NonPredMORLET;
    
    for k=1:NumPerms % perm number
        for yyy=1:size(hband_filters,1)
            randY = randperm(size(CurrCondTrialsA,1)+size(CurrCondTrialsB,1));
            
            randBUTTER = realdataBUTTER(randY,:,yyy); %randi(64,1,32);
            test1BUTTER = randBUTTER(1:size(CurrCondTrialsA,1),:); % 24 trials
            test2BUTTER = randBUTTER(size(CurrCondTrialsA,1)+1:size(CurrCondTrialsA,1)+size(CurrCondTrialsB,1),:); % 22 trials
            
            randMORLET = realdataMORLET(randY,:,yyy); %randi(64,1,32);
            test1MORLET = randMORLET(1:size(CurrCondTrialsA,1),:);
            test2MORLET = randMORLET(size(CurrCondTrialsA,1)+1:size(CurrCondTrialsA,1)+size(CurrCondTrialsB,1),:);
            
            MORLET_test1 = squeeze(abs(mean(test1MORLET./abs(test1MORLET))));
            MORLET_test2 = squeeze(abs(mean(test2MORLET./abs(test2MORLET))));
            DiffMORLET = MORLET_test1-MORLET_test2;
            
            BUTTER_test1 = squeeze(abs(mean(test1BUTTER./abs(test1BUTTER))));
            BUTTER_test2 = squeeze(abs(mean(test2BUTTER./abs(test2BUTTER))));
            DiffBUTTER = BUTTER_test1-BUTTER_test2;
            
            for j=1:size(MORLETnonPred,1)
                if real_DiffBUTTER(j,yyy)<=DiffBUTTER(1,j)
                    permdataBUTTER(yyy,j) =  permdataBUTTER(yyy,j)+1;
                end
                if real_DiffMORLET(j,yyy)<=DiffMORLET(1,j)
                    permdataMORLET(yyy,j) =  permdataMORLET(yyy,j)+1;
                end
            end
        end
    end
    
    permdataBUTTER = permdataBUTTER/NumPerms;
    permdataMORLET = permdataMORLET/NumPerms;

    PermBUTTER_ALL_Chan(chanx,:,:)=permdataBUTTER;
    PermMORLET_ALL_Chan(chanx,:,:)=permdataMORLET;
    PermBUTTER_ALL_Chan_thresh(chanx,:,:) = zeros(size(permdataBUTTER,1),size(permdataBUTTER,2));
    PermMORLET_ALL_Chan_thresh(chanx,:,:) = zeros(size(permdataBUTTER,1),size(permdataBUTTER,2));
    for jjj=1:size(permdataBUTTER,2)
        for yyy=1:size(permdataBUTTER,1)
            if PermBUTTER_ALL_Chan(chanx,yyy,jjj)<=pthresh
                PermBUTTER_ALL_Chan_thresh(chanx,yyy,jjj) =  real_DiffBUTTER(jjj,yyy);
            end
            if PermMORLET_ALL_Chan(chanx,yyy,jjj)<=pthresh
                PermMORLET_ALL_Chan_thresh(chanx,yyy,jjj) =  real_DiffMORLET(jjj,yyy);
            end
        end
    end
    
    
    
    % bootstrap standard deviation for each condition alone
    
    clear STD;
    
    for i=1:NumPerms
        for yyy=1:size(permdataBUTTER,1)
            temp = datasample(tempPredBUTTER(:,:,yyy),size(tempPredBUTTER,1),1);
            STD.tempA_Z(i,:,yyy) = squeeze(abs(mean(squeeze(temp)./abs(squeeze(temp)))));
            temp = datasample(tempNonPredBUTTER(:,:,yyy),size(tempNonPredBUTTER,1),1);
            STD.tempB_Z(i,:,yyy) = squeeze(abs(mean(squeeze(temp)./abs(squeeze(temp)))));
            temp = datasample(tempPredMORLET(:,:,yyy),size(tempPredMORLET,1),1);
            STD.tempA_ITC(i,:,yyy) = squeeze(abs(mean(squeeze(temp)./abs(squeeze(temp)))));
            temp = datasample(tempNonPredMORLET(:,:,yyy),size(tempNonPredMORLET,1),1);
            STD.tempB_ITC(i,:,yyy) = squeeze(abs(mean(squeeze(temp)./abs(squeeze(temp)))));
        end
    end
    for yyy=1:size(permdataBUTTER,1)
        CondA_BUTTER_ALL_Chan_STD(chanx,:,yyy) = std(STD.tempA_Z(:,:,yyy));
        CondB_BUTTER_ALL_Chan_STD(chanx,:,yyy) = std(STD.tempB_Z(:,:,yyy));
        CondA_MORLET_ALL_Chan_STD(chanx,:,yyy) = std(STD.tempA_ITC(:,:,yyy));
        CondB_MORLET_ALL_Chan_STD(chanx,:,yyy) = std(STD.tempB_ITC(:,:,yyy));
    end
end
save(strcat(subject,'_',task),'-v7.3', '-regexp', '^(?!(ALLCOM|ALLEEG|ALLERP|ALLERPCOM|CURRENTERP|CURRENTSET|CURRENTSTUDY|EEG|ERP|LASTCOM|PLUGINLIST|STUDY|eeglabUpdater|plotset|data|tempEEG)$).');

% zerobar = zeros(1,size(CondA_BUTTER_ALL(chanx,:),2));
% zerobar((zerotime)) = 1;

for chanx = 1:length(plot_chans) %[29 30 37 38 43:46 86:91]
    
    figure(5)
    temp = squeeze(PermBUTTER_ALL_Chan_thresh(chanx,:,:));
    imagesc(flipud(temp));figure(gcf);
    
%     bar([PermBUTTER_ALL_Chan(chanx,:)',PermBUTTER_ALL_Chan_thresh(chanx,:)'])
%     hold on
%     colormap([.9 .9 .9;0 .8 0])
%     shadedErrorBar(1:size(CondA_BUTTER_ALL(chanx,:),2),CondA_BUTTER_ALL(chanx,:)',CondA_BUTTER_ALL_Chan_STD(chanx,:)',colorarray(1),.5);
%     shadedErrorBar(1:size(CondB_BUTTER_ALL(chanx,:),2),CondB_BUTTER_ALL(chanx,:)',CondB_BUTTER_ALL_Chan_STD(chanx,:)',colorarray(2),.5);
%     plot(zeros(size(CondA_BUTTER_ALL(chanx,:),2),1)+.05);
%     plot(zerobar)
%     hold off
    title(strcat(char(chan_names_new((chanx))),' ITC BUTTERWORTH'),'FontWeight','bold');
%     set(gca, 'XTick', xtickrange)
%     set(gca, 'XTickLabel', xticklabel)
%     xlim([xaxis(1) xaxis(2)])    
%     ylim([0 1])
%     legend(legend_array);
    
    figure(6)
    temp = squeeze(PermMORLET_ALL_Chan_thresh(chanx,:,:));
    imagesc(flipud(temp));figure(gcf);
    
%     bar([PermMORLET_ALL_Chan(chanx,:)',PermMORLET_ALL_Chan_thresh(chanx,:)'])
%     hold on
%     colormap([.9 .9 .9;0 .8 0])
%     shadedErrorBar(1:size(CondA_ITC_ALL(chanx,:),2),CondA_ITC_ALL(chanx,:)',CondA_ITC_ALL_Chan_STD(chanx,:)',colorarray(1),.5);
%     shadedErrorBar(1:size(CondB_ITC_ALL(chanx,:),2),CondB_ITC_ALL(chanx,:)',CondB_ITC_ALL_Chan_STD(chanx,:)',colorarray(2),.5);
%     plot(zeros(size(CondA_ITC_ALL(chanx,:),2),1)+.05);
%     plot(zerobar)
%     hold off
    title(strcat(char(chan_names_new((chanx))),' ITC Wavelets'),'FontWeight','bold');
%     set(gca, 'XTick', xtickrange)
%     set(gca, 'XTickLabel', xticklabel)
%     xlim([xaxis(1) xaxis(2)])    
%     ylim([0 1])
%     legend(legend_array);
    
    pause();
end


