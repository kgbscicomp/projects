% ECOG_plotter

% Enter condition numbers interested in plotting

if PLOTONLY >1
else
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
    
    PermBUTTER_ALL_Chan_threshval = [];
    PermMORLET_ALL_Chan_thresh = [];
    CurrCondTrialsA = find(stim_resamp(:,2)==plot_conds(1)); %pred
    CurrCondTrialsB = find(stim_resamp(:,2)==plot_conds(2)); %nonpred
    tempPredMORLET = [];
    tempNonPredMORLET=[];
    
    for chanx = 1:length(plot_chans) %[29 30 37 38 43:46 86:91]
        disp(['Channel ' int2str(chanx) ' of ' int2str(size(plot_chans,2))])
        
        tempPredBUTTER = [];
        tempNonPredBUTTER= [];
        
        for j=1:trial_length
            for yyy=1:size(hband_filters,1)
                tempPredBUTTER(:,j,yyy) = squeeze(hband_sig_buttord_wide(CurrCondTrialsA',j,yyy,chanx));
                tempNonPredBUTTER(:,j,yyy) = squeeze(hband_sig_buttord_wide(CurrCondTrialsB',j,yyy,chanx));
            end
            tempPredMORLET(chanx,j,:) = squeeze(eeglab_sig_wide(CurrCondTrialsA',j,chanx));
            tempNonPredMORLET(chanx,j,:) = squeeze(eeglab_sig_wide(CurrCondTrialsB',j,chanx));
        end
        
        %     ITC, trials x time x freq
        BUTTERnonPred = squeeze(abs(mean(tempNonPredBUTTER./abs(tempNonPredBUTTER))));
        BUTTERPred = squeeze(abs(mean(tempPredBUTTER./abs(tempPredBUTTER))));
        
        
        CondA_BUTTER_ALL(chanx,:,:)=BUTTERPred;
        CondB_BUTTER_ALL(chanx,:,:)=BUTTERnonPred;
        
        permdataBUTTER = zeros(size(hband_filters,1),size(BUTTERnonPred,1));
        realdataBUTTER = [tempNonPredBUTTER;tempPredBUTTER];
        
        realdataBUTTERPred = [tempPredBUTTER];
        realdataBUTTERNonPred = [tempNonPredBUTTER];
        
        real_PredBUTTER =BUTTERPred;
        real_NonPredBUTTER = BUTTERnonPred;
        real_DiffBUTTER = real_PredBUTTER - real_NonPredBUTTER;
        
        for k=1:NumPerms % perm number
            for yyy=1:size(hband_filters,1)
                randY = randperm(size(CurrCondTrialsA,1)+size(CurrCondTrialsB,1));
                
                randBUTTER = realdataBUTTER(randY,:,yyy); %randi(64,1,32);
                test1BUTTER = randBUTTER(1:size(CurrCondTrialsA,1),:); % 24 trials
                test2BUTTER = randBUTTER(size(CurrCondTrialsA,1)+1:size(CurrCondTrialsA,1)+size(CurrCondTrialsB,1),:); % 22 trials
                
                BUTTER_test1 = squeeze(abs(mean(test1BUTTER./abs(test1BUTTER))));
                BUTTER_test2 = squeeze(abs(mean(test2BUTTER./abs(test2BUTTER))));
                DiffBUTTER = BUTTER_test1-BUTTER_test2;
                
                for j=1:size(DiffBUTTER,2)
                    if real_DiffBUTTER(j,yyy)<=DiffBUTTER(1,j)
                        permdataBUTTER(yyy,j) =  permdataBUTTER(yyy,j)+1;
                    end
                end
            end
        end

        permdataBUTTER = permdataBUTTER/NumPerms;
        
        PermBUTTER_ALL_Chan(chanx,:,:)=permdataBUTTER;
        PermBUTTER_ALL_Chan_threshval(chanx,:,:) = zeros(size(permdataBUTTER,1),size(permdataBUTTER,2));
        PermBUTTER_ALL_Chan_thresh(chanx,:,:) = zeros(size(permdataBUTTER,1),size(permdataBUTTER,2));
        
        for jjj=1:size(permdataBUTTER,2)
            for yyy=1:size(permdataBUTTER,1)
                if PermBUTTER_ALL_Chan(chanx,yyy,jjj)<=pthresh
                    PermBUTTER_ALL_Chan_threshval(chanx,yyy,jjj) =  real_DiffBUTTER(jjj,yyy);
                    PermBUTTER_ALL_Chan_thresh(chanx,yyy,jjj) =  PermBUTTER_ALL_Chan(chanx,yyy,jjj);
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
            end
        end
        for yyy=1:size(permdataBUTTER,1)
            CondA_BUTTER_ALL_Chan_STD(chanx,:,yyy) = std(STD.tempA_Z(:,:,yyy));
            CondB_BUTTER_ALL_Chan_STD(chanx,:,yyy) = std(STD.tempB_Z(:,:,yyy));
        end
    end
end
% save(strcat(subject,'_',task),'-v7.3', '-regexp', '^(?!(ALLCOM|ALLEEG|ALLERP|ALLERPCOM|CURRENTERP|CURRENTSET|CURRENTSTUDY|EEG|ERP|LASTCOM|PLUGINLIST|STUDY|eeglabUpdater|plotset|data|tempEEG)$).');

zerobar = zeros(1,size(CondA_BUTTER_ALL,2));
zerobar((zerotime)) = 1;

for chanx = 1:length(plot_chans) %[29 30 37 38 43:46 86:91]
    
    figure(5)
    temp = squeeze(PermBUTTER_ALL_Chan(chanx,:,:));
    imagesc(flipud(temp));figure(gcf);
    ylabel('Frequencies');
    set(gca,'YTick',1:size(hband_filters,1))
    for i=1:size(hband_filters,1),YTickNames{i}=strcat(num2str(hband_filters(i,1)),'-',num2str(hband_filters(i,2)));end
    set(gca,'YTickLabel',fliplr(YTickNames))
    set(gca,'XLim',xaxis)
    set(gca,'XTick',xaxis(1):round(srate*xticksize):xaxis(2))
    set(gca,'XTickLabel',(xaxis(1)-zerotime)/srate:xticksize:(xaxis(2)-zerotime)/srate)
    caxis([-1 1])
    
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
    xmin_temp = -(zerotime-xaxis(1))/srate;
    xmax_temp = (xaxis(2)-zerotime)/srate;
    [ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
        newtimef({tempPredMORLET(chanx,xaxis(1):xaxis(2),:) tempNonPredMORLET(chanx,xaxis(1):xaxis(2),:)}, ...
        xaxis(2), [xmin_temp xmax_temp]*1000, EEG.srate, ...
        0,'alpha',.05, 'plotphasesign','off','title',{strcat(char(chan_names_new((chanx))),' FFT ITC Pred'), ...
        strcat(char(chan_names_new((chanx))),' FFT ITC NonPred')}, ...
        'newfig','off','freqs',[hband_filters(1,1),hband_filters(end,2)],'nfreqs',size(hband_filters,1), ...
        'boottype','shuffle','naccu',NumPerms,'baseboot',0);
    
    pause();
end

if PLOTONLY==2
    for chanx = 1:length(plot_chans)
        for freqx = plot_freqs
            
            figure(10+freqx)
            bar([squeeze(PermBUTTER_ALL_Chan(chanx,freqx,:)),squeeze(PermBUTTER_ALL_Chan_thresh(chanx,freqx,:))]);
            hold on
            colormap([.9 .9 .9;0 .8 0])
            shadedErrorBar(1:size(CondA_BUTTER_ALL,2),CondA_BUTTER_ALL(chanx,:,freqx),CondA_BUTTER_ALL_Chan_STD(chanx,:,freqx),colorarray(1),.5);
            shadedErrorBar(1:size(CondB_BUTTER_ALL,2),CondB_BUTTER_ALL(chanx,:,freqx),CondB_BUTTER_ALL_Chan_STD(chanx,:,freqx),colorarray(2),.5);
            plot(zeros(size(CondA_BUTTER_ALL,2),1)+.05);
            plot(zerobar)
            hold off
            title(strcat(char(chan_names_new((chanx))),' ITC BUTTERWORTH: ',num2str(hband_filters(freqx,1)),'-',num2str(hband_filters(freqx,2)),'hz'),'FontWeight','bold');
            set(gca, 'XTick', xtickrange)
            set(gca, 'XTickLabel', xticklabel)
            xlim([xaxis(1) xaxis(2)])
            ylim([0 1])
            legend(legend_array);
            
            
            
        end
        pause();
    end 
end






