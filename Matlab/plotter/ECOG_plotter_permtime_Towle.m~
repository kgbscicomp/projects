% ECOG_plotter

% Enter condition numbers interested in plotting

NumConds = length(plot_conds);

plotbaseline = [round((plotbaseline(1)*srate)+1):round((plotbaseline(2)*srate)+1)];
zbaseline = [round((zbaseline(1)*srate)+1):round((zbaseline(2)*srate)+1)];
zerotime = round((zerotime*srate)+1);
xaxis = [round((xaxis(1)*srate)+1) round((xaxis(2)*srate)+1)];
xticklabel = epoch(1):xticksize:epoch(2);
xtickrange = 1:srate*xticksize:size(hband_data_all,3);
colorarray = {'r','g','b','y'};
if exist('bad_trials') == 0, bad_trials = [];end

% create legend
legend_array = {};

for j=1:NumConds
    temp = find(CondNums==plot_conds(j));
    legend_array(length(legend_array)+1) = Conditions(temp);
end


% compare the first two, but still plot any others

if mean(xaxis)==zerotime
else
%     beep
%     pause
    disp('Time window needs to be symmetrical around zero')
end

PLOT.STATS_C = zeros(max(plot_chans),round((xaxis(2)-xaxis(1)+1)));
PLOT.STATS_ERP = zeros(max(plot_chans),round((xaxis(2)-xaxis(1)+1)));
PLOT.STATS_C_thresh = zeros(max(plot_chans),1)-1;
PLOT.STATS_ERP_thresh = zeros(max(plot_chans),1)-1;

PLOT.PermResultsERP = [];
PLOT.PermResultsC = [];
PLOT.PermThreshERP = [];
PLOT.PermThreshC = [];

for j=1:NumConds
    PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_ALL')))) = [];
    PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_ALL')))) = [];
    PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_ALL_err')))) = [];
    PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_ALL_err')))) = [];
end
for chanx = plot_chans %[29 30 37 38 43:46 86:91]
    
    for j=1:NumConds
        
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))) = [];
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))) = [];
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_err')))) = [];
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))) = [];
        
        for i=1:numevents
            if Induced_chan_bad_array(chanx,i)==0
                if stim_resamp(i,2)==plot_conds(j)
                    temp = squeeze(hband_data_all(chanx,:,:,i));
                    PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,:,tempsizeC) = temp;
                    tempsizeC=tempsizeC+1;
                end
            end
            if Evoked_chan_bad_array(chanx,i)==0
                if stim_resamp(i,2)==plot_conds(j)
                    temp = squeeze(evoked_sig(chanx,:,i));
                    tempsize = size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)+1;
                    PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,tempsize) = temp;                
                end
            end
        end
        
        for i=1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)
            PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(plotbaseline,i),1));
        end 
        
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_err')))) = std(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))')/sqrt(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2));
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))) = std(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))')/sqrt(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2));
        
        for i=1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2)
            PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(plotbaseline,i),1));
        end
        for i=1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)
            PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(plotbaseline,i),1));
        end
        
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_ALL'))))(chanx,:) = mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2)';
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_ALL'))))(chanx,:) = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)';
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_ALL_err'))))(chanx,:) = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_err'))));
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_ALL_err'))))(chanx,:) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err'))));
        
    end
    
    % calculate one sample t-stat for each timepoint of cond 1
    
    [H,P,CI,STATS] = ttest(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))))');
    RealT_ERP = STATS.tstat;
    RealP_ERP = P;
    [H,P,CI,STATS] = ttest(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1))))))');
    RealT_C = STATS.tstat;
    RealP_C = P;
    
    tmax_ERP = zeros(1,size(RealT_ERP,2));
    tmax_C = zeros(1,size(RealT_C,2));
    temp_ERP = 0;
    temp_C = 0;
    counter_ERP = 1;
    counter_C = 1;
    
    % calc real tmass
    for i=zerotime+1:size(RealT_ERP,2)
        if ((RealT_ERP(i)<0 && RealT_ERP(i-1)>0)||(RealT_ERP(i)>0 && RealT_ERP(i-1)<0))
            temp_ERP = 0;
            counter_ERP = 1;
        end
        if RealP_ERP(i)<=pthresh
            temp_ERP = temp_ERP+RealT_ERP(i);
            counter_ERP = counter_ERP-1;
        else
            temp_ERP = 0;
            counter_ERP = 1;
        end
        tmax_ERP(i) = temp_ERP;
        if counter_ERP < 0
           for j=counter_ERP:-1
               tmax_ERP(i+j) = temp_ERP;
           end
        end
    end
    PLOT.STATS_ERP(chanx,:) = tmax_ERP;
    
    for i=zerotime+1:size(RealT_C,2)
        if ((RealT_C(i)<0 && RealT_C(i-1)>0)||(RealT_C(i)>0 && RealT_C(i-1)<0))
            temp_C = 0;
            counter_C = 1;
        end
        if RealP_C(i)<=pthresh
            temp_C = temp_C+RealT_C(i);
            counter_C = counter_C-1;
        else
            temp_C = 0;
            counter_C = 1;
        end
        tmax_C(i) = temp_C;
        if counter_C < 0
           for j=counter_C:-1
               tmax_C(i+j) = temp_C;
           end
        end
    end
    PLOT.STATS_C(chanx,:) = tmax_C;
    
    
    % Perm
    test_smpl = 1:zerotime;
    for zzz=1:NumPerms
        
        Smpl_Order = [];
        for yyy=1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2)
            randX = randi(length(test_smpl));
            Smpl_Order(:,yyy) = [test_smpl(randX:end) test_smpl(1:randX-1)];
        end
        
        [H,P,CI,STATS] = ttest(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))))(Smpl_Order)');
        RealT_ERP = STATS.tstat;
        RealP_ERP = P;
        [H,P,CI,STATS] = ttest(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1))))))(Smpl_Order)');
        RealT_C = STATS.tstat;
        RealP_C = P;
        
        tmax_ERP = zeros(1,size(RealT_ERP,2));
        tmax_C = zeros(1,size(RealT_C,2));
        temp_ERP = 0;
        temp_C = 0;
        counter_ERP = 1;
        counter_C = 1;
        
        % calc perm tmass
        for i=1:zerotime
            if i>1
                if ((RealT_ERP(i)<0 && RealT_ERP(i-1)>0)||(RealT_ERP(i)>0 && RealT_ERP(i-1)<0))
                    temp_ERP = 0;
                    counter_ERP = 1;
                end
            end
            if RealP_ERP(i)<=pthresh
                temp_ERP = temp_ERP+RealT_ERP(i);
                counter_ERP = counter_ERP-1;
            else
                temp_ERP = 0;
                counter_ERP = 1;
            end
            tmax_ERP(i) = temp_ERP;
            if counter_ERP < 0
                for j=counter_ERP:-1
                    tmax_ERP(i+j) = temp_ERP;
                end
            end
        end
        PLOT.PermResultsERP(chanx,zzz) = max(abs(tmax_ERP));
        
        for i=1:zerotime
            if i>1
                if ((RealT_C(i)<0 && RealT_C(i-1)>0)||(RealT_C(i)>0 && RealT_C(i-1)<0))
                    temp_C = 0;
                    counter_C = 1;
                end
            end
            if RealP_C(i)<=pthresh
                temp_C = temp_C+RealT_C(i);
                counter_C = counter_C-1;
            else
                temp_C = 0;
                counter_C = 1;
            end
            tmax_C(i) = temp_C;
            if counter_C < 0
                for j=counter_C:-1
                    tmax_C(i+j) = temp_C;
                end
            end
        end
        PLOT.PermResultsC(chanx,zzz) = max((tmax_C)); % not abs valued because we only want induced above baseline
        
    end
    
    CutoffVal = quantile(1:NumPerms,1-pthresh)+.5;
    temp = sort(PLOT.PermResultsERP(chanx,:));
    PLOT.PermThreshERP(chanx,1) = temp(CutoffVal);
    temp = sort(PLOT.PermResultsC(chanx,:));
    PLOT.PermThreshC(chanx,1) = temp(CutoffVal);
    
    temp = find(PLOT.STATS_C(chanx,:)>PLOT.PermThreshC(chanx,1),1); % not abs valued because we only want induced above baseline
    if ~isempty(temp)
        PLOT.STATS_C_thresh(chanx,1) = temp;
    end
    temp = find(abs(PLOT.STATS_ERP(chanx,:))>PLOT.PermThreshERP(chanx,1),1);
    if ~isempty(temp)
        PLOT.STATS_ERP_thresh(chanx,1) = temp;
    end
    disp(['Channel Number ' num2str(chanx)])
    
end

% calculate time from baseline

PLOT.STATS_C_thresh_time = PLOT.STATS_C_thresh;
PLOT.STATS_ERP_thresh_time = PLOT.STATS_ERP_thresh;
for i=1:size(PLOT.STATS_C_thresh,1),
    if PLOT.STATS_C_thresh_time(i,1) == -1,
    else
        PLOT.STATS_C_thresh_time(i,1) = (PLOT.STATS_C_thresh_time(i,1)-base2)/srate;
    end
    if PLOT.STATS_ERP_thresh_time(i,1) == -1,
    else
        PLOT.STATS_ERP_thresh_time(i,1) = (PLOT.STATS_ERP_thresh_time(i,1)-base2)/srate;
    end
end

% use the calculated std error
% first rebaseline close to onset
% % % % % rebase = [(base2-round(.1*srate)) base2];
% % % % % for chanx = plot_chans
% % % % %     PLOT.REBASE_
% % % % %     
% % % % % end
PLOT.STATS_C_thresh_time_std = zeros([size(PLOT.STATS_C_thresh)])-1;
PLOT.STATS_ERP_thresh_time_std = zeros([size(PLOT.STATS_ERP_thresh)])-1;
for i=1:size(PLOT.C1_ALL_err,1)
    foundstart = 0;
    for j=base2+1:size(PLOT.C1_ALL_err,2)
        tempAVG = PLOT.C1_ALL(i,j);
        tempSTD = PLOT.C1_ALL_err(i,j);
        if (tempAVG > 0 + abs(STD_thresh*tempSTD)) && foundstart == 0
            tempTIME = (j-base2)/srate;
            PLOT.STATS_C_thresh_time_std(i,1) = tempTIME;
            foundstart =1;
        end
    end
end
for i=1:size(PLOT.ERP1_ALL_err,1)
    foundstart = 0;
    for j=base2+1:size(PLOT.ERP1_ALL_err,2)
        tempAVG = PLOT.ERP1_ALL(i,j);
        tempSTD = PLOT.ERP1_ALL_err(i,j);
        if (tempAVG > 0 + abs(STD_thresh*tempSTD)) && foundstart == 0
            tempTIME = (j-base2)/srate;
            PLOT.STATS_ERP_thresh_time_std(i,1) = tempTIME;
            foundstart =1;
        end
    end
end





zerobarInduced = zeros(1,size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1));
zerobarERP = zeros(1,size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1));
zerobarInduced(zerotime) = min(mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2)); % 0 ms from aud onset
zerobarERP(zerotime) = min(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2))/2; % 0 ms from aud onset
save(strcat(subject,'_',task),'-v7.3', '-regexp', '^(?!(ALLCOM|ALLEEG|ALLERP|ALLERPCOM|CURRENTERP|CURRENTSET|CURRENTSTUDY|EEG|ERP|LASTCOM|PLUGINLIST|STUDY|eeglabUpdater|plotset|data|tempEEG)$).');


for chanx = plot_chans %[29 30 37 38 43:46 86:91]
    
% % %     figure(6)
% % %     if PLOT.STATS_C_thresh(chanx,1)>0
% % %         temp = zeros(1,size(PLOT.STATS_C,2));
% % %         temp(1:PLOT.STATS_C_thresh(chanx,1)) = 1.1*max((PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)),'_ALL'))))(chanx,:)));
% % %         bar(temp)
% % %         colormap([.9 1 .9])
% % %         hold on
% % %     end
% % %     
% % % 
% % %     shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1),(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)),'_ALL'))))(chanx,:)),PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)),'_err')))),colorarray(1),.5);
% % %     hold on
% % %     if NumConds>1
% % %         for j=2:NumConds
% % %             shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1),(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_ALL'))))(chanx,:)),PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_err')))),colorarray(2),.5);
% % %         end
% % %     end
% % %     
% % %     plot(zerobarInduced)
% % %     hold off
% % %     title(strcat('Induced CH:   ', chan_names(chanx)),'FontWeight','bold');
% % %     set(gca, 'XTick', xtickrange)
% % %     set(gca, 'XTickLabel', xticklabel)
% % %     xlim([xaxis(1) xaxis(2)])
% % %     legend(legend_array,'Location','SouthWest');
% % %     


    figure(6)
    if PLOT.STATS_C_thresh(chanx,1)>0
        temp = zeros(1,size(PLOT.STATS_C,2));
        temp(1:PLOT.STATS_C_thresh(chanx,1)) = 1.1*max((PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)),'_ALL'))))(chanx,:)));
        bar(temp)
        colormap([.8 1 .8])
        hold on
    end
    
    h(1) = plot((PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)),'_ALL'))))(chanx,:)),'Color',char(colorarray(1)));
    hold on
    shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1),(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)),'_ALL'))))(chanx,:)),PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)),'_ALL_err'))))(chanx,:),colorarray(1),.5);
    if NumConds>1
        for j=2:NumConds
            h(j) = plot((PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_ALL'))))(chanx,:)),'Color',char(colorarray(j)));
            shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1),(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_ALL'))))(chanx,:)),PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_ALL_err'))))(chanx,:),colorarray(j),.5);
        end
    end
    
    plot(zerobarInduced)
    hold off
    title(strcat('Induced CH:   ', chan_names(chanx)),'FontWeight','bold');
    set(gca, 'XTick', xtickrange)
    set(gca, 'XTickLabel', xticklabel)
    xlim([xaxis(1) xaxis(2)])
    legend(h,legend_array,'Location','SouthWest');
    
    
    
    figure(7)
    if PLOT.STATS_ERP_thresh(chanx,1)>0
        temp = zeros(1,size(PLOT.STATS_ERP,2));
        temp(1:PLOT.STATS_ERP_thresh(chanx,1)) = 1.1*max((PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_ALL'))))(chanx,:)));
        bar(temp)
        colormap([.8 1 .8])
        hold on
    end
    
    h(1) = plot((PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_ALL'))))(chanx,:)),'Color',char(colorarray(1)));
    hold on
    shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),1),(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_ALL'))))(chanx,:)),PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_ALL_err'))))(chanx,:),colorarray(1),.5);
    if NumConds>1
        for j=2:NumConds
            h(j) = plot((PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_ALL'))))(chanx,:)),'Color',char(colorarray(j)));
            shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),1),(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_ALL'))))(chanx,:)),PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_ALL_err'))))(chanx,:),colorarray(j),.5);
        end
    end
    
    plot(zerobarInduced)
    hold off
    title(strcat('Induced CH:   ', chan_names(chanx)),'FontWeight','bold');
    set(gca, 'XTick', xtickrange)
    set(gca, 'XTickLabel', xticklabel)
    xlim([xaxis(1) xaxis(2)])
    legend(h,legend_array,'Location','SouthWest');
    
    
    
% % % %     
% % % %     figure(7)
% % % %     if PLOT.STATS_ERP_thresh(chanx,1)>0
% % % %         temp = zeros(1,size(PLOT.STATS_ERP,2));
% % % %         temp(1:PLOT.STATS_ERP_thresh(chanx,1)) = 1.1*max((PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_ALL'))))(chanx,:)));
% % % %         bar(temp)
% % % %         colormap([0 .8 0])
% % % %         hold on
% % % %     end
% % % %     
% % % % 
% % % %     shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),1),(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_ALL'))))(chanx,:)),PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_err')))),colorarray(1),.5);
% % % %     hold on
% % % %     if NumConds>1
% % % %         for j=2:NumConds
% % % %             shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),1),(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_ALL'))))(chanx,:)),PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))),colorarray(2),.5);
% % % %         end
% % % %     end
% % % %     
% % % %     plot(zerobarERP)
% % % %     hold off
% % % %     title(strcat('Evoked CH:   ', chan_names(chanx)),'FontWeight','bold');
% % % %     set(gca, 'XTick', xtickrange)
% % % %     set(gca, 'XTickLabel', xticklabel)
% % % %     xlim([xaxis(1) xaxis(2)])
% % % %     if change_yaxis ==1
% % % %         ylim([yaxis_ERP(1) yaxis_ERP(2)])
% % % %     end
% % % %     legend(legend_array);
    
    pause
end

