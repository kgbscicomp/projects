% ECOG_plotter

% Enter condition numbers interested in plotting
xaxis=SPECS.xaxis;
zerotime=SPECS.zerotime;
plotbaseline=SPECS.plotbaseline;
zbaseline=SPECS.zbaseline;
xticksize=SPECS.xticksize;
change_yaxis=SPECS.change_yaxis;
yaxis_ERP=SPECS.yaxis_ERP;
yaxis_C=SPECS.yaxis_C;
db_on=SPECS.db_on;
plot_conds=SPECS.plot_conds;
plot_chans=SPECS.plot_chans;
pthresh=SPECS.pthresh;
NumPerms=SPECS.NumPerms;
NPerm=SPECS.NPerm;



NumConds = length(plot_conds);
zerotime = zerotime-xaxis(1);
plotbaseline = plotbaseline-xaxis(1);
zbaseline = zbaseline-xaxis(1);

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

PLOT.STATS_ERP = zeros(max(plot_chans),round((xaxis(2)-xaxis(1)+1)));
PLOT.STATS_ERP_thresh = [];

PLOT.PermResultsERP = [];
PLOT.PermResultsERPmax = [];
PLOT.PermResultsERPmin = [];
PLOT.PermThreshERP = [];

for j=1:NumConds
    PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_ALL')))) = [];
    PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_ALL_err')))) = [];
end
for chanx = plot_chans %[29 30 37 38 43:46 86:91]
    
    for j=1:NumConds
        
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))) = [];
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))) = [];
        
        for i=1:numevents
            if Evoked_chan_bad_array(chanx,i)==0
                if stim_resamp(i,2)==plot_conds(j)
                    temp = (squeeze(evoked_sig(chanx,[xaxis(1):xaxis(2)],i)));
                    tempsize = size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)+1;
                    PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,tempsize) = temp;                
                end
            end
        end
        
        for i=1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)
            PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(plotbaseline,i),1));
        end 
        
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))) = std(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))')/sqrt(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2));        
        for i=1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)
            PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(plotbaseline,i),1));
        end
        
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_ALL'))))(chanx,:) = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)';
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_ALL_err'))))(chanx,:) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err'))));        
    end
    
    % time varying null distribution and pixel based multiple comp corr
    % 1. Save M trial by X time array
    % 2. For N perms, for M trials reorder the timeseries at 1 pt (define cut pt and move data after cut to start of timeseries)
    % 3. For N perms, save the max and min (2 vals) during the epoch
    % 4. Find lower and upper 2.5% of these max/min distributions
    % 5. Any real ERP value above the max or below the min is sig
    
    % Perm
    test_smpl = xaxis(1):xaxis(2);
    if pthresh==-1, NumPerms=1000;end
    for zzz=1:NumPerms
        Smpl_Order = [];
        for yyy=1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2)
            randX = randi(length(test_smpl));
            Smpl_Order(:,yyy) = [test_smpl(randX:end) test_smpl(1:randX-1)];
        end
        tempAVG = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))))(Smpl_Order);
        tempAVG = mean(tempAVG,2);
        tempAVG = tempAVG-mean(tempAVG(plotbaseline,1));
        PLOT.PermResultsERPmax(zzz,chanx) = max(tempAVG(zerotime-xaxis(1)+1:xaxis(2)-xaxis(1)+1,1));
        PLOT.PermResultsERPmin(zzz,chanx) = min(tempAVG(zerotime-xaxis(1)+1:xaxis(2)-xaxis(1)+1,1));
    end
    if pthresh==-1
        CutoffValmax = [quantile(1:1000,1-.05)+.5 quantile(1:1000,1-.01)+.5 quantile(1:1000,1-.005)+.5];
        CutoffValmin = [quantile(1:1000,.05)-.5 quantile(1:1000,.01)-.5 quantile(1:1000,.005)-.5];
        tempA = sort(PLOT.PermResultsERPmax(:,chanx));
        tempA = tempA(CutoffValmax');
        tempB = sort(PLOT.PermResultsERPmin(:,chanx));
        tempB = tempB(CutoffValmin');
        for i=1:3
            tempMAX = find(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_ALL'))))(chanx,zerotime-xaxis(1)+1:xaxis(2)-xaxis(1)+1)>tempA(i),1);
            tempMIN = find(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_ALL'))))(chanx,zerotime-xaxis(1)+1:xaxis(2)-xaxis(1)+1)<tempB(i),1);
            if isempty(tempMIN) && isempty(tempMAX)
                PLOT.STATS_ERP_thresh(chanx,i) = 0;
            else
                PLOT.STATS_ERP_thresh(chanx,i) = 1;
            end
        end
        else
        CutoffVal = quantile(1:NumPerms,1-pthresh)+.5;
        temp = sort(PLOT.PermResultsERPmax(:,chanx));
        temp = temp(CutoffVal);
        tempMAX = find(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_ALL'))))(chanx,zerotime-xaxis(1)+1:xaxis(2)-xaxis(1)+1)>temp,1);

        CutoffVal = quantile(1:NumPerms,pthresh)-.5;
        temp = sort(PLOT.PermResultsERPmin(:,chanx));
        temp = temp(CutoffVal);
        tempMIN = find(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_ALL'))))(chanx,zerotime-xaxis(1)+1:xaxis(2)-xaxis(1)+1)<temp,1);

        if isempty(tempMIN) && isempty(tempMAX)
            PLOT.STATS_ERP_thresh(chanx,1) = 0;
        else
            PLOT.STATS_ERP_thresh(chanx,1) = 1;
        end
    end
     disp(['Channel Number ' num2str(chanx)])
end
    

% % % % % 
% % % % % 
% % % % % 
% % % % %     % using cluster thresholds by testing tmass of epoch against
% % % % %     % 0 w one-sample ttest. then building null dist of tmss from baseline w
% % % % %     % one sample ttests.
% % % % %     
% % % % %     
% % % % %     % calculate one sample t-stat for each timepoint of cond 1
% % % % %     
% % % % %     [H,P,CI,STATS] = ttest(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))))');
% % % % %     RealT_ERP = STATS.tstat;
% % % % %     RealP_ERP = P;
% % % % %     
% % % % %     tmax_ERP = zeros(1,size(RealT_ERP,2));
% % % % %     temp_ERP = 0;
% % % % %     counter_ERP = 1;
% % % % %     
% % % % %     % calc real tmass
% % % % %     for i=zerotime+1:size(RealT_ERP,2)
% % % % %         if ((RealT_ERP(i)<0 && RealT_ERP(i-1)>0)||(RealT_ERP(i)>0 && RealT_ERP(i-1)<0))
% % % % %             temp_ERP = 0;
% % % % %             counter_ERP = 1;
% % % % %         end
% % % % %         if RealP_ERP(i)<=pthresh
% % % % %             temp_ERP = temp_ERP+RealT_ERP(i);
% % % % %             counter_ERP = counter_ERP-1;
% % % % %         else
% % % % %             temp_ERP = 0;
% % % % %             counter_ERP = 1;
% % % % %         end
% % % % %         tmax_ERP(i) = temp_ERP;
% % % % %         if counter_ERP < 0
% % % % %            for j=counter_ERP:-1
% % % % %                tmax_ERP(i+j) = temp_ERP;
% % % % %            end
% % % % %         end
% % % % %     end
% % % % %     PLOT.STATS_ERP(chanx,:) = tmax_ERP;
% % % % %     
% % % % %     % Perm
% % % % %     test_smpl = 1:zerotime;
% % % % %     for zzz=1:NumPerms
% % % % %         
% % % % %         Smpl_Order = [];
% % % % %         for yyy=1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2)
% % % % %             randX = randi(length(test_smpl));
% % % % %             Smpl_Order(:,yyy) = [test_smpl(randX:end) test_smpl(1:randX-1)];
% % % % %         end
% % % % %         
% % % % %         [H,P,CI,STATS] = ttest(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))))(Smpl_Order)');
% % % % %         RealT_ERP = STATS.tstat;
% % % % %         RealP_ERP = P;
% % % % %         
% % % % %         tmax_ERP = zeros(1,size(RealT_ERP,2));
% % % % %         temp_ERP = 0;
% % % % %         counter_ERP = 1;
% % % % %         
% % % % %         % calc perm tmass
% % % % %         for i=1:zerotime
% % % % %             if i>1
% % % % %                 if ((RealT_ERP(i)<0 && RealT_ERP(i-1)>0)||(RealT_ERP(i)>0 && RealT_ERP(i-1)<0))
% % % % %                     temp_ERP = 0;
% % % % %                     counter_ERP = 1;
% % % % %                 end
% % % % %             end
% % % % %             if RealP_ERP(i)<=pthresh
% % % % %                 temp_ERP = temp_ERP+RealT_ERP(i);
% % % % %                 counter_ERP = counter_ERP-1;
% % % % %             else
% % % % %                 temp_ERP = 0;
% % % % %                 counter_ERP = 1;
% % % % %             end
% % % % %             tmax_ERP(i) = temp_ERP;
% % % % %             if counter_ERP < 0
% % % % %                 for j=counter_ERP:-1
% % % % %                     tmax_ERP(i+j) = temp_ERP;
% % % % %                 end
% % % % %             end
% % % % %         end
% % % % %         PLOT.PermResultsERP(chanx,zzz) = max(abs(tmax_ERP));
% % % % %         
% % % % %     end
% % % % %     
% % % % %     CutoffVal = quantile(1:NumPerms,1-pthresh)+.5;
% % % % %     temp = sort(PLOT.PermResultsERP(chanx,:));
% % % % %     PLOT.PermThreshERP(chanx,1) = temp(CutoffVal);
% % % % %     temp = find(abs(PLOT.STATS_ERP(chanx,:))>PLOT.PermThreshERP(chanx,1),1);
% % % % %     if ~isempty(temp)
% % % % %         PLOT.STATS_ERP_thresh(chanx,1) = temp;
% % % % %     end
% % % % %     disp(['Channel Number ' num2str(chanx)])
% % % % %     
% % % % % end

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
%     xlim([xaxis(1) xaxis(2)])
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
%     xlim([xaxis(1) xaxis(2)])
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

