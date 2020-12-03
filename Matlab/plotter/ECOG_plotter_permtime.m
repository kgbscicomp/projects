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
xtickrange = 1:srate*xticksize:size(hband_sig_buttord_wide,3);
colorarray = {'b','r','g'};

% create legend
legend_array = {};

for j=1:NumConds
    temp = find(CondNums==plot_conds(j));
    legend_array(length(legend_array)+1) = Conditions(temp);
end


% compare the first two, but still plot any others

PLOT.STATS_C = zeros(max(plot_chans),size(test_window,1));
PLOT.STATS_ERP = zeros(max(plot_chans),size(test_window,1));
PLOT.PermResultsERP = [];
PLOT.PermResultsC = [];

for chanx = plot_chans %[29 30 37 38 43:46 86:91]
    
    for j=1:NumConds
        
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))) = [];
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))) = [];
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_err')))) = [];
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))) = [];
        
        for i=1:numevents
            if Induced_chan_bad_array(chanx,i)==0
                if stim_resamp(i,2)==plot_conds(j)
                    temp = squeeze(hband_sig_buttord_wide(chanx,i,:));
                    tempsize = size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2)+1;
                    PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,tempsize) = temp;
                end
            end
            if Evoked_chan_bad_array(chanx,i)==0
                if stim_resamp(i,2)==plot_conds(j)
                    temp = squeeze(evoked_downsmpl_induced(chanx,i,:));
                    tempsize = size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)+1;
                    PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,tempsize) = temp;
                end
            end
        end
        
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_err')))) = std(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))')/sqrt(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2));
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))) = std(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))')/sqrt(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2));
        
        for i=1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2)
            PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(plotbaseline,i),1));
        end
        for i=1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)
            PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(plotbaseline,i),1));
        end
        
    end
    
    % calculate t-stat between first two condiitions for each timepoint
    
    test_smpl = [];
    for i=1:size(test_window,1)
        test_smpl = round(test_window(i,1)*srate):round(test_window(i,2)*srate);
        
        % permutation testing between 2 conditions
        % Evoked

        tempdataA = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))));
        tempdataB = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(2))))));
        PLOT.RealCombined = [mean(tempdataA(test_smpl,:),1),mean(tempdataB(test_smpl,:),1)];
        PLOT.RealDiff = mean2(tempdataA(test_smpl,:))-mean2(tempdataB(test_smpl,:));
        
        NumCondA = size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2);
        NumCondB = size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(2)))))),2);
        for zzz=1:NumPerms
            LabelRand = randperm(NumCondA+NumCondB);
            test1 = PLOT.RealCombined(:,LabelRand(1:NumCondA));
            test2 = PLOT.RealCombined(:,LabelRand(NumCondA+1:NumCondA+NumCondB));
            permdiff = mean(test1,2)-mean(test2,2);
            if abs(permdiff)>abs(PLOT.RealDiff)
                PLOT.STATS_ERP(chanx,i) = PLOT.STATS_ERP(chanx,i) +1;
            end
        end
        PLOT.STATS_ERP(chanx,i) = PLOT.STATS_ERP(chanx,i)/NumPerms;
        if PLOT.STATS_ERP(chanx,i)<pthresh
            PLOT.PermResultsERP(chanx,i) = 50;
        else
            PLOT.PermResultsERP(chanx,i) = 0;
        end
        
        % Induced
        tempdataA = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1))))));
        tempdataB = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(2))))));
        PLOT.RealCombined = [mean(tempdataA(test_smpl,:),1),mean(tempdataB(test_smpl,:),1)];
        PLOT.RealDiff = mean2(tempdataA(test_smpl,:))-mean2(tempdataB(test_smpl,:));
       
        NumCondA = size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2);
        NumCondB = size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(2)))))),2);
        for zzz=1:NumPerms
            LabelRand = randperm(NumCondA+NumCondB);
            test1 = PLOT.RealCombined(:,LabelRand(1:NumCondA));
            test2 = PLOT.RealCombined(:,LabelRand(NumCondA+1:NumCondA+NumCondB));
            permdiff = mean(test1,2)-mean(test2,2);
            if abs(permdiff)>abs(PLOT.RealDiff)
                PLOT.STATS_C(chanx,i) = PLOT.STATS_C(chanx,i) +1;
            end
        end
        PLOT.STATS_C(chanx,i) = PLOT.STATS_C(chanx,i)/NumPerms;
        if PLOT.STATS_C(chanx,i)<pthresh
            PLOT.PermResultsC(chanx,i) = 50;
        else
            PLOT.PermResultsC(chanx,i) = 0;
        end
    end
end

zerobarInduced = zeros(1,size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1));
zerobarERP = zeros(1,size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1));
zerobarInduced(zerotime) = min(mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2)); % 0 ms from aud onset
zerobarERP(zerotime) = min(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2))/2; % 0 ms from aud onset



for chanx = plot_chans %[29 30 37 38 43:46 86:91]
    
    figure(chanx+200)
    bar([PLOT.PermResultsHighInduced,PLOT.PermResultsLowInduced],'stacked')
    colormap([.8 .8 .8])
    hold on
    
    shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1),mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2)',PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)),'_err')))),colorarray(1),.5);
    hold on
    
    if NumConds>1
        for j=2:NumConds
            temp = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))));
            shadedErrorBar(1:size(temp,1),mean(temp,2)',PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_err')))),colorarray(j),.5); %{'color',[1,.647,0]},.5);
        end
    end
    
    plot(zerobarInduced)
    hold off
    title(strcat('Induced CH:   ', chan_names(chanx)),'FontWeight','bold');
    set(gca, 'XTick', xtickrange)
    set(gca, 'XTickLabel', xticklabel)
    xlim([xaxis(1) xaxis(2)])
    legend(legend_array);
    
    
    
    
    
    
    
    
    figure(chanx+100)
    bar([PLOT.PermResultsHighERP,PLOT.PermResultsLowERP],'stacked')
    colormap([.8 .8 .8])
    hold on
    
    shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),1),mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2)',PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_err')))),colorarray(1),.5);
    
    if NumConds>1
        for j=2:NumConds
            temp = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))));
            shadedErrorBar(1:size(temp,1),mean(temp,2)',PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))),colorarray(j),.5); %{'color',[1,.647,0]},.5);
        end
    end
    
    plot(zerobarERP)
    hold off
    title(strcat('Evoked CH:   ', chan_names(chanx)),'FontWeight','bold');
    set(gca, 'XTick', xtickrange)
    set(gca, 'XTickLabel', xticklabel)
    xlim([xaxis(1) xaxis(2)])
    if change_yaxis ==1
        ylim([yaxis_ERP(1) yaxis_ERP(2)])
    end
    legend(legend_array);
    
    % pause
end

