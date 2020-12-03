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
    
    % permutation testing between 2 conditions
    
    PLOT.RealCombined = [PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(2))))))];
    PLOT.RealDiff = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2)-mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(2)))))),2);
    
    PLOT.PermA = zeros(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))));
    PLOT.PermB = zeros(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))));
    PLOT.PermResults = zeros(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),1),1);
    
    PLOT.PermResultsLowERP = zeros(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),1),1);
    PLOT.PermResultsHighERP = zeros(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),1),1);
    
    % Evoked
    NumCondA = size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2);
    NumCondB = size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(2)))))),2);
    for zzz=1:NumPerms
        LabelRand = randperm(NumCondA+NumCondB);
        test1 = PLOT.RealCombined(:,LabelRand(1:NumCondA));
        test2 = PLOT.RealCombined(:,LabelRand(NumCondA+1:NumCondA+NumCondB));
        permdiff = mean(test1,2)-mean(test2,2);
        for i=1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),1)
           if permdiff(i,1)>PLOT.RealDiff(i,1)
               PLOT.PermResults(i,1) = PLOT.PermResults(i,1) +1;
           end
        end
    end
    PLOT.PermResults = PLOT.PermResults./NumPerms;
    temp = find(PLOT.PermResults<pthresh);
    PLOT.PermResultsLowERP(temp) = 50;
    temp = find(PLOT.PermResults>1-pthresh);
    PLOT.PermResultsHighERP(temp) = -50;
    

    % Induced
    
    PLOT.RealCombinedInduced = [PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),PLOT.(genvarname(char(strcat('C',num2str(plot_conds(2))))))];
    PLOT.RealDiffInduced = mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2)-mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(2)))))),2);
    PLOT.PermResultsLowInduced = zeros(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1),1);
    PLOT.PermResultsHighInduced = zeros(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1),1);
    
    PLOT.PermA = zeros(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))));
    PLOT.PermB = zeros(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))));
    PLOT.PermResults = zeros(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1),1);
    
    NumCondA = size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2);
    NumCondB = size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(2)))))),2);
    for zzz=1:NumPerms
        LabelRand = randperm(NumCondA+NumCondB);
        test1 = PLOT.RealCombinedInduced(:,LabelRand(1:NumCondA));
        test2 = PLOT.RealCombinedInduced(:,LabelRand(NumCondA+1:NumCondA+NumCondB));
        permdiff = mean(test1,2)-mean(test2,2);
        for i=1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1)
           if permdiff(i,1)>PLOT.RealDiffInduced(i,1)
               PLOT.PermResults(i,1) = PLOT.PermResults(i,1) +1;
           end
        end
    end
    PLOT.PermResults = PLOT.PermResults./NumPerms;
    temp = find(PLOT.PermResults<pthresh);
    PLOT.PermResultsLowInduced(temp) = 1;
    temp = find(PLOT.PermResults>1-pthresh);
    PLOT.PermResultsHighInduced(temp) = -.5;
    
    
    
    
    zerobarInduced = zeros(1,size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1));
    zerobarERP = zeros(1,size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1));
    zerobarInduced(zerotime) = min(mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2)); % 0 ms from aud onset
    zerobarERP(zerotime) = min(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2))/2; % 0 ms from aud onset
    
    
    
    
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

