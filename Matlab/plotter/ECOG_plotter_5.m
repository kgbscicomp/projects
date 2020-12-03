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
clear AVG
for chanx = plot_chans %[29 30 37 38 43:46 86:91]
     
    for j=1:NumConds
        
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))) = [];
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))) = [];
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_err')))) = [];
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))) = [];
        PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(j)))))) = [];
        PLOT.(genvarname(char(strcat('avgERP',num2str(plot_conds(j)))))) = [];
        tempsizeC=1;
        tempsizeERP=1;
        for i=1:numevents
            if Induced_chan_bad_array(chanx,i)==0 && ~ismember(i,bad_trials)
                if stim_resamp(i,2)==plot_conds(j)
                    temp = squeeze(hband_data_all(chanx,:,[xaxis(1):xaxis(2)],i));
                    PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,:,tempsizeC) = temp;
                    tempsizeC=tempsizeC+1;
                end
            end
            if Evoked_chan_bad_array(chanx,i)==0 && ~ismember(i,bad_trials)
                if stim_resamp(i,2)==plot_conds(j)
                    temp = (squeeze(evoked_sig(chanx,[xaxis(1):xaxis(2)],i)));
%                     tempsize = size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)+1;
%                     PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,tempsize) = temp; 
                    PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,tempsizeERP) = temp;
                    tempsizeERP=tempsizeERP+1;
                end
            end
            
        end
        
        for i=1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)
            PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(plotbaseline,i),1));
        end        
        
    end
    
    % Decible Power calculated on freq averages
    for j=1:NumConds
        for i=1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1) % calculate power for each freq
            tempdata = squeeze(mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(i,:,:),3));
            meanx = mean(tempdata(zbaseline));
            stdx = std(tempdata(zbaseline));
            if db_on==1
                PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(j))))))(:,i) = 10*log10(tempdata/meanx);
            else
                PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(j))))))(:,i) = (tempdata-meanx)/stdx;
            end
        end
        % run perms for stdev
        tempPerm = [];
        for k=1:NPerm
            permAvg = [];
            permX = randi(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),3),[1 size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),3)]);
            for i=1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1) % calculate power for each freq
                % pull the rand freqs many times then repeat the avging above
                permCurr = mean(squeeze(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(i,:,permX)),2);
                meanx = mean(permCurr(zbaseline));
                stdx = std(permCurr(zbaseline));
                if db_on==1
                    permAvg(:,i) = 10*log10(permCurr/meanx);
                else
                    permAvg(:,i) = (permCurr-meanx)/stdx;
                end
            end
            tempPerm(:,k) = mean(permAvg,2);
%             disp(['Perm: ' int2str(k) ' of ' int2str(NPerm)])
        end
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_err')))) = std(tempPerm');
    end
    for j=1:NumConds
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))) = std(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))')/sqrt(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2));
    end
    
    zerobarInduced = zeros(1,size(PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(1)))))),1));
    zerobarERP = zeros(1,size(PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(1)))))),1));
    zerobarInduced(zerotime) = min(mean(PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(1)))))),2)); % 0 ms from aud onset
    zerobarERP(zerotime) = min(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2))/2; % 0 ms from aud onset
    
%     AVG.C1(:,chanx) = mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2);
%     AVG.ERP1(:,chanx) = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2);
%     
    AVG.(genvarname(char(strcat('C',num2str(plot_conds(1))))))(:,chanx,:) = (PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(1)))))));
    AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))))(:,chanx) = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2);
            
    figure(5)
    h(1) = plot(mean(PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(1)))))),2)','Color',char(colorarray(1)));
    hold on
    shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(1)))))),1),mean(PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(1)))))),2)',PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)),'_err')))),colorarray(1),.5);
    

    
    if NumConds>1
        for j=2:NumConds
             
            temp = PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(j))))));
            h(j) = plot(mean(temp,2)','Color',char(colorarray(j)));
            shadedErrorBar(1:size(temp,1),mean(temp,2)',PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_err')))),colorarray(j),.5); %{'color',[1,.647,0]},.5);
            
            AVG.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,chanx,:) = (PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(j)))))));
            AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,chanx) = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2);
        end
    end
    
    plot(zerobarInduced)
    hold off
    title(strcat('Induced CH:   ', chan_names(chanx)),'FontWeight','bold');
    set(gca, 'XTick', xtickrange)
    set(gca, 'XTickLabel', xticklabel)
%     xlim([xaxis(1) xaxis(2)])
xlim([1 size(PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(1)))))),1)])
    if change_yaxis ==1
        ylim([yaxis_C(1) yaxis_C(2)])
    end
    legend(h,legend_array,'Location','NorthWest');
    
    
    
    figure(6) %chanx+100)
    h(1) = plot(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2)','Color',char(colorarray(1)));
    hold on
    shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),1),mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2)',PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_err')))),colorarray(1),.5);
    
    
    if NumConds>1
        for j=2:NumConds
            temp = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))));
            h(j) = plot(mean(temp,2)','Color',char(colorarray(j)));
            shadedErrorBar(1:size(temp,1),mean(temp,2)',PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))),colorarray(j),.5); %{'color',[1,.647,0]},.5);
        end
    end
  
    plot(zerobarERP)
    hold off
    title(strcat('Evoked CH:   ', chan_names(chanx)),'FontWeight','bold');
    set(gca, 'XTick', xtickrange)
    set(gca, 'XTickLabel', xticklabel)
%     xlim([xaxis(1) xaxis(2)])
    xlim([1 size(PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(1)))))),1)])
    if change_yaxis ==1
        ylim([yaxis_ERP(1) yaxis_ERP(2)])
    end
    legend(h,legend_array,'Location','NorthWest');
    
    pause
end

