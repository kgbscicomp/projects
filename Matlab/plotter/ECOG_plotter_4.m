% ECOG_plotter

% Enter condition numbers interested in plotting

NumConds = length(plot_conds);

plotbaseline = [round((plotbaseline(1)*srate)+1):round((plotbaseline(2)*srate)+1)];
zbaseline = [round((zbaseline(1)*srate)+1):round((zbaseline(2)*srate)+1)];
zerotime = round((zerotime*srate)+1);
xaxis = [round((xaxis(1)*srate)+1) round((xaxis(2)*srate)+1)];
xticklabel = epoch(1):xticksize:epoch(2);
xtickrange = 1:srate*xticksize:size(hband_sig_morlet,3);
colorarray = {'r','g','b','y'};
if exist('bad_trials') == 0, bad_trials = [];end

% create legend
legend_array = {};

for j=1:NumConds
    temp = find(CondNums==plot_conds(j));
    legend_array(length(legend_array)+1) = Conditions(temp);
end

AVG.C1 = [];
AVG.C2 = [];
AVG.ERP1 = [];
AVG.ERP2 = [];

for chanx = plot_chans %[29 30 37 38 43:46 86:91]
     
    for j=1:NumConds
        
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))) = [];
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))) = [];
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_err')))) = [];
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))) = [];

        for i=1:numevents
            if Induced_chan_bad_array(chanx,i)==0 && ~ismember(i,bad_trials)
                if stim_resamp(i,2)==plot_conds(j)
                    temp = squeeze(hband_sig_morlet(chanx,i,:));
                    tempsize = size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2)+1;
                    PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,tempsize) = temp;
                end
            end
            if Evoked_chan_bad_array(chanx,i)==0 && ~ismember(i,bad_trials)
                if stim_resamp(i,2)==plot_conds(j)
                    temp = squeeze(evoked_sig(chanx,i,:));
                    tempsize = size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)+1;
                    PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,tempsize) = temp;                
                end
            end
        end
        
        for i=1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2)
            PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(plotbaseline,i),1));
        end
        for i=1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)
            PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(plotbaseline,i),1));
        end
        
        
    end
    
    % z-norm each trial then channel (avg/std calculated over all trials in tested conditions.
    alldata = [];
    for j=1:NumConds
    for i=1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2)
        alldata(:,size(alldata,2+1)) = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(zbaseline,i);
        meanx = mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(zbaseline,i));
        stdx = std(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(zbaseline,i));
        if z_on==1
            PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i) = (PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i)-meanx)/stdx;
        else
            PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i) = (PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i)-meanx);%/stdx;
        end
    end
    end
%     alldata = mean(alldata,2);
%     meanx = mean(alldata);
%     stdx = std(alldata);
    for j=1:NumConds
    for i=1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2)
%         PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i) = (PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i)-meanx)/stdx;
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(plotbaseline,i),1));
    end
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_err')))) = std(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))')/sqrt(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2));
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))) = std(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))')/sqrt(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2));
    end
    
    zerobarInduced = zeros(1,size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1));
    zerobarERP = zeros(1,size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1));
    zerobarInduced(zerotime) = min(mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2)); % 0 ms from aud onset
    zerobarERP(zerotime) = min(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2))/2; % 0 ms from aud onset
    
%     AVG.C1(:,chanx) = mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2);
%     AVG.ERP1(:,chanx) = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2);
%     
    AVG.(genvarname(char(strcat('C',num2str(plot_conds(1))))))(:,chanx) = mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2);
    AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))))(:,chanx) = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2);
            
    figure(5)
    shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1),mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2)',PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)),'_err')))),colorarray(1),.5);
    hold on
    

    
    if NumConds>1
        for j=2:NumConds
             
            temp = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))));
            shadedErrorBar(1:size(temp,1),mean(temp,2)',PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_err')))),colorarray(j),.5); %{'color',[1,.647,0]},.5);
            
            AVG.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,chanx) = mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2);
            AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,chanx) = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2);
        end
    end
    
    plot(zerobarInduced)
    hold off
    title(strcat('Induced CH:   ', chan_names(chanx)),'FontWeight','bold');
    set(gca, 'XTick', xtickrange)
    set(gca, 'XTickLabel', xticklabel)
    xlim([xaxis(1) xaxis(2)])
    if change_yaxis ==1
        ylim([yaxis_C(1) yaxis_C(2)])
    end
    legend(legend_array);
    
    figure(6) %chanx+100)
    shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),1),mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2)',PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_err')))),colorarray(1),.5);
    hold on
    
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
    
    pause
end

