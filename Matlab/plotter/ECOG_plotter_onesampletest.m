% ECOG_plotter

% Enter condition numbers interested in plotting

NumConds = length(plot_conds);

plotbaseline = [round((plotbaseline(1)*srate)+1):round((plotbaseline(2)*srate)+1)];
zerotime = round((zerotime*srate)+1);
xaxis = [round((xaxis(1)*srate)+1) round((xaxis(2)*srate)+1)];
xticklabel = epoch(1):xticksize:epoch(2);
xtickrange = 1:srate*xticksize:size(hband_sig_buttord_wide,3);
colorarray = {'r','g','b','y'};

% create legend
legend_array = {};

for j=1:NumConds
    temp = find(CondNums==plot_conds(j));
    legend_array(length(legend_array)+1) = Conditions(temp);
end

for chanx = plot_chans %[29 30 37 38 43:46 86:91]
    
    OneSmplStatC = [];
    OneSmplStatC_thresh = [];
    OneSmplStatERP = [];
    OneSmplStatERP_thresh = [];
    
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
        
        
        for i=1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2)
            PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(plotbaseline,i),1));
        end
        for i=1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)
            PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i) = (PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(plotbaseline,i),1)))/squeeze(std(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(plotbaseline,i),1));
        end
        
        PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_err')))) = std(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))')/sqrt(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2));
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))) = std(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))')/sqrt(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2));
        
        % one sample ttests for each condition and separately for Induced
        % and Evoked Conditions.
        for i=1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1)
            tempval = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(i,:);
            [h,p]=ttest(tempval);
            OneSmplStatC(i,j) = p;
            tempval = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(i,:);
            [h,p]=ttest(tempval);
            OneSmplStatERP(i,j) = p;
        end
        OneSmplStatC_thresh(:,j) = zeros(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1),1);
        OneSmplStatERP_thresh(:,j) = zeros(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1),1);
        
        for i=zerotime:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1)
            if OneSmplStatC(i,j) <= pthresh
                OneSmplStatC_thresh(i,j) = .25+((j-1)*.25);
            end 
            if OneSmplStatERP(i,j) <= pthresh
                OneSmplStatERP_thresh(i,j) = .25+((j-1)*.25);
            end 
        end
    end
    
    zerobarInduced = zeros(1,size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1));
    zerobarERP = zeros(1,size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1));
    zerobarInduced(zerotime) = min(mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2)); % 0 ms from aud onset
    zerobarERP(zerotime) = min(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2))/2; % 0 ms from aud onset
    
    figure(5)
    
    bar(OneSmplStatC_thresh)
    hold on
    tempsize=size(OneSmplStatC_thresh,2);
    if tempsize == 1, colormap([1 .4 .4]);
    elseif tempsize == 2, colormap([1 .4 .4;.4 1 .4]);
    elseif tempsize == 3, colormap([1 .4 .4;.4 1 .4;.4 .4 1]);
    else colormap([1 .4 .4;.4 1 .4;.4 .4 1;1 1 .4]);end        
    shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1),mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2)',PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)),'_err')))),colorarray(1),.5);
    
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
    
    figure(6) %chanx+100)
    bar(OneSmplStatERP_thresh)
    hold on
    tempsize=size(OneSmplStatERP_thresh,2);
    if tempsize == 1, colormap([1 .4 .4]);
    elseif tempsize == 2, colormap([1 .4 .4;.4 1 .4]);
    elseif tempsize == 3, colormap([1 .4 .4;.4 1 .4;.4 .4 1]);
    else colormap([1 .4 .4;.4 1 .4;.4 .4 1;1 1 .4]);end   
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
    
    pause
end

