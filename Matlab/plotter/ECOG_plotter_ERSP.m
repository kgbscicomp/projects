% ECOG_plotter

% Enter condition numbers interested in plotting

Freq = SPECS.Freq;
FreqWidth = SPECS.FreqWidth;
FreqCycles = SPECS.FreqCycles;
if length(FreqCycles)==1
    variablefreq = 1;
elseif length(FreqCycles)==2
    variablefreq = 0;
end
plotaxis=SPECS.plotaxis;
xaxis=SPECS.xaxis;
zerotime=SPECS.zerotime;
zbaseline=SPECS.zbaseline;
xticksize=SPECS.xticksize;
change_yaxis=SPECS.change_yaxis;
yaxis_ERP=SPECS.yaxis_ERP;
yaxis_C=SPECS.yaxis_C;
db_on=SPECS.db_on;
plot_conds=SPECS.plot_conds;
plot_chans=SPECS.plot_chans;

NumConds = length(plot_conds);
zerotime = zerotime-xaxis(1);
zbaseline = zbaseline;

zbaseline = [round((zbaseline(1)*srate)+1):round((zbaseline(2)*srate)+1)];
zerotime = round((zerotime*srate)+1);
xaxis = [round((xaxis(1)*srate)+1) round((xaxis(2)*srate)+1)];
xticklabel = epoch(1):xticksize:epoch(2);
xtickrange = 1:srate*xticksize:size(hband_data_all,3);
colorarray = {'r','g','b','y'};

times = (epoch(1)*srate:epoch(2)*srate)./srate;
frex = Freq(1):FreqWidth:Freq(2);

if exist('bad_trials') == 0, bad_trials = [];end

% create legend
legend_array = {};

for j=1:NumConds
    temp = find(CondNums==plot_conds(j));
    legend_array(length(legend_array)+1) = Conditions(temp);
end

clear PLOT AVG
for j=1:NumConds
    PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))=[];
    AVG.(genvarname(char(strcat('C',num2str(plot_conds(j)))))) = [];
%     PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(j))))))=[];
end
for chanx = plot_chans %[29 30 37 38 43:46 86:91]
     
    for j=1:NumConds
        PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))) = [];
        PLOT.(genvarname(char(strcat('avgERP',num2str(plot_conds(j)))))) = [];
        tempsizeERP=1;
        for i=1:numevents
            if Evoked_chan_bad_array(chanx,i)==0 && ~ismember(i,bad_trials)
                if stim_resamp(i,2)==plot_conds(j)
                    temp = (squeeze(evoked_sig(chanx,:,i)));
                    PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,tempsizeERP) = temp;
                    tempsizeERP=tempsizeERP+1;
                end
            end
            
        end
        
        for i=1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2)
            PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,i),1));
        end        
        
    end
    
    % freq analysis
    for j=1:NumConds
        temp = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))));
        [~,aa] = bbpower_phase(temp, srate, Freq,FreqWidth,FreqCycles, 0, srate, variablefreq);
        aa = aa.^2;
        for i=1:size(aa,2) % calculate power for each freq
            temp = squeeze(mean(aa(:,i,:),1));
            meanx = mean(temp(zbaseline));
            if db_on==1
                PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(i,:) = 10*log10(temp/meanx);
            elseif SPECS.z_on==1
                PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(i,:) = 10*log10(temp/meanx);
            else
                PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(i,:) = (temp-meanx);
            end
        end
        
        figure(10+j)
        contourf(times,frex,PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),40,'linecolor','none')
        tempmax = max(max(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,xaxis(1):xaxis(2))));
        tempmin = min(min(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,xaxis(1):xaxis(2))));
        templim = max([tempmax abs(tempmin)]);
%         set(gca,'clim',[-2 2],'xlim',plotaxis)
        set(gca,'xlim',plotaxis,'clim',[-templim templim])
        title(strcat('CH:   ', chan_names(chanx),      ' Cond: ',legend_array(j)),'FontWeight','bold');
        
        AVG.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(chanx,:,:) = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))));
    end
    pause();
end

