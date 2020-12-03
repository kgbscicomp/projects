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

xticklabel = SPECS.xaxis(1)-SPECS.zerotime:xticksize:SPECS.xaxis(2)-SPECS.zerotime;
xticklabel = round(xticklabel*10)/10;
temp = ceil((SPECS.xaxis(2)-SPECS.xaxis(1))*srate);
xtickrange = 1:srate*xticksize:temp;
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

STD.C1 = [];
STD.C2 = [];
STD.ERP1 = [];
STD.ERP2 = [];

for j=1:NumConds
    tempsizeC2=1;
    tempsizeERP2=1;
    for i=1:numevents
        
        AVG.tempC= [];
        AVG.tempERP = [];
        
        tempsizeC=1;
        tempsizeERP=1;
        tempChanx=1;
        for chanx = plot_chans %[29 30 37 38 43:46 86:91]
            if Induced_chan_bad_array(chanx,i)==0 && ~ismember(i,bad_trials)
                if stim_resamp(i,2)==plot_conds(j)
                    temp = squeeze(hband_data_all(chanx,:,[xaxis(1):xaxis(2)],i));
                    AVG.tempC(:,:,tempsizeC) = temp;
                    tempsizeC=tempsizeC+1;
                end
            end
            if Evoked_chan_bad_array(chanx,i)==0 && ~ismember(i,bad_trials)
                if stim_resamp(i,2)==plot_conds(j)
                    temp = (squeeze(evoked_sig(chanx,[xaxis(1):xaxis(2)],i)));
                    AVG.tempERP(:,tempsizeERP) = temp;
                    tempsizeERP = tempsizeERP+1;
                end
            end
            
        end
        if ~isempty(AVG.tempC)
            AVG.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,tempsizeC2) = squeeze(mean(mean(AVG.tempC,1),3))';
            tempsizeC2=tempsizeC2+1;
        end
        if ~isempty(AVG.tempERP)
            AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,tempsizeERP2) = squeeze((mean(AVG.tempERP,2)))';
            
            %baseline erp
            AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,tempsizeERP2) = AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(:,tempsizeERP2) -mean(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))(plotbaseline,tempsizeERP2));
            tempsizeERP2=1+tempsizeERP2;
        end
        
        
    end
end

% Cond switching perm tests
% ERP
PLOT_ERP.RealCombined = [AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(2))))))];
PLOT_ERP.RealDiff = mean(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2)-mean(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(2)))))),2);

PLOT_ERP.PermA = zeros(size(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1))))))));
PLOT_ERP.PermB = zeros(size(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(2))))))));
PLOT_ERP.PermResults = zeros(size(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),1),1);
PLOT_ERP.DiffResults = zeros(size(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),1),NumPerms);

NumCondA = size(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2);
NumCondB = size(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(2)))))),2);
for zzz=1:NumPerms
    LabelRand = randperm(NumCondA+NumCondB);
    test1 = PLOT_ERP.RealCombined(:,LabelRand(1:NumCondA));
    test2 = PLOT_ERP.RealCombined(:,LabelRand(NumCondA+1:NumCondA+NumCondB));
    permdiff = mean(test1,2)-mean(test2,2);
    PLOT_ERP.DiffResults(:,zzz)=permdiff;
    for i=1:size(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),1)
        if permdiff(i,1)>PLOT_ERP.RealDiff(i,1)
            PLOT_ERP.PermResults(i,1) = PLOT_ERP.PermResults(i,1) +1;
        end
    end
end
PLOT_ERP.PermResults = PLOT_ERP.PermResults./NumPerms;
PLOT_ERP.STDerr = std(PLOT_ERP.DiffResults')';

% Induced
PLOT_C.RealCombined = [AVG.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),AVG.(genvarname(char(strcat('C',num2str(plot_conds(2))))))];
PLOT_C.RealDiff = mean(AVG.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2)-mean(AVG.(genvarname(char(strcat('C',num2str(plot_conds(2)))))),2);

PLOT_C.PermA = zeros(size(AVG.(genvarname(char(strcat('C',num2str(plot_conds(1))))))));
PLOT_C.PermB = zeros(size(AVG.(genvarname(char(strcat('C',num2str(plot_conds(2))))))));
PLOT_C.PermResults = zeros(size(AVG.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1),1);
PLOT_C.DiffResults = zeros(size(AVG.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1),NumPerms);

NumCondA = size(AVG.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2);
NumCondB = size(AVG.(genvarname(char(strcat('C',num2str(plot_conds(2)))))),2);
for zzz=1:NumPerms
    LabelRand = randperm(NumCondA+NumCondB);
    test1 = PLOT_C.RealCombined(:,LabelRand(1:NumCondA));
    test2 = PLOT_C.RealCombined(:,LabelRand(NumCondA+1:NumCondA+NumCondB));
    permdiff = mean(test1,2)-mean(test2,2);
    PLOT_C.DiffResults(:,zzz)=permdiff;
    for i=1:size(AVG.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1)
        if permdiff(i,1)>PLOT_C.RealDiff(i,1)
            PLOT_C.PermResults(i,1) = PLOT_C.PermResults(i,1) +1;
        end
    end
end
PLOT_C.PermResults = PLOT_C.PermResults./NumPerms;
PLOT_C.STDerr = std(PLOT_C.DiffResults')';

for j=1:2
PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(j)))))) = mean(AVG.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2);
PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(j)))))) = mean(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2);
PLOT_C.(genvarname(char(strcat('std',num2str(plot_conds(j)))))) = std(AVG.(genvarname(char(strcat('C',num2str(plot_conds(j))))))')'/sqrt(size(AVG.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2));
PLOT_ERP.(genvarname(char(strcat('std',num2str(plot_conds(j)))))) = std(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))')'/sqrt(size(AVG.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),2));
end

zerobarInduced = zeros(1,size(PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),1));
zerobarERP = zeros(1,size(PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),1));
zerobarInduced(zerotime) = min([PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(1))))));PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(2))))))]); % 0 ms from aud onset
zerobarERP(zerotime) = min([PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(1))))));PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(2))))))]); % 0 ms from aud onset

figure(5)
h(1) = plot(PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),'Color',char(colorarray(1)));
hold on
shadedErrorBar(1:size(PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),1),PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),PLOT_C.STDerr,colorarray(1),.5);
h(2) = plot(PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(2)))))),'Color',char(colorarray(2)));
shadedErrorBar(1:size(PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(2)))))),1),PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(2)))))),PLOT_C.STDerr,colorarray(2),.5);
    
plot(zerobarInduced)
hold off
title('Induced Perm','FontWeight','bold');
set(gca, 'XTick', xtickrange)
set(gca, 'XTickLabel', xticklabel)
if change_yaxis ==1
    ylim([yaxis_C(1) yaxis_C(2)])
end
legend(h,legend_array,'Location','NorthWest');

    
    
figure(6)
h(1) = plot(PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),'Color',char(colorarray(1)));
hold on
shadedErrorBar(1:size(PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),1),PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),PLOT_ERP.STDerr,colorarray(1),.5);
h(2) = plot(PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(2)))))),'Color',char(colorarray(2)));
shadedErrorBar(1:size(PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(2)))))),1),PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(2)))))),PLOT_ERP.STDerr,colorarray(2),.5);

plot(zerobarERP)
hold off
title('Evoked Perm','FontWeight','bold');
set(gca, 'XTick', xtickrange)
set(gca, 'XTickLabel', xticklabel)
%     xlim([xaxis(1) xaxis(2)])
if change_yaxis ==1
    ylim([yaxis_ERP(1) yaxis_ERP(2)])
end
legend(h,legend_array,'Location','NorthWest');


figure(7)
h(1) = plot(PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),'Color',char(colorarray(1)));
hold on
shadedErrorBar(1:size(PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),1),PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),PLOT_C.(genvarname(char(strcat('std',num2str(plot_conds(1)))))),colorarray(1),.5);
h(2) = plot(PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(2)))))),'Color',char(colorarray(2)));
shadedErrorBar(1:size(PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(2)))))),1),PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(2)))))),PLOT_C.(genvarname(char(strcat('std',num2str(plot_conds(2)))))),colorarray(2),.5);
    
plot(zerobarInduced)
hold off
title('Induced std err','FontWeight','bold');
set(gca, 'XTick', xtickrange)
set(gca, 'XTickLabel', xticklabel)
xlim([1 size(PLOT_C.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),1)])
if change_yaxis ==1
    ylim([yaxis_C(1) yaxis_C(2)])
end
legend(h,legend_array,'Location','NorthWest');

    
    
figure(8)
h(1) = plot(PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),'Color',char(colorarray(1)));
hold on
shadedErrorBar(1:size(PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),1),PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),PLOT_ERP.(genvarname(char(strcat('std',num2str(plot_conds(1)))))),colorarray(1),.5);
h(2) = plot(PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(2)))))),'Color',char(colorarray(2)));
shadedErrorBar(1:size(PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(2)))))),1),PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(2)))))),PLOT_ERP.(genvarname(char(strcat('std',num2str(plot_conds(2)))))),colorarray(2),.5);

plot(zerobarERP)
hold off
title('Evoked std err','FontWeight','bold');
set(gca, 'XTick', xtickrange)
set(gca, 'XTickLabel', xticklabel)
xlim([1 size(PLOT_ERP.(genvarname(char(strcat('avg',num2str(plot_conds(1)))))),1)])
if change_yaxis ==1
    ylim([yaxis_ERP(1) yaxis_ERP(2)])
end
legend(h,legend_array,'Location','NorthWest');




