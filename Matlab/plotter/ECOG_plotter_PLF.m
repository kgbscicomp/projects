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

% first run the calculations then plot in 2nd for loop

CondA_Z_ALL = [];
CondB_Z_ALL = [];
CondA_ITC_ALL = [];
CondB_ITC_ALL = [];
PermZ_ALL_Chan = [];
PermITC_All_Chan = [];
PermZ_ALL_Chan_thresh = [];
PermITC_ALL_Chan_thresh = [];
CurrCondTrialsA = find(stim_resamp(:,2)==CondNums(1)); %pred
CurrCondTrialsB = find(stim_resamp(:,2)==CondNums(2)); %nonpred

for chanx = plot_chans %[29 30 37 38 43:46 86:91]
    disp(['Channel ' int2str(chanx) ' of ' int2str(size(plot_chans,2))])
    
    tempPredZ = [];
    tempPredZ_Std = [];
    tempNonPredZ= [];
    tempNonPredZ_Std= [];
    tempPredITC = [];
    tempNonPredITC=[];
    for j=1:trial_length
            [pval,zval] = circ_rtest(angle(squeeze(hband_sig_buttord_wide(CurrCondTrialsA,j,1,chanx))));
            tempPredZ(:,j,1) = zval;
            [pval,zval] = circ_rtest(angle(squeeze(hband_sig_buttord_wide(CurrCondTrialsB,j,1,chanx))));
            tempNonPredZ(:,j,1) = zval;
            tempPredITC(:,j,1) = squeeze(hband_sig_buttord_wide(CurrCondTrialsA',j,1,chanx));
            tempNonPredITC(:,j,1) = squeeze(hband_sig_buttord_wide(CurrCondTrialsB',j,1,chanx));
    end
    CondA_Z_ALL(chanx,:)=tempPredZ;
    CondB_Z_ALL(chanx,:)=tempNonPredZ;
    
    %     ITC, trials x time x freq
    ITCnonPred = squeeze(abs(mean(tempNonPredITC./abs(tempNonPredITC))));
    ITCPred = squeeze(abs(mean(tempPredITC./abs(tempPredITC))));
    
    CondA_ITC_ALL(chanx,:)=ITCPred;
    CondB_ITC_ALL(chanx,:)=ITCnonPred;
    
  
% % % % %     tempPredZ = [];
% % % % %     tempPredZ_Std = [];
% % % % %     tempNonPredZ= [];
% % % % %     tempNonPredZ_Std= [];
% % % % %     tempPredITC = [];
% % % % %     tempNonPredITC=[];
% % % % %     for j=1:trial_length
% % % % %         for k=chanX %1:num_el
% % % % %             for m=1:size(hband_sig_buttord_wide,3)
% % % % %                 [pval,zval] = circ_rtest(angle(squeeze(hband_sig_elliptic_wide(CurrCondTrialsA,j,m,k))));
% % % % %                 tempPredZ(:,j,m) = zval;
% % % % %                 tempPredZ_Std(:,j,m) = circ_std(angle(squeeze(hband_sig_elliptic_wide(CurrCondTrialsA,j,m,k))));
% % % % %                 [pval,zval] = circ_rtest(angle(squeeze(hband_sig_elliptic_wide(CurrCondTrialsB,j,m,k))));
% % % % %                 tempNonPredZ(:,j,m) = zval;
% % % % %                 tempNonPredZ_Std(:,j,m) = circ_std(angle(squeeze(hband_sig_elliptic_wide(CurrCondTrialsB,j,m,k))));
% % % % %                 tempPredITC(:,j,m) = squeeze(hband_sig_elliptic_wide(CurrCondTrialsA',j,m,k));
% % % % %                 tempNonPredITC(:,j,m) = squeeze(hband_sig_elliptic_wide(CurrCondTrialsB',j,m,k));
% % % % %             end
% % % % %         end
% % % % %     end
% % % % %     tempPredZ_Std = tempPredZ_Std/sqrt(32);
% % % % %     tempNonPredZ_Std = tempNonPredZ_Std/sqrt(32);
    
% % % % %     figure(7)
% % % % %     shadedErrorBar(1:size(tempPredZ,2),tempPredZ,tempPredZ_Std,'b',.5);
% % % % %     hold on
% % % % %     shadedErrorBar(1:size(tempNonPredZ,2),tempNonPredZ,tempNonPredZ_Std,'r',.5);
% % % % %     hold off
% % % % %     title(strcat(char(chan_names(chanstoplot(zzz))),' Elliptic'));
% % % % %     
    
    
%     ITC, trials x time x freq
% % % % %     ITCnonPred = squeeze(abs(mean(tempNonPredITC./abs(tempNonPredITC))));
% % % % %     ITCPred = squeeze(abs(mean(tempPredITC./abs(tempPredITC))));
% % % % %     figure(8)
% % % % %     plot([ITCPred',ITCnonPred']);figure(gcf);ylim([0 1])
% % % % %     title(strcat(char(chan_names(chanstoplot(zzz))),' ITC Elliptic'));
% % % % %     
   


    
    permdataZ = zeros(1,size(ITCnonPred,2));
    permdataITC= zeros(1,size(ITCnonPred,2));
    realdataITC = [tempNonPredITC;tempPredITC];
    realdataZ = [];
    realdataZPred = [];
    realdataZNonPred = [];
    for i=1:size(ITCnonPred,2)
        realdataZ(:,i) = [angle(squeeze(hband_sig_buttord_wide(CurrCondTrialsA,i,1,chanX)));angle(squeeze(hband_sig_buttord_wide(CurrCondTrialsB,i,1,chanX)))];
        realdataZPred(:,i) = [angle(squeeze(hband_sig_buttord_wide(CurrCondTrialsA,i,1,chanX)))];
        realdataZNonPred(:,i) = [angle(squeeze(hband_sig_buttord_wide(CurrCondTrialsB,i,1,chanX)))];
    end
    
    real_PredZ =tempPredZ;
    real_NonPredZ = tempNonPredZ;
    real_DiffZ = real_PredZ - real_NonPredZ;
    
    real_PredITC =ITCPred;
    real_NonPredITC = ITCnonPred;
    real_DiffITC = real_PredITC - real_NonPredITC;
    
        
    
    PermZ_ALL = [];
    CircR_Pred_All = [];
    CircR_NonPred_All = [];
    permtestPredZ = [];
    permtestNonPredZ = [];
    for k=1:NumPerms % perm number
        randY = randperm(size(CurrCondTrialsA,1)+size(CurrCondTrialsB,1));
        
        randZ = realdataZ(randY,:); %randi(64,1,32);
        test1Z = randZ(1:size(CurrCondTrialsA,1),:); % 24 trials
        test2Z = randZ(size(CurrCondTrialsA,1)+1:size(CurrCondTrialsA,1)+size(CurrCondTrialsB,1),:); % 22 trials
        
        randITC = realdataITC(randY,:); %randi(64,1,32);
        test1ITC = randITC(1:size(CurrCondTrialsA,1),:);
        test2ITC = randITC(size(CurrCondTrialsA,1)+1:size(CurrCondTrialsA,1)+size(CurrCondTrialsB,1),:);
        
        ITC_test1 = squeeze(abs(mean(test1ITC./abs(test1ITC))));
        ITC_test2 = squeeze(abs(mean(test2ITC./abs(test2ITC))));
        DiffITC = ITC_test1-ITC_test2;
        
        for j=1:size(ITCnonPred,2)
            [pval,zval] = circ_rtest((squeeze(test1Z(:,j))));
            CircR_Pred = zval;
            [pval,zval] = circ_rtest((squeeze(test2Z(:,j))));
            CircR_NonPred = zval;
            DiffZ = CircR_Pred-CircR_NonPred;
            
            CircR_Pred_All(k,j) = CircR_Pred;
            CircR_NonPred_All(k,j) = CircR_NonPred;
            DiffZ_ALL(j) = DiffZ;
            
            temp = realdataZPred(:,j)';
            permtestPred = temp(ceil(length(temp) * rand(length(temp),1)));
            [pval,zval] = circ_rtest(permtestPred);
            permtestPredZ(k,j) = zval;
            
            temp = realdataZNonPred(:,j)';
            permtestNonPred = temp(ceil(length(temp) * rand(length(temp),1)));
            [pval,zval] = circ_rtest(permtestNonPred);
            permtestNonPredZ(k,j) = zval;
            
            if real_DiffZ(1,j)<=DiffZ
                permdataZ(1,j) =  permdataZ(1,j)+1;
            end
            if real_DiffITC(1,j)<=DiffITC(1,j)
                permdataITC(1,j) =  permdataITC(1,j)+1;
            end
        end
% % % % %         
% % % % %         % downsample the diff data
% % % % %         DiffZ_rs = resample(DiffZ_ALL,100,srate);
% % % % %         DiffITC_rs = resample(DiffITC,100,srate);
% % % % %         for j=1:151
% % % % %             if real_DiffZ_rs(1,j)<=DiffZ_rs(1,j)
% % % % %                 permdataZ_rs(j,1) =  permdataZ_rs(j,1)+1;
% % % % %             end
% % % % %             if real_DiffITC_rs(1,j)<=DiffITC_rs(1,j)
% % % % %                 permdataITC_rs(j,1) =  permdataITC_rs(j,1)+1;
% % % % %             end
% % % % %             if real_DiffZ_rs(1,j)>=DiffZ_rs(1,j)
% % % % %                 permdataZ_rs(j,2) =  permdataZ_rs(j,2)+1;
% % % % %             end
% % % % %             if real_DiffITC_rs(1,j)>=DiffITC_rs(1,j)
% % % % %                 permdataITC_rs(j,2) =  permdataITC_rs(j,2)+1;
% % % % %             end
% % % % %         end
% % % % %         
        PermZ_ALL(k,:) = DiffZ_ALL;
        
        
    end
    
    permdataZ = permdataZ/NumPerms;
    permdataITC = permdataITC/NumPerms;
% % % % %     permdataZ_rs = permdataZ_rs/NumPerms;%figure(9),plot(permdataZ_rs(:,1))
% % % % %     permdataITC_rs = permdataITC_rs/NumPerms;%figure(10),plot(permdataITC_rs(:,1))
    
    
% % % % %     figure(zzz+10)
% % % % %     shadedErrorBar(1:size(tempPredZ,2),tempPredZ,permtestPredZ_Std,'b',1);
% % % % %     hold on
% % % % %     shadedErrorBar(1:size(tempNonPredZ,2),tempNonPredZ,permtestNonPredZ_Std,'r',1);
% % % % %     hold off
% % % % %     title(char(chan_names(chanstoplot(zzz))));
    
% % % % %     bar([permdataZ,PLOT.PermResultsLowInduced],'stacked')
% % % % %     colormap([.8 .8 .8])
% % % % %     hold on
% % % % % 


    PermZ_ALL_Chan(chanx,:)=permdataZ;
    PermITC_ALL_Chan(chanx,:)=permdataITC;
    PermZ_ALL_Chan_thresh(chanx,:) = zeros(1,size(permdataZ,2));
    PermITC_ALL_Chan_thresh(chanx,:) = zeros(1,size(permdataZ,2));
    for jjj=1:size(permdataZ,2)
        if PermZ_ALL_Chan(chanx,jjj)<=.1
            PermZ_ALL_Chan_thresh(chanx,jjj) =  PermZ_ALL_Chan(chanx,jjj);
        end
        if PermITC_ALL_Chan(chanx,jjj)<=.1
            PermITC_ALL_Chan_thresh(chanx,jjj) =  PermITC_ALL_Chan(chanx,jjj);
        end
    end
end


for chanx = plot_chans %[29 30 37 38 43:46 86:91]
    
    figure(5)
    bar([PermZ_ALL_Chan(chanx,:)',PermZ_ALL_Chan_thresh(chanx,:)'])
    hold on
    colormap([.95 .95 .95;.8 0 0])
    plot([CondA_Z_ALL(chanx,:)',CondB_Z_ALL(chanx,:)']);figure(gcf);
    plot(zeros(size(CondA_Z_ALL(chanx,:),2),1)+.05);
    title(strcat(char(chan_names(chanstoplot(chanx))),' Z'));
    hold off
    
    figure(6)
    bar([PermITC_ALL_Chan(chanx,:)',PermITC_ALL_Chan_thresh(chanx,:)'])
    hold on
    colormap([.95 .95 .95;.8 0 0])
    plot([CondA_ITC_ALL(chanx,:)',CondB_ITC_ALL(chanx,:)']);figure(gcf);
    plot(zeros(size(CondA_ITC_ALL(chanx,:),2),1)+.05);
    title(strcat(char(chan_names(chanstoplot(chanx))),' ITC'));
    hold off
    
    pause()
end


% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % %         
% % %             for j=1:NumConds
% % % 
% % %         
% % %         PLOT.(genvarname(char(strcat('BUTT',num2str(plot_conds(j)))))) = [];
% % %         PLOT.(genvarname(char(strcat('ELLIP',num2str(plot_conds(j)))))) = [];
% % %         
% % %         for i=1:numevents
% % %             if Induced_chan_bad_array(chanx,i)==0
% % %                 if stim_resamp(i,2)==plot_conds(j)
% % %                     temp = squeeze(hband_sig_buttord_wide(chanx,i,:));
% % %                     tempsize = size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),2)+1;
% % %                     PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,tempsize) = temp;
% % %                 end
% % %             end
% % %         end
% % %         
% % %         for i=1:size(PLOT.(genvarname(char(strcat('BUTT',num2str(plot_conds(j)))))),2)
% % %             PLOT.(genvarname(char(strcat('BUTT',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))(plotbaseline,i),1));
% % %         end
% % %         for i=1:size(PLOT.(genvarname(char(strcat('ELLIP',num2str(plot_conds(j)))))),2)
% % %             PLOT.(genvarname(char(strcat('ELLIP',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ELLIP',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ELLIP',num2str(plot_conds(j))))))(plotbaseline,i),1));
% % %         end
% % %         
% % %     end
% % %     
% % %     % permutation testing between 2 conditions
% % %     
% % %     PLOT.RealCombined = [PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(2))))))];
% % %     PLOT.RealDiff = mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2)-mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(2)))))),2);
% % %     
% % %     PLOT.PermA = zeros(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))));
% % %     PLOT.PermB = zeros(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))))));
% % %     PLOT.PermResults = zeros(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),1),1);
% % %     
% % %     PLOT.PermResultsLowERP = zeros(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),1),1);
% % %     PLOT.PermResultsHighERP = zeros(size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))),1),1);
% % %     
% % %     % Evoked
% % %     NumCondA = size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2);
% % %     NumCondB = size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(2)))))),2);
% % %     for zzz=1:NumPerms
% % %         LabelRand = randperm(NumCondA+NumCondB);
% % %         test1 = PLOT.RealCombined(:,LabelRand(1:NumCondA));
% % %         test2 = PLOT.RealCombined(:,LabelRand(NumCondA+1:NumCondA+NumCondB));
% % %         permdiff = mean(test1,2)-mean(test2,2);
% % %         for i=1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),1)
% % %            if permdiff(i,1)>PLOT.RealDiff(i,1)
% % %                PLOT.PermResults(i,1) = PLOT.PermResults(i,1) +1;
% % %            end
% % %         end
% % %     end
% % %     PLOT.PermResults = PLOT.PermResults./NumPerms;
% % %     temp = find(PLOT.PermResults<pthresh);
% % %     PLOT.PermResultsLowERP(temp) = 50;
% % %     temp = find(PLOT.PermResults>1-pthresh);
% % %     PLOT.PermResultsHighERP(temp) = -50;
% % %     
% % % 
% % %     % Induced
% % %     
% % %     PLOT.RealCombinedInduced = [PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),PLOT.(genvarname(char(strcat('C',num2str(plot_conds(2))))))];
% % %     PLOT.RealDiffInduced = mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2)-mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(2)))))),2);
% % %     PLOT.PermResultsLowInduced = zeros(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1),1);
% % %     PLOT.PermResultsHighInduced = zeros(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1),1);
% % %     
% % %     PLOT.PermA = zeros(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))));
% % %     PLOT.PermB = zeros(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))));
% % %     PLOT.PermResults = zeros(size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)))))),1),1);
% % %     
% % %     NumCondA = size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2);
% % %     NumCondB = size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(2)))))),2);
% % %     for zzz=1:NumPerms
% % %         LabelRand = randperm(NumCondA+NumCondB);
% % %         test1 = PLOT.RealCombinedInduced(:,LabelRand(1:NumCondA));
% % %         test2 = PLOT.RealCombinedInduced(:,LabelRand(NumCondA+1:NumCondA+NumCondB));
% % %         permdiff = mean(test1,2)-mean(test2,2);
% % %         for i=1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1)
% % %            if permdiff(i,1)>PLOT.RealDiffInduced(i,1)
% % %                PLOT.PermResults(i,1) = PLOT.PermResults(i,1) +1;
% % %            end
% % %         end
% % %     end
% % %     PLOT.PermResults = PLOT.PermResults./NumPerms;
% % %     temp = find(PLOT.PermResults<pthresh);
% % %     PLOT.PermResultsLowInduced(temp) = 1;
% % %     temp = find(PLOT.PermResults>1-pthresh);
% % %     PLOT.PermResultsHighInduced(temp) = -.5;
% % %     
% % %     
% % %     
% % %     
% % %     zerobarInduced = zeros(1,size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1));
% % %     zerobarERP = zeros(1,size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1));
% % %     zerobarInduced(zerotime) = min(mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2)); % 0 ms from aud onset
% % %     zerobarERP(zerotime) = min(mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2))/2; % 0 ms from aud onset
% % %     
% % %     
% % %     
% % %     
% % %     figure(chanx+200)
% % %     bar([PLOT.PermResultsHighInduced,PLOT.PermResultsLowInduced],'stacked')
% % %     colormap([.8 .8 .8])
% % %     hold on
% % % 
% % %     shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),1),mean(PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)))))),2)',PLOT.(genvarname(char(strcat('C',num2str(plot_conds(1)),'_err')))),colorarray(1),.5);
% % %     hold on
% % %     
% % %     if NumConds>1
% % %         for j=2:NumConds
% % %             temp = PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))));
% % %             shadedErrorBar(1:size(temp,1),mean(temp,2)',PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j)),'_err')))),colorarray(j),.5); %{'color',[1,.647,0]},.5);
% % %         end
% % %     end
% % %     
% % %     plot(zerobarInduced)
% % %     hold off
% % %     title(strcat('Induced CH:   ', chan_names(chanx)),'FontWeight','bold');
% % %     set(gca, 'XTick', xtickrange)
% % %     set(gca, 'XTickLabel', xticklabel)
% % %     xlim([xaxis(1) xaxis(2)])
% % %     legend(legend_array);
% % %     
% % %     
% % %     
% % %     
% % %     
% % %     
% % %     
% % %     
% % %     figure(chanx+100)
% % %     bar([PLOT.PermResultsHighERP,PLOT.PermResultsLowERP],'stacked')
% % %     colormap([.8 .8 .8])
% % %     hold on
% % %    
% % %     shadedErrorBar(1:size(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),1),mean(PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)))))),2)',PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(1)),'_err')))),colorarray(1),.5);
% % %     
% % %     if NumConds>1
% % %         for j=2:NumConds
% % %             temp = PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j))))));
% % %             shadedErrorBar(1:size(temp,1),mean(temp,2)',PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)),'_err')))),colorarray(j),.5); %{'color',[1,.647,0]},.5);
% % %         end
% % %     end
% % %     
% % %     plot(zerobarERP)
% % %         hold off
% % %     title(strcat('Evoked CH:   ', chan_names(chanx)),'FontWeight','bold');
% % %     set(gca, 'XTick', xtickrange)
% % %     set(gca, 'XTickLabel', xticklabel)
% % %     xlim([xaxis(1) xaxis(2)])
% % %     if change_yaxis ==1
% % %         ylim([yaxis_ERP(1) yaxis_ERP(2)])
% % %     end
% % %     legend(legend_array);
% % %     
% % %     % pause
% % % end
% % % 
