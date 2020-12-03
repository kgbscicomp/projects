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
plot_chansA=SPECS.plot_chansA;
plot_chansB=SPECS.plot_chansB;

NumConds = length(plot_conds);
zerotime = zerotime-xaxis(1);
zbaseline = zbaseline;

zbaseline = [round((zbaseline(1)*srate)+1):round((zbaseline(2)*srate)+1)];
zerotime = round((zerotime*srate)+1);
xaxis = [round((xaxis(1)*srate)+1) round((xaxis(2)*srate)+1)];
[~,a]=min(abs(times-(times(xaxis(1))+(SPECS.ISPCwindow/2))));
[~,b]=min(abs(times-(times(xaxis(2))-(SPECS.ISPCwindow/2))));
ISPCwindow = [a:b];
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
% for j=1:NumConds
%     PLOT.(genvarname(char(strcat('C',num2str(plot_conds(j))))))=[];
%     AVG.(genvarname(char(strcat('C',num2str(plot_conds(j)))))) = [];
% %     PLOT.(genvarname(char(strcat('avgC',num2str(plot_conds(j))))))=[];
% end
for chanxA = plot_chansA %[29 30 37 38 43:46 86:91]
    for chanxB = plot_chansB %[29 30 37 38 43:46 86:91]
     
    for j=1:NumConds
%         PLOT.(genvarname(char(strcat('ERP',num2str(plot_conds(j)))))) = [];
%         PLOT.(genvarname(char(strcat('avgERP',num2str(plot_conds(j)))))) = [];
        tempsizeERP=1;
        for i=1:numevents
            if Evoked_chan_bad_array(chanxA,i)==0 && ~ismember(i,bad_trials) && Evoked_chan_bad_array(chanxB,i)==0
                if stim_resamp(i,2)==plot_conds(j)
                    temp = (squeeze(evoked_sig(chanxA,:,i)));
                    PLOT.(genvarname(char(strcat('ERP_A',num2str(plot_conds(j))))))(:,tempsizeERP) = temp;
                    temp = (squeeze(evoked_sig(chanxB,:,i)));
                    PLOT.(genvarname(char(strcat('ERP_B',num2str(plot_conds(j))))))(:,tempsizeERP) = temp;
                    tempsizeERP=tempsizeERP+1;
                end
            end
            
        end
        
        for i=1:size(PLOT.(genvarname(char(strcat('ERP_A',num2str(plot_conds(j)))))),2)
            PLOT.(genvarname(char(strcat('ERP_A',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP_A',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP_A',num2str(plot_conds(j))))))(:,i),1));
            PLOT.(genvarname(char(strcat('ERP_B',num2str(plot_conds(j))))))(:,i) = PLOT.(genvarname(char(strcat('ERP_B',num2str(plot_conds(j))))))(:,i)-squeeze(mean(PLOT.(genvarname(char(strcat('ERP_B',num2str(plot_conds(j))))))(:,i),1));
        end        
        
    end
    
    % freq analysis
    for j=1:NumConds
        
        temp = PLOT.(genvarname(char(strcat('ERP_A',num2str(plot_conds(j))))));
        [mphaseA,mpowA] = bbpower_phase(temp, srate, Freq,FreqWidth,FreqCycles, 0, srate, variablefreq);
        temp = PLOT.(genvarname(char(strcat('ERP_B',num2str(plot_conds(j))))));
        [mphaseB,mpowB] = bbpower_phase(temp, srate, Freq,FreqWidth,FreqCycles, 0, srate, variablefreq);
        mphaseA = angle(mphaseA);
        mphaseB = angle(mphaseB);
        
        for i=1:size(mphaseA,2)
            
            % calculate Inter Signals Phase Coh - TRIALS
            tempA = (squeeze(mphaseA(:,i,:)));
            tempB = (squeeze(mphaseB(:,i,:)));
            temp = abs(mean(exp(1i*(tempA-tempB)'),2));
            meanx = mean(temp(zbaseline));
            PLOT.(genvarname(char(strcat('ISPCtrials',num2str(plot_conds(j))))))(i,:) = (temp-meanx);
            
            % calculate Inter Signals Phase Coh - TIME
            tempISPCtime = zeros(1,size(mphaseB,3));
            for k=ISPCwindow
                timex = [round(k-((SPECS.ISPCwindow/2)*srate)):round(k+((SPECS.ISPCwindow/2)*srate))];
                tempA = (squeeze(mphaseA(:,i,timex)));
                tempB = (squeeze(mphaseB(:,i,timex)));
                temp = abs(mean(exp(1i*(tempA-tempB)),2));
                tempISPCtime(k) = mean(temp);
            end
            meanx = mean(tempISPCtime(zbaseline));
            PLOT.(genvarname(char(strcat('ISPCtime',num2str(plot_conds(j))))))(i,:) = (tempISPCtime-meanx); 
            
            % calulate Spectral Coherence
            tempAphase = (squeeze(mphaseA(:,i,:)));
            tempBphase = (squeeze(mphaseB(:,i,:)));
            tempApow = (squeeze(mpowA(:,i,:)));
            tempBpow = (squeeze(mpowB(:,i,:)));
            
            %spec1 = mean(abs(sig1).^2,2);
            %spec2 = mean(abs(sig2).^2,2);
            %specX = abs(mean( abs(sig1).*abs(sig2) .* exp(1i*(angle(sig1)-angle(sig2))) ,2)).^2;
            
            % compute spectral coherence, using only requested time points
            spectcoher(fi,:) = specX(times2saveidx)./(spec1(times2saveidx).*spec2(times2saveidx));
            
            
            
            temp = abs(mean( abs(tempApow)'.*abs(tempBpow)' .* exp(1i*(tempAphase-tempBphase)'),2)).^2;
            meanx = mean(temp(zbaseline));
            PLOT.(genvarname(char(strcat('ISPCtrials',num2str(plot_conds(j))))))(i,:) = (temp-meanx);
            
            
            %specX = abs(mean( abs(sig1).*abs(sig2) .* exp(1i*(angle(sig1)-angle(sig2))) ,2)).^2;
    
            % compute spectral coherence, using only requested time points
            spectcoher(fi,:) = specX(times2saveidx)./(spec1(times2saveidx).*spec2(times2saveidx));
            
            
            
        end
        AVG.(genvarname(char(strcat('ISPCtrials',num2str(plot_conds(j))))))(chanxA,:,:) = PLOT.(genvarname(char(strcat('ISPCtrials',num2str(plot_conds(j))))));
    end
    end
end

% Process then plot

for chanxA = plot_chansA %[29 30 37 38 43:46 86:91]
    for chanxB = plot_chansB %[29 30 37 38 43:46 86:91]
    for j=1:NumConds

        figure(10+j)
        subplot(2,1,1)
        contourf(times,frex,PLOT.(genvarname(char(strcat('ISPCtrials',num2str(plot_conds(j)))))),40,'linecolor','none')
        set(gca,'xlim',plotaxis,'clim',[-.5 .5])
        title(strcat('CHANS:   ', chan_names(chanxA), ' to ', chan_names(chanxB),     ' Cond: ',legend_array(j),' ISPC TRIALS'),'FontWeight','bold');
        
        subplot(2,1,2)
        contourf(times,frex,PLOT.(genvarname(char(strcat('ISPCtime',num2str(plot_conds(j)))))),40,'linecolor','none')
        set(gca,'xlim',plotaxis,'clim',[-.5 .5])
        title(strcat('CHANS:   ', chan_names(chanxA), ' to ', chan_names(chanxB),     ' Cond: ',legend_array(j),' ISPC TIME'),'FontWeight','bold');
        
    end
    pause();
    end
end




