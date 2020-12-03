%% Fits a GLM and calculates R-squared and SE and saves them out and displays a plot. 
%% Set show_plots = 0 in line 30 if you want to turn off the plots (total of 4 figures). 
% clear
% clc
%% Folder with subjects' electrodes and p/d values- the parent folder with all time-points
dataDir = 'C:\Users\gkarthik\Desktop\619\manus\all_conds\';

% Subject IDs
subjids = {'1187HF',...
            '1167HF',...
            '1164UM'...
            '1162HF'...
            '1151HF'...
            '1132UC'...
            '1125UM'...
            '1124UC'...
            '1119UC'...
            '1118UC'...
            '1111UC'...
            '1110UC'...
            '1108UC'};
analysis_window = 'audonset\';
% analysis_window = 'prestim\';
% analysis_window = 'faceonset\';
% analysis_window = 'facemove\';

show_plots = 1; %% Set show_plots=0 if you don't want GLM plots. 
elec_dir = dataDir;

%% Identify electrodes and coordinates to be plotted
% Should be able to do computationally. Here I've just pulled the elec
% names manually for 3 subjects
if(strcmp(analysis_window,'audonset\'))
    an_window = 'audonset\';
elseif(strcmp(analysis_window,'prestim\'))
    an_window = 'prestim\';
elseif(strcmp(analysis_window,'faceonset\'))
    an_window = 'faceonset\';
elseif(strcmp(analysis_window,'facemove\'))
    an_window = 'facemove\';
end

for i = 1:length(subjids)
    subj = num2str(str2double(regexp(subjids{i},'[\d.]+','match')));
    ROI1.(['elec' subj]) = load([elec_dir,subjids{i},'\',an_window,'Elec_STG.mat']);
    d_and_p = ROI1.(['elec' subj]).electrodes_in_STG;
    remove_rows = any(cellfun(@isempty, d_and_p), 2); 
    d_and_p(remove_rows,:) = [];  
    ROI1.(['elec' subj]) = d_and_p;
end

subjs = fieldnames(ROI1);
all_subjs_data = {};
all_subjs_data = arrayfun(@(n) [all_subjs_data;ROI1.(subjs{n})], 1:length(subjs),...
                'UniformOutput', false);
all_subjs_data = vertcat(all_subjs_data{:});
all_subjs_data(:,1) = strrep(all_subjs_data(:,1),'''','');

Rsquared_matrix = {};
Rsquared_matrix(:,1) = subjids;

counter = 1;
for i = 1:length(subjids)+1
    if(i == length(subjids)+1)
        corr_data = all_subjs_data;
        Rsquared_matrix{i,1} = 'All_subs';
        subj = 'AllSubs';
    else
        subj = num2str(str2double(regexp(subjids{i},'[\d.]+','match')));
        corr_data = ROI1.(['elec' subj]);
    end
    congruent = cell2mat(corr_data(:,7));
    audio = cell2mat(corr_data(:,9));
    visual = cell2mat(corr_data(:,11));
    if(rem(i+3,4) == 0)
        figure(counter)
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        counter = counter + 1;
        col = 1;
    end
    %%Congruent and audio
    congaud = fitlm(congruent,audio);
    Rsq2_congaud = congaud.Rsquared.Adjusted;
    err_congaud = congaud.Coefficients.SE(2);
    subplot(4,3,col);plot(congaud);
    col = col+1;
    title(['GLM b/w cong and aud cohens d' ': ' subj])
    xlabel('Congruent cohens-d')
    ylabel('Audio cohens-d')
    legend(['rsq:' num2str(Rsq2_congaud)])
    %%Congruent and visual
    congvis = fitlm(congruent,visual);
    Rsq2_congvis = congvis.Rsquared.Adjusted;
    err_congvis = congvis.Coefficients.SE(2);
    subplot(4,3,col);plot(congvis);
    col = col+1;
    title(['GLM b/w cong and vis cohens d' ': ' subj])
    xlabel('Congruent cohens-d')
    ylabel('Visual cohens-d')
    legend(['rsq:' num2str(Rsq2_congvis)])
    %%Audio and visual
    audvis = fitlm(audio,visual);
    Rsq2_audvis = audvis.Rsquared.Adjusted;
    err_audvis = audvis.Coefficients.SE(2);
    subplot(4,3,col);plot(audvis);
    title(['GLM b/w aud and vis cohens d' ': ' subj])
    xlabel('Audio cohens-d')
    ylabel('Visual cohens-d')
    legend(['rsq:' num2str(Rsq2_audvis)])
    col = col + 1;
    if(col > 12)
        col = 1;
        print('-dpng','-r300',[dataDir analysis_window 'rsquared' analysis_window(1:end-1) num2str(counter-1) '.png'])
    end
    if(i == length(subjids)+1)
        print('-dpng','-r300',[dataDir analysis_window 'rsquared' analysis_window(1:end-1) num2str(counter-1) '.png'])
    end
    
    Rsquared_matrix{i,2} = Rsq2_congaud;
    Rsquared_matrix{i,3} = Rsq2_congvis;
    Rsquared_matrix{i,4} = Rsq2_audvis;
    Rsquared_matrix{i,5} = err_congaud;
    Rsquared_matrix{i,6} = err_congvis;
    Rsquared_matrix{i,7} = err_audvis;
end

if(show_plots==0)
    close all;
end

header = {'SubID','r2-Cong-Aud','r2-Cong-Vis', 'r2-Aud-Vis','SE-congaud','SE-congvis','SE-audvis'};
Rsquared_matrix = [header;Rsquared_matrix];

xlswrite([dataDir analysis_window ['rsquared' analysis_window(1:end-1) '.xls']],Rsquared_matrix)
