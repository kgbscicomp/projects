m %% Analyse R-squared for all subjects across the conditions across all subjects 
% clear all
% clc

dataDir = 'C:\Users\gkarthik\Desktop\619\manus\all_conds\';

prestim = xlsread([dataDir 'prestim\' ['rsquared' 'prestim' '.xls']]);
faceonset = xlsread([dataDir 'faceonset\' ['rsquared' 'faceonset' '.xls']]);
facemove = xlsread([dataDir 'facemove\' ['rsquared' 'facemove' '.xls']]);
audonset = xlsread([dataDir 'audonset\' ['rsquared' 'audonset' '.xls']]);


allsubs_data = [prestim(end,:);faceonset(end,:);facemove(end-1,:);audonset(end,:)];


names = {'prestim'; 'faceonset'; 'facemove'; 'audonset'};
x = [];
congaud = allsubs_data(:,1);
err_congaud = allsubs_data(:,4);
congvis = allsubs_data(:,2);
err_congvis = allsubs_data(:,5);
audvis = allsubs_data(:,6);
err_audvis = allsubs_data(:,6);
A = shadedErrorBar(x,congaud,err_congaud,{'-or','markerfacecolor','red'},1);
hold on
B = shadedErrorBar(x,congvis,err_congvis,{'-og','markerfacecolor','green'},1);
hold on
C = shadedErrorBar(x,audvis,err_audvis,{'-ob','markerfacecolor','blue'},1);
legend([A.mainLine,B.mainLine,C.mainLine],'cong-aud','cong-vis','aud-vis','Location','northwest')
ylabel('r-squared cohen-d')
set(gca,'xtick',1:4,'xticklabel',names)
title('GLM between different test conditions');



