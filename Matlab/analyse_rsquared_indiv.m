%%Plot errorbars for individual subjects for all four conditions (prestim,
%%faceonset, facemove, audonset).

dataDir = 'C:\Users\gkarthik\Desktop\619\manus\all_conds\';

audonset = xlsread([dataDir 'audonset\' ['rsquared' 'audonset' '.xls']]);
audonset_data = audonset(1:end-1,1);
audonset_err = audonset(1:end-1,4);

faceonset = xlsread([dataDir 'faceonset\' ['rsquared' 'faceonset' '.xls']]);
faceonset_data = faceonset(1:end-1,1);
faceonset_err = faceonset(1:end-1,4);

facemove = xlsread([dataDir 'facemove\' ['rsquared' 'facemove' '.xls']]);
facemove_data = facemove(1:end-1,1);
facemove_err = facemove(1:end-1,4);

prestim = xlsread([dataDir 'prestim\' ['rsquared' 'prestim' '.xls']]);
prestim_data = prestim(1:end-1,1);
prestim_err = prestim(1:end-1,4);

errorbar(audonset_data,audonset_err,'r');hold on;
errorbar(faceonset_data,faceonset_err,'g');hold on;
errorbar(facemove_data,facemove_err,'b');hold on;
errorbar(prestim_data,prestim_err,'k')
legend('audonset','faceonset','facemove','prestim')