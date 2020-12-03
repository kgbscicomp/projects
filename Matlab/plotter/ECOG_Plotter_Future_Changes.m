%% Needs fixing

% In ECOG_filter_analysis_v7 line 172
% Length of bad_chans1 instead? Maybe not, but not plotting all listed
% chans. Only a prob when chans named a number.

% ITPC stats doesn't incorp new changes yet (from other stats)


%% Future Changes

% Cluster-wise pvals

% Ability to select a portion of the data. Useful for hold-out analyses.
% E.g., analyze 75% of data in individual subject design to identify
% electrodes where a difference is present, then pass the remaining
% (untested) 25% to group analysis. Easy to just modify bad trials array in
% the short term. Better to use stratified selections so equal balance of
% conditions across training and test sets.

% Add causal frequency filters

% Simplify non-parametric analysis

% Analyses for phase-based-measures

% Need to test with artifical data (raw ECoG + signal) to see when each
% specification optimal
% Maybe even include options for artifical data testing

% Create option to run stats as a localizer
% This way it samples even numbers of trials before doing one-samp ttest
% Use Cluster for localizing times and freqs?

% Possible to average data across channels and condition then run a search
% to optimize the parameters using the localizer? E.g. Baseline, Freq,
% Cycles

% add local min-max detector (peak-to-peak in 200 ms windows to chan detector)

% make plotter a function

% when doing laplacian, double check that it's pulling the correct chan
% numbers in instances channels are out of order (e.g., 1128UC)

% even if not doing laplacian, still create the matlab array w matching
% vertices to make it easier to plot fMRI-like data

% for noisy ERP baselines with time-locked intrisic componenents (e.g.,
% alpha), use a baseilne that is a multiple of the freq of the dominant
% component.

% add in option to lowpass at .1, prob instead of detrend
% add both low pass and high pass specs

% change sig boxes to go from bottom to half-way to data (or max half to
% xaxis 0 line)

% Since we're doing parametric stats now, might make sense to exclue
% outliers from ERP and ERSP somehow.

% enter key bypassing trial rejection instead of 0

% when plotting bad channels, maybe just overlap all epochs from same chan
% Can also do one figure per channel

% what about taking the absolute min value across all trials If powerbaseline is empty

% % ******
% Global field power

% Total Spectral Power
% Separately, can just take raw power of ERP as an option.
% Would be the exact same thing as the filtered EEG response, but not filt
% Mean baseline raw power. no longer care about polarity. can combine data
% across electrodes then. still keeps ms precision but only in terms of
% power. Downside is we loose phase informationv (and so some statistical power if ITPC contributing to ERP).
% Could show ITPC separately; but would we need to filter for ITPC?


%% Create group analysis script

% give it electrodes from each subject
% As long as the plotter is a function, will be easy for it to batch
% use multi-level modeling

% 1D array of subject matfile/locs
% 2D array of channels to include (same column = mult elec/sub, ignore 0s)



%% RAS Plotter

% make script to pull channel labels and electrodes (plus grid array) and
% can click on electrodes to renumber them if never numbered and then
% generate laplacian from that point. would remove need for the box-clicker
% program. If numbers already in csv, would just make it's own grid from
% the plane of the RAS coordinates.

% Take midway point between contacts
% Maybe find the closest pial vertex if non-depth

% Display brain w elec loc on each plot


%% New stats
% conjuction analysis: require both sig
% then use orig method on non-sig ones
% 
% * do simultaneously, and choose if either true
% 
% or 
% 
% 
% use pooled standard deviatio
% 
% 
% 
% treat it as diff of conds relative to order of trials
% 
% 
% 
% 
% both pos, both neg, or one above or below
% 
% 
% ???
% Also can test non sig electrodes (non ROI selected electrodes) with bonf corrections, maybe elec grouping.