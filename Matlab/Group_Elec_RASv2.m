%% Choose analysis window in lines 27-30. New directories will be created 
%% under data_dir for writing out image files

clear
clc
% Dropbox dir
DatDir = 'C:\Users\gkarthik\Dropbox\Electrode_Registration\';

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
% Plot left hemi
Hemi = 'LH';

% analysis_window = 'audonset\';
% analysis_window = 'prestim\';
% analysis_window = 'faceonset\';
analysis_window = 'facemove\';
% plot_cond = 'Congruent';
% plot_cond = 'Audio';
plot_cond = 'Visual';
% Will create an 'image' folder in the following directory
mkdir([dataDir analysis_window])
cd([dataDir analysis_window])

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

% Elec names
% ROI.elec = {'LT3','LT4','LT5','LT6','LT7','LT10','LT11','LT13','LT14','LT15','LT17','LT22','LT24','LT30','LT31','LT32','LT33','LT34','LT35','LT36','LT37','LT38','LT39','LT40','LT41','LT42','LT43','LT44','LT45','LT46','LT47','LT48','LT49','LT50','LT51','LT52','PT6','PT7','PT8'};
% ROI.elec = [ROI.elec {'LAM5','LAM6','LAM7','LAM8','LAM9','LMH6','LMH7','LMH8','LMH9','RPB2','RPB3','RPB4','RPB5','RPB8','RPB9','RPB10','RPB11','RPB13','RPB14','RPB15','RPB16','STG4','STG5'}];
% ROI.elec = [ROI.elec {'LT1_b','LT2_b','LT3_b','LT4_b','LT5_b','LT7_b','LT8_b','LT9_b','LT10_b','LT11_b','LT12_b','LT13_b','LT14_b','LT15_b','LT16_b','LT17_b','LT18_b','LT19_b','LT20_b','LT21_b','LT22_b','LT23_b','LT24_b','LT25_b','LT26_b','LT27_b','LT28_b','LT29_b','LT30_b','LT31_b','LT32_b','LT33_b','LT34_b','LT35_b','LT36_b','LT37_b','LT38_b','LT39_b','LT40_b','LT41_b','LT42_b','LT43_b','LT44_b','LT45_b','LT46_b','LT47_b','LT48_b','LT49_b','LT50_b','LT51_b','LT52_b','LT53_b','LT54_b','LT55_b','LT56_b','LT57_b','LT58_b','PT9_b','PT12_b','PT13_b'}];

subjs = fieldnames(ROI1);
elecs_to_plot = {};
elecs_to_plot = arrayfun(@(n) [elecs_to_plot;ROI1.(subjs{n})], 1:length(subjs),...
                'UniformOutput', false);
elecs_to_plot = vertcat(elecs_to_plot{:});
elecs_to_plot(:,1) = strrep(elecs_to_plot(:,1),'''','');

ROI.elec = elecs_to_plot(:,1);
ROI.Hemi = Hemi;

subj = num2str(str2double(regexp(subjids{1},'[\d.]+','match')));
ROI.SUB = repmat({subjids{1}},1,size(ROI1.(['elec' subj]),1));
for i = 2:length(subjids)
    subj = num2str(str2double(regexp(subjids{i},'[\d.]+','match')));
    ROI.SUB = [ROI.SUB repmat({subjids{i}},1,size(ROI1.(['elec' subj]),1))];
end

% subjid_parts = regexp(subjids{i}, '(?<subj>\d+)(?<hosp>[A-Za-z]+)', 'names')
%% Identify the electrode coodinates

% separate electrode group and numbers
ROI.Group = regexp(ROI.elec,'\D*','match');
ROI.Contact = regexp(ROI.elec,'\d*','match');
ROI.Contact = str2double(cellfun(@(x) cell2mat(x),ROI.Contact,'un',0));

% Iterate through subjects
ROI.SubNames = unique(ROI.SUB);

electrode_RAS_ALL = [];
counter = 1;
for i=1:length(ROI.SubNames)
    
   % load in the Electrode_Labels for sub
   load([DatDir,ROI.SubNames{i},'\Electrodes\Electrode_Labels.mat'])
   
   % Find current electrodes
   curr_el = find(strcmp(ROI.SubNames{i},ROI.SUB));
   
   % Iterate through the electrodes for the subject
   % Pull the RAS coords for all electrodes
   electrode_RAS = [];
   for j=1:length(curr_el)
       curr_el_name = ROI.Group{curr_el(j)}{:};
       curr_el_number = ROI.Contact(curr_el(j));
       
       % See if name case needs to be changed
       if isfield(Electrode_Labels.(genvarname([ROI.Hemi,'_MNI'])).Bipolar,curr_el_name)==0
           curr_el_name = lower(curr_el_name);
           if isfield(Electrode_Labels.(genvarname([ROI.Hemi,'_MNI'])).Bipolar,curr_el_name)==0
               curr_el_name = upper(curr_el_name);
           end
       end
       
       % Pull the RAS coords
       index = find(Electrode_Labels.(genvarname([ROI.Hemi,'_MNI'])).Bipolar.(genvarname(curr_el_name))(:,1)==curr_el_number);
       temp = Electrode_Labels.(genvarname([ROI.Hemi,'_MNI'])).Bipolar.(genvarname(curr_el_name))(index,5:7);
       electrode_RAS(j,:) = temp;
       electrode_RAS_ALL(counter,:) = temp;
       counter = counter+1;
   end
end

%% Plot data
% This section will prep the data to be passed along to the electrode
% plotter script.

run_analysis_num = 1;

% Compare anatomically auditory electrodes that are either sig and non-sig
% in response to congruent AV speech
% Elec_Data.Cong_Sig_Values_Uncorr = [3.77000000000000e-11;1.17000000000000e-22;0.00152267400000000;4.46000000000000e-26;1.51000000000000e-39;0.164382267000000;0.923350403000000;3.51000000000000e-10;8.35000000000000e-11;4.90000000000000e-30;0.000307092000000000;3.19000000000000e-06;0.00235397300000000;1.68000000000000e-11;0.423291754000000;0.000295771000000000;0.0102802500000000;0.0439336490000000;0.542602885000000;0.000195136000000000;0.00104172600000000;0.108299427000000;0.0121675990000000;0.913208390000000;0.0426124670000000;0.399974854000000;0.000990583000000000;6.16000000000000e-05;0.0118574970000000;0.0745146360000000;0.00139212000000000;0.0361010490000000;0.114327048000000;0.0140081530000000;0.00211849600000000;0.237301137000000;0.0224718000000000;8.28000000000000e-05;0.217183121000000;7.94000000000000e-12;3.90000000000000e-06;4.23000000000000e-07;0.000120612000000000;1.48000000000000e-06;0.0363524470000000;0.0146609860000000;0.539583165000000;1.94000000000000e-05;1.25000000000000e-27;2.18000000000000e-09;1.13000000000000e-16;4.29000000000000e-12;2.75000000000000e-32;3.31000000000000e-07;0.263521733000000;5.47000000000000e-15;1.87000000000000e-13;2.07000000000000e-28;0.0558934630000000;1.34000000000000e-05;3.15000000000000e-25;7.08000000000000e-65;0.239773526000000;0.0466217450000000;0.954514846000000;0.299943719000000;0.297785522000000;0.994337552000000;0.585158120000000;0.682981893000000;0.0120086100000000;0.969409583000000;0.708352254000000;0.979719385000000;0.238932421000000;0.775979498000000;0.465901997000000;0.999892229000000;0.260387982000000;0.437241624000000;0.666408187000000;0.991655366000000;0.357273546000000;0.0980983360000000;0.00166276800000000;0.0340144700000000;0.000974924000000000;0.997529341000000;0.999976691000000;0.00768311600000000;1.33000000000000e-08;0.130113886000000;1.41000000000000e-05;0.313814115000000;0.766712700000000;0.583858349000000;4.75000000000000e-05;2.31000000000000e-11;0.0412872500000000;0.00263066600000000;0.810221832000000;0.0178435170000000;0.00283233200000000;1.30000000000000e-05;0.00413884000000000;0.0545339440000000;0.297181577000000;0.957840170000000;0.0306685310000000;0.702634263000000;0.911433986000000;0.00282678400000000;0.727395595000000;0.788777880000000;0.959569822000000;0.165886748000000;0.963296171000000;0.965435487000000;0.181947196000000;0.618134744000000;0.568145035000000;9.34000000000000e-08];
% Elec_Data.Cong_Sig_Color = repmat([204 32 60]/255,size(Elec_Data.Cong_Sig_Values_Uncorr,1),1);
% Elec_Data.Cong_NonSig_Color = repmat([0 0 0]/255,size(Elec_Data.Cong_Sig_Values_Uncorr,1),1);

if(strcmp(plot_cond,'Congruent'))
    plot_col = 6;
elseif(strcmp(plot_cond,'Audio'))
    plot_col = 8;
elseif(strcmp(plot_cond,'Visual'))
    plot_col = 10;
end
Elec_Data.Cong_Sig_Values_Uncorr = cell2mat(elecs_to_plot(:,plot_col));
Elec_Data.Cong_Sig_Color = repmat([204 32 60]/255,size(Elec_Data.Cong_Sig_Values_Uncorr,1),1);
Elec_Data.Cong_NonSig_Color = repmat([0 0 0]/255,size(Elec_Data.Cong_Sig_Values_Uncorr,1),1);

% Add other data below...

% 1. Exp1: Electrode Locs Fig - All red or black
if run_analysis_num==1
    analysis_array = 'ROIs';
    Auditory_El = electrode_RAS_ALL;
    
    % correct p-vals for mult comparisons
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(Elec_Data.Cong_Sig_Values_Uncorr);
    
    % Find sig and non-sig
    Auditory_Color = zeros(size(Elec_Data.Cong_Sig_Values_Uncorr,1),3);
    Auditory_Color(find(adj_p<.05),:) = Elec_Data.Cong_Sig_Color(find(adj_p<.05),:);
    Auditory_Color(find(adj_p>=.05),:) = Elec_Data.Cong_NonSig_Color(find(adj_p>=.05),:);
end

% 2. Plot P-vals from cong condition only
if run_analysis_num==2
    analysis_array = 'Sig_Cong';
    Auditory_El = electrode_RAS_ALL;
    
    % correct p-vals for mult comparisons
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(Elec_Data.Cong_Sig_Values_Uncorr);
    
    % Color elecs by sig
    Auditory_Color = adj_p;

end


% radius of elecs
sphere_rad = 2; % mm

filenames = analysis_array;
alphaLevel = .3; % 2nd run

lambda = 1; %% Lower lambda values = more transparency
[nearest_elec,idx] = max(electrode_RAS_ALL(:,3));
all_elecs_bkup = electrode_RAS_ALL;
all_elecs_bkup(:,4) = arrayfun(@(n) sqrt(sum(all_elecs_bkup(idx,:)-all_elecs_bkup(n,:)).^2), 1:length(all_elecs_bkup));
all_elecs_bkup(:,4) = all_elecs_bkup(:,4)/max(all_elecs_bkup(:,4));
all_elecs_bkup(:,4) = arrayfun(@(n) 1/(1+exp(-lambda*all_elecs_bkup(n,4))), 1:length(all_elecs_bkup));

% Update colors if pval
for alphaLevel_temp=[1 alphaLevel]
    if size(Auditory_Color,2)==1
        
        n_colors=10000000;
        c = autumn(round(n_colors));
        c = flipud(c);
        
        % Add Zero
        c = [c;zeros(1,3)];
        
        sigThreshold = 0.05;
        
        %map p-values to colors
        stats=Auditory_Color;
        temp = find(stats>sigThreshold);
        stats(temp) = zeros(length(temp),1)+sigThreshold;
        scaledStats = rescale([stats; sigThreshold; 0],1, n_colors+1);                  %scale values; concat min and max p to enforce correct color range
        scaledStats = scaledStats(1:end-2);                                         %remove min and max p
        Auditory_Color = c(round(scaledStats),:);
        
    end
    
% Electrode_Plotter_Group_McGurk_v2([filenames,'-alpha-',num2str(alphaLevel_temp),'-Img_'],Auditory_Color,Auditory_El,alphaLevel_temp,sphere_rad,DatDir)
Electrode_Plotter_Group_McGurk_v3([filenames,'-alpha-',plot_cond,num2str(alphaLevel_temp),'-Img_'],Auditory_Color,all_elecs_bkup,alphaLevel_temp,sphere_rad,DatDir)

end



%% plot color bar
x = [0.5 1 ] ;
y = linspace(0,10) ;
[X,Y] = meshgrid(x,y) ;
Z = X+Y ;
figure
subplot(1,10,1)
surf(X,Y,Z,'edgecolor','none')
view(2)

% add zeros
c_new = [c;zeros(round(n_colors*.18),3)];

colormap(flipud(c_new));
% shading interp ;
hAxes = gca;
% hAxes.XRuler.Axle.LineStyle = 'none';  
axis off
set(gca, 'YTick', []);
set(gca,'XTickLabel',{' '})
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';   % remove x-axis



% temp = linspace(1.7,10,54);
% for i=1:4:54
%     text(1.05,temp(i),[num2str(i+1),' Hz'])
% end
% 
fontsizeZ = 36;
text(1.05,.4,'p = N.S.','FontSize',fontsizeZ,'fontname','Helvetica Neue','FontWeight','bold')
text(1.05,1.7,'p = .05','FontSize',fontsizeZ,'fontname','Helvetica Neue','FontWeight','bold')
text(1.05,10,'p = .0001','FontSize',fontsizeZ,'fontname','Helvetica Neue','FontWeight','bold')


