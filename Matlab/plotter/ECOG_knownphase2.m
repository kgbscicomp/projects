% pull data from eeglab, downsample, filter, push back
EEG = pop_resample( EEG, 1000);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'setname','resampled','gui','off'); 
DATA = double(EEG.data');
Freq = 10; % in hz
Rs = 1000;
DATAfilt = [];
for el=1:size(DATA,2)
    DATAfilt(el,:) =elliptic_eeglabdata(DATA(:,el),[Freq-.1 Freq+.1],EEG.srate,3,120);
end
EEG.data = DATAfilt;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

% extract epochs

DATA = double(mean(EEG.data,3)');
R2_input = [];
for i=1:size(EEG.data,1)
    for j=1:size(EEG.data,3)
        for k=1:size(EEG.data,2)
            R2_input(j,1,k,i) = (EEG.data(i,k,j));
        end
    end
end


EEGtimes = EEG.times;
timerange = [1:500]; %[5121:15361];
angle_result = [];
rms_result = [];
maxmin_result = [];
ITC_data = [];
z = [];
for el=1:size(DATA,2)
    [params,yest,yred] = sinefit(DATA(timerange,el),EEGtimes(timerange)/1000,10,0,0,0);
%     plot([DATAfilt(el,timerange)',yest])
    angle_result(el,:) = angle(hilbert(yest))';
    rms_result(el,1) = rms(DATA(timerange,el))';
    maxmin_result(el, 1) = max(yest)-min(yest);
    z_raw = R2_input(:,:,:,el);
    for i=1:size(z_raw,1)
        z(i,1,:) = hilbert(squeeze(z_raw(i,1,:)));
    end
    size_z = size(z);
    z = z - repmat(mean(z), [size_z(1) 1 1]);
    ang = angle(z);
    R2 = squeeze(abs(mean(exp(1i*(2*ang)))));
    ITC_data(el,1) = mean(R2);
end



% load data from the distance array

[vertex_coords_sphere_lh] = read_surf('lh.sphere');
spherevertices = [];
for i=1:size(tempvertices,1) 
    tempvertex = vertex_coords_sphere_lh(tempvertices(i,1)+1,:);
    [azimuth,elevation,r] = cart2sph(tempvertex(1),tempvertex(2),tempvertex(3));
     spherevertices(i,:) = [azimuth,elevation,r];
end

% build an 1xN array of distances from aud ctx

AudCtxEl = 39; %LTMG23 for KellyA 40 or 41
% AudCtxEl = 42; %LTMG23 for RL 34 or 42
% AudCtxEl = 39; %LTMG23 for ER 16 or 39
% AudCtxEl = 50; %LTMG23 for Tamsen

sphereDistAud = zeros(size(tempvertices,1),1)-1;
for i=1:size(sphereDistAud,1)
    [arclen,az] = distance([spherevertices(AudCtxEl,1:2)],[spherevertices(i,1:2)],'radians');
    sphereDistAud(i,1) = arclen;
end

[rho pval] = circ_corrcl(angle_result(:,1), sphereDistAud)
[rho pval] = corr(rms_result, sphereDistAud)
[rho pval] = corr(maxmin_result, sphereDistAud)
[rho pval] = corr(ITC_data, sphereDistAud)



% Save data in mgh




% extract -2 2 epochs


DATA = double(mean(EEG.data,3)');
EEGtimes = EEG.times;
angle_result = [];
angle_result_mean = [];
angle_result_stdev = [];
yest_array = [];
FreqX = 10;
TestTimes = min(EEGtimes):1000/FreqX:max(EEGtimes);
EEGTestTimes = [];
for i=1:length(TestTimes)
    EEGTestTimes(i) = find(round(EEGtimes)==TestTimes(i));
end
EEGTestTimes = [EEGTestTimes length(EEGtimes)];

for el=1:size(DATA,2)
    [params,yest] = sinefit(DATA(:,el),EEGtimes(:)/1000,FreqX,0,0,0);
    angle_result(el,:) = angle(hilbert(yest))';
    yest_array(el,:) = yest;
    yest_array_dist = [];
% % % % %     for i=1:length(EEGTestTimes)-1
% % % % %         [params,yest] = sinefit(DATA(EEGTestTimes(i):EEGTestTimes(i+1),el),EEGtimes(EEGTestTimes(i):EEGTestTimes(i+1))/1000,FreqX,0,0,0);
% % % % %         % plot([yest';DATA(EEGTestTimes(i):EEGTestTimes(i+1),el)']')
% % % % %         temp = angle(hilbert(yest))';
% % % % %         yest_array_dist(i) = temp(1);
% % % % %     end
% % % % %     angle_result_mean(el,1) = circ_mean(yest_array_dist');
% % % % %     angle_result_stdev(el,1) = circ_std(yest_array_dist');
end
angle_result2 = angle_result(:,find(EEGtimes==0));

chan_names = {};
for i=1:size(EEG.chanlocs,2)
chan_names{i} = EEG.chanlocs(1, i).labels;
end

datatimes = [find(round(EEGtimes)==-100):find(round(EEGtimes)==100)];
el=40;
plot([yest_array(el,datatimes);DATA(datatimes,el)']')




for el=1:size(DATA,2)
    
    h = rose(angle_result2(el,1),360);
    set(h,'linewidth',10)
    mycmd = ['print -painters -dpdf ',char(strcat(chan_names(el),'.pdf'))];
    eval(mycmd)
end










% load data from RAS script
setenv SUBJECTS_DIR /Applications/freesurfer/subjects/
cd /Applications/freesurfer/subjects/KellyA
cd surf
load 10HZ_KellyA.mat

% sol_lh  = angle_result(:,1);
sol_lh  = angle_result2;
vertices2 = [];
sol2 = [];
tpoints = size(sol_lh,2); % number of timepoints

for k=1:size(sol_lh,2) % k is for time
    solcount = 0;
    for i=1:size(leftelectrodes,1);
        for j=1:size(vertex_coords_pial,1)
            if tempvertexdist(j,i)<=electroderadius
                solcount = solcount+1;
                vertices2(solcount,1) = j-1;
                sol2(solcount,k) = sol_lh(i,k);
            end
        end
    end
    k
end
doublevertex = [];
counter = 1;
for i=1:size(vertices2,1)

    test1 =find(vertices2(i,1)==vertices2);
    if size(test1,1)>1
        doublevertex(counter,1:length(test1')) = test1';
        counter = counter+1;
    end
end
doublevertex2 = unique(doublevertex,'rows');
sol2_lh = rad2deg(sol2);
vertices2_lh = vertices2;

% average shared vertices
for k=1:size(sol_lh,2) % k is for time
for i=1:size(doublevertex2,1)
    tempvtx = nonzeros(doublevertex2(i,:))';
    for j=1:length(tempvtx)
        sol2_lh(doublevertex2(i,j),k) = rad2deg(circ_mean(sol2(tempvtx,k)));
    end
end
end

vertexsize = size(vertex_coords_pial_lh,1)-1;
thresh_vertex = [1:vertexsize]';
thresh_sol = zeros(vertexsize,size(sol2_lh,2))+500;
counter = 1;
for i=[vertices2_lh(:,1)]'
    thresh_sol(i,:) = sol2_lh(counter,:);
    counter=counter+1;
end


data2mgh(thresh_vertex, thresh_sol, 'Cond1_RL_Angle-Resid.mgh','KellyA','lh',[],[],[],[]);



