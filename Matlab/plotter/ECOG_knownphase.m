
chan_names = {};
for i=1:size(EEG.chanlocs,2)
chan_names{i} = EEG.chanlocs(1, i).labels;
end

DATA = double(mean(EEG.data,3)');
EEGtimes = EEG.times;
zerotime = 1025;
Freq = 10; % in hz
angle_result = [];
for el=1:size(DATA,2)
    DATAellip =butterpass_eeglabdata(DATA(:,el),[9.99 10.01],EEG.srate,3,60);
    [params,yest] = sinefit(DATAellip,EEGtimes/1000,10,0,0,0);
    % plot([yest,DATAellip])
    test = angle(hilbert(yest));
    angle_result(el,1) = test(zerotime);
end
angle_result_deg = rad2deg(angle_result);



DATA = double(mean(EEG.data,3)');
EEGtimes = EEG.times;
zerotime = 1025;
Freq = 10; % in hz
angle_result = [];
for el=1:size(DATA,2)
    [params,yest] = sinefit(DATA(:,el)',EEGtimes/1000,10,0,0,0);
    test = angle(hilbert(yest));
    angle_result(el,1) = test(zerotime);
end
angle_result_deg = rad2deg(angle_result);
% this will put the values from -180 to 180
% use these max/mins at the extremes in freesurfer colors

save('10hzAngleData.mat','DATA','EEGtimes');

% show the plot for a representative channel as a figure


% load data from RAS script
setenv SUBJECTS_DIR /Applications/freesurfer/subjects/
cd /Applications/freesurfer/subjects/KellyA
cd surf
load 10HZ_KellyA.mat
sol_lh  = angle_result_deg;
vertices2 = [];
sol2 = [];

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
sol2_lh = sol2;
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

threshold = 0;
vertexsize = size(vertex_coords_pial_lh,1)-1;
thresh_vertex = [1:vertexsize]';
thresh_sol = zeros(vertexsize,1)+500;
counter = 1;
for i=[vertices2_lh(:,1)]'
    thresh_sol(i,1) = sol2_lh(counter,1);
    counter=counter+1;
end



data2mgh(thresh_vertex, thresh_sol, '10HZ_KellyA_phase_long.mgh','KellyA','lh',[],[],[],[]);








%
%
%
%
% Pick some range of phases to be shown, thresh, and show a movie of the
% changing phase.

% DIDN"T WORK WELL BECAUSE OF FREEVIEW


chan_names = {};
for i=1:size(EEG.chanlocs,2)
chan_names{i} = EEG.chanlocs(1, i).labels;
end

DATA = double(mean(EEG.data,3)');
EEGtimes = EEG.times;
timerange = [5121:15361];
Freq = 10; % in hz
angle_result = [];
rms_result = [];
maxmin_result = [];
for el=1:size(DATA,2)
    DATAellip =butterpass_eeglabdata(DATA(:,el),[Freq-.1 Freq+.1],EEG.srate,3,150);
    [params,yest] = sinefit(DATAellip(timerange,:),EEGtimes(timerange)/1000,10,0,0,0);
    % plot([yest,DATAellip])
    test = angle(hilbert(yest));
    angle_result(el,:) = test';
    rms_result(el,1) = rms(DATAellip(timerange,:));
    maxmin_result(el, 1) = max(yest)-min(yest);
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

AudCtxEl = 41; %LTMG23 for KellyA 40 or 41
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


















% angle_result_deg = rad2deg(angle_result);

% load data from RAS script
setenv SUBJECTS_DIR /Applications/freesurfer/subjects/
cd /Applications/freesurfer/subjects/KellyA
cd surf
load 10HZ_KellyA.mat
sol_lh  = angle_result;
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
sol2_lh = sol2;
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


threshold = [87.5 92.5]; % number of degs symmetric around 0
randX = randi(10,size(sol_lh,2),1);
thresh_sol = zeros([size(sol2_lh)])+.1;
for k=1:size(sol_lh,2) % k is for time
for i=1:size(vertices2_lh,1)
    if sol2_lh(i,k) > threshold(1) && sol2_lh(i,k) < threshold(2)
        thresh_sol(i,k) = sol2_lh(i,k);
    else
        thresh_sol(i,k) = randX(k,1);
    end
end
end


data2mgh(vertices2_lh, thresh_sol, '10HZ_KellyA_phase_thresh.mgh','KellyA','lh',[],[],[],[]);














