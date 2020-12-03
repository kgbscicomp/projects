%conditions = {}; %% copy paste phoneme dictionary from ecog conditions excel sheet
temp = 1:size(OUTPUT.stim,1);temp(find(OUTPUT.stim(:,2)==5))=[];

IndexC = strfind(conditions(:,1), 'BA');
IndexC = find(not(cellfun('isempty', IndexC)));
OUTPUT.stim(temp(IndexC),2) = OUTPUT.stim(temp(IndexC),2)+100;
C = cell(1, length(IndexC));C(:) = {'999'};conditions(IndexC,1) = C;conditions(IndexC,2) = C;
 

IndexC = strfind(conditions(:,1), 'DA');
IndexC = find(not(cellfun('isempty', IndexC)));
OUTPUT.stim(temp(IndexC),2) = OUTPUT.stim(temp(IndexC),2)+200;
conditions(IndexC,1) = num2cell(repmat(999,1,length(IndexC)));
C = cell(1, length(IndexC));C(:) = {'999'};conditions(IndexC,1) = C;conditions(IndexC,2) = C;
 

IndexC = strfind(conditions(:,1), 'TA');
IndexC = find(not(cellfun('isempty', IndexC)));
OUTPUT.stim(temp(IndexC),2) = OUTPUT.stim(temp(IndexC),2)+300;
conditions(IndexC,1) = num2cell(repmat(999,1,length(IndexC)));
C = cell(1, length(IndexC));C(:) = {'999'};conditions(IndexC,1) = C;conditions(IndexC,2) = C;
 

IndexC = strfind(conditions(:,1), 'THA');
IndexC = find(not(cellfun('isempty', IndexC)));
OUTPUT.stim(temp(IndexC),2) = OUTPUT.stim(temp(IndexC),2)+400;
conditions(IndexC,1) = num2cell(repmat(999,1,length(IndexC)));
C = cell(1, length(IndexC));C(:) = {'999'};conditions(IndexC,1) = C;conditions(IndexC,2) = C;
 

% visual
IndexC = strfind(conditions(:,2), 'BA');
IndexC = find(not(cellfun('isempty', IndexC)));
OUTPUT.stim(temp(IndexC),2) = OUTPUT.stim(temp(IndexC),2)+100;
 

IndexC = strfind(conditions(:,2), 'DA');
IndexC = find(not(cellfun('isempty', IndexC)));
OUTPUT.stim(temp(IndexC),2) = OUTPUT.stim(temp(IndexC),2)+200;
 

IndexC = strfind(conditions(:,2), 'TA');
IndexC = find(not(cellfun('isempty', IndexC)));
OUTPUT.stim(temp(IndexC),2) = OUTPUT.stim(temp(IndexC),2)+300;
 

IndexC = strfind(conditions(:,2), 'THA');
IndexC = find(not(cellfun('isempty', IndexC)));
OUTPUT.stim(temp(IndexC),2) = OUTPUT.stim(temp(IndexC),2)+400;
 

 

OUTPUT.Conditions = {'Aud','Cong','Incong','Vis-Alone','FaceOnset','A_BA','C_BA','I_BA','V_BA','A_DA','C_DA','I_DA','V_DA','A_TA','C_TA','I_TA','V_TA','A_THA','C_THA','I_THA','V_THA'};
OUTPUT.CondNums = [1:5 101:104 201:204 301:304 401:404];