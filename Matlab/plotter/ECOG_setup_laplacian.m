% Modify the data structure to create laplacian referenced virtual
% electrodes. Need to first report the arrangement of the electrodes, then
% virtually fill in all the gaps.

% Use the simple_gui3.m as a backbone.
% Input the number of rows and and columns. Click the corners in order (4
% clicks) or 2 clicks if 1xN and extrapolate the remaining numbers

% Interpolate the laplacian referenced activity between the rows and
% colums, which will actually yield a larger number of virtual electrodes
% than real electrodes. It's important to apply this procedure uniformly in
% the localizer and task-of-interest.

%% Find the number of unique electrode groups

grid_names = {};
for i=1:size(chan_names,2)
    temp = isletter(char((chan_names(i))));
    temp2 = cell2mat(chan_names(i));
    grid_names{i} = temp2(temp);
end

[C,~,IC] = unique(grid_names);
grid_names = C;
grid_groups = IC;

%% Iterate through each group to prompt about the number of rows and columns, then click the corners and extrapolate

% Find the number of electrodes in this group in the orig dataset
orig_ch_names = regexp(chan_names_orig,'\D*','match');
orig_ch_numbers = regexp(chan_names_orig,'\d*','match');
temparray = zeros(1,length(orig_ch_numbers));
for xxx=1:length(orig_ch_numbers)
    try
        temparray(xxx) = str2double(cell2mat(orig_ch_numbers{xxx}));
    catch
    end
end
orig_ch_numbers = temparray;

for i=1:length(grid_names)
    
    CurrGridName = cell2mat(grid_names(i));
    
    
    % At some point add in a check if the map already exists. If it does,
    % load it and update the channels
    
    
    [dim_w,dim_h] = ECOG_grid_maker(CurrGridName);
    
    % Find the number of electrodes in this group in the orig dataset
    grid_el_counter = 0;
    grid_el_numbers = []; % el number in orig_ch_label
    for xxx=1:length(orig_ch_names)
        if strcmp(cell2mat(orig_ch_names{xxx}),CurrGridName)
            grid_el_counter = grid_el_counter+1;
            grid_el_numbers = [grid_el_numbers xxx];
        end
    end
    
    % Find contact numbers within grid
    % Usually in order but to prevent weirdness if they needed to be
    % renumbered after data was collected (e.g., tech plugged in wrong
    % cable)
    grid_contact_numbers = orig_ch_numbers(grid_el_numbers);
    
    % Check number of electrodes matches?
    % size(find(grid_groups==i),1);
    while grid_el_counter ~= dim_w*dim_h
        disp('Size of Grid does not match')
        [dim_w,dim_h] = ECOG_grid_maker(CurrGridName);
        
    end
    
    % load in mat
    load([pwd '/Laplacian_Map_',CurrGridName,'.mat']);

    
    % Order electrodes according to the given map
    % determine if left to right or right to left
    % determine if bottom to top or top to bottom
    
    % Gui starts in the lower left corner and moves right
    
    El_map = zeros(dim_h,dim_w);
    
    if dim_h==1 || dim_w==1
        num_corners = 2;
    else
        num_corners = 4;
    end
    
    EL_locs = zeros(3,num_corners);
    for xxx=1:num_corners
        EL_locs(1,xxx) = find(outputdata==xxx);
        
        % find row
        EL_locs(2,xxx) = dim_h - ceil(EL_locs(1,xxx)/dim_w)+1;
        
        % find col
        EL_locs(3,xxx) = EL_locs(1,xxx) - ((ceil(EL_locs(1,xxx)/dim_w)-1)*dim_w);
        
    end
    
    % add row and col numbers under each EL_locs el
    if EL_locs(1,1)>EL_locs(1,2)
        if num_corners==1
            for xxx=1
                temp =  EL_locs(1,xxx);
                EL_locs(1,xxx) =  EL_locs(1,xxx+1);
                EL_locs(1,xxx+1) =  temp;
            end
        else
            for xxx=[1 dim_h-1]
                temp =  EL_locs(1,xxx);
                EL_locs(1,xxx) =  EL_locs(1,xxx+1);
                EL_locs(1,xxx+1) =  temp;
            end
        end
    end
    %leftright = 1;else leftright = 0;end
    
    for xxx=1:num_corners
        El_map(EL_locs(2,xxx),EL_locs(3,xxx)) = EL_locs(1,xxx);
    end
    
    % Fill in 4 corners
    El_map(1,:) = linspace(El_map(1,1),El_map(1,end),dim_w); % Top Row
    El_map(end,:) = linspace(El_map(end,1),El_map(end,end),dim_w); % Bottom Row
    for xxx=1:dim_w
        El_map(:,xxx) = linspace(El_map(1,xxx),El_map(end,xxx),dim_h); % Bottom Row
    end

    % Prob best to initially use the chan_names_orig var then corr later
    % Actually, I think this will be run once. Can use chan_names

    % Find bad chans in orig array
    Grid_bad_chans = intersect(grid_el_numbers,input_bad_chans);
    grid_el_numbers_good = grid_el_numbers;
    grid_el_numbers_good(Grid_bad_chans) = -1;

    % Replace map numbers with orig chan numbers (matches EEG data index)
    GridMaps.Real.(genvarname(CurrGridName)) = grid_el_numbers_good(El_map);
        
    % Create virtual laplacian referenced electrodes
    % Need to make sure at some stage that we're accounting for bad_chans
    % chan_lbls is our map the channel index in data
    % We want to exclude channels that are bad from Laplacian
    
    % Create laplacian data structure with extra dimension for what
    % difference between electrodes to take
    % Size is ((n*2)-1) x ((m*2)-1) - (n x m)
    % for 2 x 4 grid, it'd be 
    num_lap_el = (((dim_w*2)-1) * ((dim_h*2)-1)) - (dim_w * dim_h) - ((dim_w-1)*(dim_h-1));
    Lap_array_temp = zeros((dim_h*2)-1,(dim_w*2)-1);
    
    % Find coords in Lap array that need to be filled in
    % In
    El_counter = -1;
    for xxx=1:size(Lap_array_temp,1)
        for yyy=1:size(Lap_array_temp,2)
            Lap_array_temp(xxx,yyy) = El_counter;
            El_counter = El_counter*-1;
        end
    end
    %imagesc(Lap_array_temp)
    
    % input the real data electrodes (after bad chans removed)
    % Map to GridMaps.Real locations
    x_counter = 1;
    y_counter = 1;
    for xxx=1:2:size(Lap_array_temp,1)
        for yyy=1:2:size(Lap_array_temp,2)          
            Lap_array_temp(xxx,yyy) = GridMaps.Real.(genvarname(CurrGridName))(x_counter,y_counter);
            y_counter = y_counter+1;
        end
        x_counter = x_counter+1;
        y_counter = 1;
    end
    
    % Create Laplacian mapping by adding extra dimension
    GridMaps.Lap.(genvarname(CurrGridName)) = zeros((dim_h*2)-1,(dim_w*2)-1,2);
    for xxx=1:size(Lap_array_temp,1)
        for yyy=1:size(Lap_array_temp,2)
            if Lap_array_temp(xxx,yyy)==1
                % if column is even, sub (yyy+1) - (yyy-1)
                % if column is odd, sub (xxx+1) - (xxx-1)
                if mod(yyy,2)==0 %if 0, then even, left - right
                    temp1 = Lap_array_temp(xxx,yyy-1);
                    temp2 = Lap_array_temp(xxx,yyy+1);
                    if min([temp1 temp2])>0 % Check if value is one. if it is, run calculation. Else leave 0.
                        GridMaps.Lap.(genvarname(CurrGridName))(xxx,yyy,1) = temp1;
                        GridMaps.Lap.(genvarname(CurrGridName))(xxx,yyy,2) = temp2;
                    end
                else %then odd, top - bottom
                    temp1 = Lap_array_temp(xxx-1,yyy);
                    temp2 = Lap_array_temp(xxx+1,yyy);
                    if min([temp1 temp2])>0 % Check if value is one. if it is, run calculation. Else leave 0.
                        GridMaps.Lap.(genvarname(CurrGridName))(xxx,yyy,1) = temp1;
                        GridMaps.Lap.(genvarname(CurrGridName))(xxx,yyy,2) = temp2;
                    end
                end
            end
        end
    end
    
%     % Create a simplified map of only the Lap electrodes (including bad)
%     x_counter = 1;
%     y_counter = 1;
%     GridMaps.LapCorrds.(genvarname(CurrGridName)) = [];
%     for xxx=1:2:size(Lap_array_temp,1)
%         for yyy=1:2:size(Lap_array_temp,2)
%             GridMaps.LapCorrds.(genvarname(CurrGridName))(x_counter,y_counter,:) = GridMaps.Lap.(genvarname(CurrGridName))(xxx,yyy,:);
%             y_counter = y_counter+1;
%         end
%         x_counter = x_counter+1;
%         y_counter = 1;
%     end
    
    
    
    
end

% Merge maps together and delete the multiple .mats output by grid_maker







%% Plot modified grids with newly numbered electrodes on the grid to be saved/printed
% Also plot the original names of the elelectrodes on a matching map for
% easy comparison




%% Save the laplacian maps to be loaded next time







