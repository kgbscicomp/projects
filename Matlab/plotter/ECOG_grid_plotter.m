%% Find the number of unique electrode groups

grid_names = fieldnames(GridMaps.Real);


%% Plot all the grids
for i=1:length(grid_names)
    
    CurrGridName = cell2mat(grid_names(i));
    
    % Elec spacing
    EL_spacing = 5;
    EL_radius = 2;
    
    
    h(i) = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    
%     % Plot once with both, then subplot 2 with only lap
%     subplot(1,2,1)
    
    % draw orig electrodes first
    for xxx=1:size(GridMaps.Real.(genvarname(CurrGridName)),1) % DOWN
        for yyy=1:size(GridMaps.Real.(genvarname(CurrGridName)),2) % ACROSS
            
            % Plot circles
            tempx = ((yyy-1)*EL_spacing*2); % x and y flipped on matlab matrix
            tempy = ((xxx-1)*(-EL_spacing)*2); % x and y flipped on matlab matrix
            viscircles([ tempx tempy],EL_radius,'Color',[.8 .8 .8]);
            
            % Pull numbers from array and plot
            t = text(tempx,tempy,(chan_names_orig(GridMaps.Real.(genvarname(CurrGridName))(xxx,yyy))),'HorizontalAlignment','center');
            t.Color = [.7 .7 .7];
            t.FontSize = 10;

        end
    end
    
    % Lap
    counter = 1;
    for xxx=1:size(GridMaps.Lap.(genvarname(CurrGridName)),1) % DOWN
        for yyy=1:size(GridMaps.Lap.(genvarname(CurrGridName)),2) % ACROSS
            
            % Plot circles
            if GridMaps.Lap.(genvarname(CurrGridName))(xxx,yyy)>0
                tempx = ((yyy-1)*EL_spacing); % x and y flipped on matlab matrix
                tempy = ((xxx-1)*(-EL_spacing)); % x and y flipped on matlab matrix
                viscircles([ tempx tempy],EL_radius,'Color','k')
                
                % Pull numbers from array
                t = text(tempx,tempy,Virtual_El.(genvarname(CurrGridName)).Name(counter),'HorizontalAlignment','center');
                %t.Color = [.7 .7 .7];
                t.FontSize = 10;
                
                counter = counter+1;
            end
        end
    end
    axis equal;
%     hold off
%     
%     subplot(1,2,2)
%     hold on
%     
%     % Lap
%     counter = 1;
%     for xxx=1:size(GridMaps.Lap.(genvarname(CurrGridName)),1) % DOWN
%         for yyy=1:size(GridMaps.Lap.(genvarname(CurrGridName)),2) % ACROSS
%             
%             % Plot circles
%             if GridMaps.Lap.(genvarname(CurrGridName))(xxx,yyy)>0
%                 tempx = ((yyy-1)*EL_spacing); % x and y flipped on matlab matrix
%                 tempy = ((xxx-1)*(-EL_spacing)); % x and y flipped on matlab matrix
%                 viscircles([ tempx tempy],EL_radius,'Color','k')
%                 
%                 % Pull numbers from array
%                 t = text(tempx,tempy,Virtual_El.(genvarname(CurrGridName)).Name(counter),'HorizontalAlignment','center');
%                 %t.Color = [.7 .7 .7];
%                 t.FontSize = 10;
%                 
%                 counter = counter+1;
%             end
%         end
%     end
%     axis equal;
%     hold off
    
end