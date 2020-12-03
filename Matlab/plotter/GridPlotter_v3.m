% Plot data in grid of circles
% function GridPlotter(AvgTime,MinMaxScale,condlist,GridOrder)
% AvgTime = [0 .5]; % in seconds, time to average data
% MinMaxScale = [-3 3];

% May 8 2017
% Updated for new scripts

condlist = SPECS.plot_conds(:,1)';


if OUTPUT.reference==2 % Modify for laplacian
    
    for xxx = 1:length(condlist)
        currcond = condlist(xxx);
        
        [c TimeSmpl1] = min(abs(time-AvgTime(1)));
        [c TimeSmpl2] = min(abs(time-AvgTime(2)));
        if strcmp(GridStat,'ERSP-line')
            currdata = AVG.(genvarname(['ERSP' num2str(currcond) '_avg']));
        elseif strcmp(GridStat,'ERSPP')
            currdata = AVG.ERSPP;
            currdata = 1-abs(currdata);
        elseif strcmp(GridStat,'ERSPpvalues')
            currdata = AVG.ERSPpvalues;
            currdata = 1-currdata;
        else
            currdata = AVG.(genvarname([GridStat]));
        end
        
        GridGroup = regexp(chan_names,'\D*','match');
        ChanNumbers = regexp(chan_names,'\d*','match');
        temp = [];
        tempb = {};
        for i=1:size(ChanNumbers,2)
            temp = [temp str2double(ChanNumbers{1,i})];
            tempb = [tempb GridGroup{1,i}{1}];
        end
        ChanNumbers = temp;
        GridGroup = tempb;
        temp = strcmp(GridGroup,'A');
        ChanGroup = [];
        ChanGroup(temp) = 1;
        if length(ChanNumbers)~=sum(ChanGroup)
            temp = strcmp(GridGroup,'B');
            temp = find(temp);
            ChanGroup(temp) = 2;
        end
        
        Grid.A1 = find(ChanGroup==1);
        Grid.A2 = find(ChanGroup==2);
        
        % size of each small box. Final image will be 5Nx4N: 5 down x 4 across
        N = 100;
        
        % Arrangement of grid: 5 x 4, 1 is on top right 5 on bottom right
        SizeDown = 9;
        SizeAcross = 7;
        if grid_arrangement==0
            GridArrgment = [1:7;8:14;15:21;22:28;29:35;36:42;43:49;50:56;57:63];
        elseif grid_arrangement==1
            GridArrgment = fliplr([1:9;10:18;19:27;28:36;37:45;46:54;55:63])';
        end
        
        % create a circle mask
        t = linspace(0,2*pi,50);   % approximated by 100 lines
        r = (N-10)/2;              % circles will be separated by a 10 pixels border
        emptycircle = poly2mask(r*cos(t)+N/2+0.5, r*sin(t)+N/2+0.5, N, N);
        
        
        % Create zero matrix for each grid
        Grid.Image1 = zeros(N*SizeDown,N*SizeAcross);
        Grid.Image2 = zeros(N*SizeDown,N*SizeAcross);
        Grid.Image1alpha = zeros(N*SizeDown,N*SizeAcross);
        Grid.Image2alpha = zeros(N*SizeDown,N*SizeAcross);
        
        figure
        
        % need to create array to adjust the positioning of the electrodes
        % from the 4x5 grid to the 7x9 laplacian grid
        
        lap_grid = zeros(1,7*9);
        temp = 0;
        counter = 1;
        for i=1:length(lap_grid)
            if temp==1
            lap_grid(i) = counter;
            counter = counter+1;
            temp = 0;
            else
                temp = 1;
            end
        end
        
        
        % Run for loop for each grid with data
        for j=1:2
            currentgrid = Grid.(genvarname(['A' num2str(j)]));
            for i=currentgrid
                
                % current electrode
                el = i;% + ;
                % gridel = ChanNumbers(i) - ((SizeDown*SizeAcross)*(j-1));
                gridel = ChanNumbers(i);
                
                % fill in data from ERSP
                datacircle = emptycircle*mean(currdata(el,TimeSmpl1:TimeSmpl2));
                
                % find the location on the grid for circle
                [x,y] = find(GridArrgment==find(gridel==lap_grid));
                currX = ((x-1)*N)+1:((x-1)*N)+N;
                currY = ((y-1)*N)+1:((y-1)*N)+N;
                
                Grid.(genvarname(['Image' num2str(j)]))(currX,currY) = datacircle;
                
            end
            subplot(str2num(strcat(num2str(12),num2str(j))))
            h = imagesc(Grid.(genvarname(['Image' num2str(j)])));
            caxis([MinMaxScale(1) MinMaxScale(2)])
            
            % Make the background clear
            [ tempy] = find(Grid.(genvarname(['Image' num2str(j)])));
            Grid.(genvarname(['Image' num2str(j) 'alpha']))(Grid.(genvarname(['Image' num2str(j)]))~=0)=1;
                set(h, 'AlphaData', logical(Grid.(genvarname(['Image' num2str(j) 'alpha']))));
% % % % %             
% % % % %             % Plot black circles where no data exists
% % % % %             hold on
% % % % %             currentgrid = ChanNumbers(Grid.(genvarname(['A' num2str(j)]))) - ((SizeDown*SizeAcross)*(j-1));
% % % % %             emptygrid = (1:(SizeDown*SizeAcross))+((SizeDown*SizeAcross)*(j-1));
% % % % %             %emptygrid(currentgrid) = [];
% % % % %             for i=1:length(emptygrid)
% % % % %                 
% % % % %                 currel = emptygrid(i);
% % % % %                 
% % % % %                 % Keep or delete electrode
% % % % %                                 
% % % % %                 if ismember(emptygrid(i)-((SizeDown*SizeAcross)*(j-1)),currentgrid)
% % % % %                     el_color = [1 1 1];
% % % % %                 else
% % % % %                     el_color = [0 0 0];
% % % % %                 end
% % % % %                 % current electrode
% % % % %                 el = currel;% + ;
% % % % %                 gridel = currel - ((SizeDown*SizeAcross)*(j-1));
% % % % %                 
% % % % %                 % find the location on the grid for circle
% % % % %                 [x,y] = find(GridArrgment==find(gridel==lap_grid));
% % % % %                 % [y,x] = find(GridArrgment==gridel);
% % % % %                 currX = ((x-1)*N)+(N/2);
% % % % %                 currY = ((y-1)*N)+(N/2);
% % % % %                 viscircles([currX currY], N/2,'EdgeColor',el_color,'LineWidth',10);
% % % % %             end
            hold off
        end
        suptitle(Conditions(currcond))
    end
else
    
    for xxx = 1:length(condlist)
        currcond = condlist(xxx);
        
        [c TimeSmpl1] = min(abs(time-AvgTime(1)));
        [c TimeSmpl2] = min(abs(time-AvgTime(2)));
        if strcmp(GridStat,'ERSP-line')
            currdata = AVG.(genvarname(['ERSP' num2str(currcond) '_avg']));
        elseif strcmp(GridStat,'ERSPP')
            currdata = AVG.ERSPP;
            currdata = 1-abs(currdata);
        elseif strcmp(GridStat,'ERSPpvalues')
            currdata = AVG.ERSPpvalues;
            currdata = 1-currdata;
        else
            currdata = AVG.(genvarname([GridStat]));
        end
        
        
        if isnumeric(chan_names)
            ChanNumbers = str2double(chan_names);
        else
            ChanNumbers = [];
            k = strfind(chan_names,'A');
            newStr = erase(chan_names,'A');
            %index = [];
            for yyy = 1:length(k)
                if cell2mat(k(yyy))==1
                    %index = [index yyy];
                    ChanNumbers = [ChanNumbers str2double(newStr(yyy))];
                end
            end
            
            k = strfind(chan_names,'B');
            newStr = erase(chan_names,'B');
            %index = [];
            for yyy = 1:length(k)
                if cell2mat(k(yyy))==1
                    %index = [index yyy];
                    ChanNumbers = [ChanNumbers str2double(newStr(yyy))+20];
                end
            end
        end
        if GridOrder==1
            Grid.A1 = find(ChanNumbers<=20);
            Grid.A2 = find(ChanNumbers>20);
        elseif GridOrder==2
            Grid.A2 = find(ChanNumbers<=20);
            Grid.A1 = find(ChanNumbers>20);
        end
        
        
        
%         
%         if isnumeric(chan_names)
%             ChanNumbers = str2double(chan_names);
%             if GridOrder==1
%                 Grid.A1 = find(ChanNumbers<=20);
%                 Grid.A2 = find(ChanNumbers>20);
%             elseif GridOrder==2
%                 Grid.A2 = find(ChanNumbers<=20);
%                 Grid.A1 = find(ChanNumbers>20);
%             end
%         else
%             k = strfind(chan_names,'A');
%             newStr = erase(chan_names,'A');
%             ChanNumbers = cellfun(@(x)str2double(x), newStr);
%             Grid.A1 = find(ChanNumbers<=20);
%             
%             k = strfind(chan_names,'B');
%             index = find([k{:}] == 1);
%             temp = chan_names(index);
%             k(~cellfun('isempty',k))
%             newStr = erase(chan_names,'B');
%             temp = cellfun(@(x)str2double(x), newStr);
%             ChanNumbers = [ChanNumbers temp+20];
%             Grid.A2 = find(ChanNumbers>20);
%             
%         end
%         
        % size of each small box. Final image will be 5Nx4N: 5 down x 4 across
        N = 100;
        
        % Arrangement of grid: 5 x 4, 1 is on top right 5 on bottom right
        SizeDown = 5;
        SizeAcross = 4;
        GridArrgment = fliplr([1:5;6:10;11:15;16:20]');
        
        
        % create a circle mask
        t = linspace(0,2*pi,50);   % approximated by 100 lines
        r = (N-10)/2;              % circles will be separated by a 10 pixels border
        emptycircle = poly2mask(r*cos(t)+N/2+0.5, r*sin(t)+N/2+0.5, N, N);
        
        
        % Create zero matrix for each grid
        Grid.Image1 = zeros(N*SizeDown,N*SizeAcross);
        Grid.Image2 = zeros(N*SizeDown,N*SizeAcross);
        Grid.Image1alpha = zeros(N*SizeDown,N*SizeAcross);
        Grid.Image2alpha = zeros(N*SizeDown,N*SizeAcross);
        
        figure
        
        % Run for loop for each grid with data
        for j=1:2
            currentgrid = Grid.(genvarname(['A' num2str(j)]));
            for i=currentgrid
                
                % current electrode
                el = i;% + ;
                %gridel = ChanNumbers(i) - ((SizeDown*SizeAcross)*(j-1));
                gridel = ChanNumbers(i);
                if gridel>20,gridel=gridel-20;end
                
                % fill in data from ERSP
                datacircle = emptycircle*mean(currdata(el,TimeSmpl1:TimeSmpl2));
                
                % find the location on the grid for circle
                [x,y] = find(GridArrgment==gridel);
                currX = ((x-1)*N)+1:((x-1)*N)+N;
                currY = ((y-1)*N)+1:((y-1)*N)+N;
                
                Grid.(genvarname(['Image' num2str(j)]))(currX,currY) = datacircle;
                
            end
            subplot(str2num(strcat(num2str(12),num2str(j))))
            h = imagesc(Grid.(genvarname(['Image' num2str(j)])));
            caxis([MinMaxScale(1) MinMaxScale(2)])
            
            % Make the background clear
            [ tempy] = find(Grid.(genvarname(['Image' num2str(j)])));
            Grid.(genvarname(['Image' num2str(j) 'alpha']))(Grid.(genvarname(['Image' num2str(j)]))~=0)=1;
            %     set(h, 'AlphaData', logical(Grid.(genvarname(['Image' num2str(j) 'alpha']))));
            
            % Plot black circles where no data exists
            hold on
            currentgrid = ChanNumbers(Grid.(genvarname(['A' num2str(j)]))) - ((SizeDown*SizeAcross)*(j-1));
            emptygrid = (1:(SizeDown*SizeAcross))+((SizeDown*SizeAcross)*(j-1));
            %emptygrid(currentgrid) = [];
            for i=1:length(emptygrid)
                
                currel = emptygrid(i);
                
                % Keep or delete electrode
                if ismember(emptygrid(i)-((SizeDown*SizeAcross)*(j-1)),currentgrid)
                    el_color = [1 1 1];
                else
                    el_color = [0 0 0];
                end
                % current electrode
                el = currel;% + ;
                gridel = currel - ((SizeDown*SizeAcross)*(j-1));
                
                % find the location on the grid for circle
                [y,x] = find(GridArrgment==gridel);
                currX = ((x-1)*N)+(N/2);
                currY = ((y-1)*N)+(N/2);
                viscircles([currX currY], N/2,'EdgeColor',el_color,'LineWidth',10);
            end
            hold off
        end
        suptitle(Conditions(xxx))
    end
    % end
    
end







