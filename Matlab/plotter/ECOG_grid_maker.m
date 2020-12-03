
function [dim_w,dim_h] = ECOG_grid_maker(CurrGridName)
% SIMPLE_GUI2 Select a data set from the pop-up menu, then
% click one of the plot-type push buttons. Clicking the button
% plots the selected data in the axes.

% dim_w = input(['Enter number of electrodes along the width of the grid ',CurrGridName,' (el numbers 1 2 3 etc.): ']);
% dim_h = input(['Enter number of electrodes along the height of the grid ',CurrGridName,' (next group of el): ']);
dim_w = input(['Enter number of electrodes along the width of the grid ',CurrGridName,': ']);
dim_h = input(['Enter number of electrodes along the height of the grid ',CurrGridName,': ']);
% dim_w = 5;
% dim_h = 4;
total_el = dim_w*dim_h;

% Expected corner clicks
% For 5 x 2, [1 5 6 10];
clicknumbers = [];
for j=[1 dim_h]
    clicknumbers = [clicknumbers (((j-1)*dim_w)+1)];
    clicknumbers = [clicknumbers (((j-1)*dim_w)+dim_w)];
end

% This will order the numbers in the following manner for 5 rows x 2 height
% 10 9  8  7  6
% 5  4  3  2  1

el_box_size = 50;
input_el = zeros(1,total_el);

global counter;
counter =1;


%  Create and then hide the UI as it is being constructed.
f = figure('Position',[0,0,dim_w*el_box_size,(dim_h*el_box_size)+50]);

% hsurf    = uicontrol('Style','pushbutton',...
%     'String','Surf','Position',[315,220,70,25],...
%     'Callback',@surfbutton_Callback);


for curr_w=1:dim_w
    for curr_h=1:dim_h
        k=curr_w+((curr_h-1)*dim_w);
        app.bit(k) = uicontrol('Style','pushbutton',...
            'String','','Position',[(curr_w-1)*el_box_size,(curr_h-1)*el_box_size,el_box_size,el_box_size],...
            'Callback',{@surfbutton_Callback,k});
    end
end

exit_bttn    = uicontrol('Style','pushbutton',...
    'String','Done','Position',[0,(curr_h)*el_box_size,el_box_size,el_box_size/2],...
    'Callback',{@exitbutton_Callback,input_el});

clear_bttn    = uicontrol('Style','pushbutton',...
    'String','Clear','Position',[el_box_size,(curr_h)*el_box_size,el_box_size,el_box_size/2],...
    'Callback',{@clearbutton_Callback,input_el});


% Initialize the UI.
% Change units to normalized so components resize automatically.
f.Units = 'normalized';
% %     hsurf.Units = 'normalized';

% Assign the a name to appear in the window title.
f.Name = 'Click the 4 corners of the grid from smallest to largest electrode';

disp('************************************************************************');
disp('************************************************************************');
disp('** Click the 4 corners of the grid from smallest to largest electrode **');
disp('************************************************************************');
disp('************************************************************************');

% Move the window to the center of the screen.
movegui(f,'center')

% Make the window visible.
f.Visible = 'on';

% Push button callbacks. Each callback plots current_data in the
% specified plot type.

% Once clicked, change to number clicked
% Have them input RAS into function, in order clicked.


    function [out] = surfbutton_Callback(source,eventdata,x)
        % Display surf plot of the currently selected data.
        %surf(current_data);
        %input_el(x) = 1;
        out = counter;
        set( source , 'UserData' , out )
        set(app.bit(x), 'BackgroundColor', [1,1,0])
        set(app.bit(x),'string',clicknumbers(counter),'enable','off');
        counter = counter+1;
    end
    
    function clearbutton_Callback(source,eventdata,input_el)
        counter = 1;
        for curr_w=1:dim_w
            for curr_h=1:dim_h
                k=curr_w+((curr_h-1)*dim_w);        
                app.bit(k) = uicontrol('Style','pushbutton',...
            'String','','Position',[(curr_w-1)*el_box_size,(curr_h-1)*el_box_size,el_box_size,el_box_size],...
            'Callback',{@surfbutton_Callback,k});
            end
        end
    end

    
    function exitbutton_Callback(source,eventdata,input_el)
        outputdata = zeros(1,total_el);
        % Display surf plot of the currently selected data.
        %surf(current_data);
        
        %set(app.bit(x), 'BackgroundColor', [1,1,0])
        
        for xxx=1:total_el
            myValues = get(app.bit(xxx),'UserData');
            temp = myValues;
            if ~isempty(temp)
                outputdata(xxx) = temp;
            end
        end
        if dim_w > 1 && dim_h > 1
            if counter == 5
                save([pwd '/Laplacian_Map_',CurrGridName,'.mat'],'outputdata')
            else
                error('Did not click the 4 corners only');
            end
        else
           if counter == 3
                save([pwd '/Laplacian_Map_',CurrGridName,'.mat'],'outputdata')
            else
                error('Did not click the 2 corners only');
           end
        end
        close(f)
    end
end
