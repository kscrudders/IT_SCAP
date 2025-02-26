function [roi_corners, previous_ROI_count] = IRCE_CellROIs(Ch1_corr_IRM, Save_individual_acq_dir, Add_more_ROIs_flag, IRM_thres, IRM_LUT)
    % Updated 20240822: KLS - Built in check for x/y positions outside
        % image dimensions
    n = 2; % Setpoint two corners (Top left, bottom right)

    %---------------------------------------------------------%
    % View Min Projection of IRM
    %---------------------------------------------------------%
    img = min(Ch1_corr_IRM,[],3);
    
    %imshow(img,[min(img,[],'all') prctile(img,99.5,'all')])
    imshow(img, IRM_LUT);
    g = gcf;
    g.WindowState = 'maximized';  
    
       %---------------------------------------------------------%
    % Let do a median projection to find stationary cells
    %---------------------------------------------------------%
    if size(Ch1_corr_IRM,3) > 15
        img = min(movmedian(Ch1_corr_IRM, 15, 3), [], 3);
    else
        img = min(movmedian(Ch1_corr_IRM, round(size(Ch1_corr_IRM,3)*.1)), [], 3);
    end
    median_mask = img < IRM_thres;
    
    [rows, cols] = size(median_mask);

    L = median_mask;
    %{
    [B,~] = bwboundaries(L>0,'noholes');
    hold on
    for k = 1:length(B)
        boundary = B{k};

        % Check if the boundary touches the edge of the image
        if any(boundary(:,1) == 1) || any(boundary(:,1) == rows) || any(boundary(:,2) == 1) || any(boundary(:,2) == cols)
            % Boundary touches the edge, color red
            plot(boundary(:,2), boundary(:,1), 'Color', 'red', 'LineStyle', ':', 'LineWidth', 2)
        else
            % Boundary does not touch the edge, color yellow
            plot(boundary(:,2), boundary(:,1), 'Color', 'blue', 'LineStyle', ':', 'LineWidth', 2)
        end
    end
    %}

    % Get boundaries from the binary image
    [B, ~] = bwboundaries(L > 0, 'noholes');
    
    % Initialize arrays to hold all x and y coordinates for red and blue boundaries
    x_all_red = [];
    y_all_red = [];
    x_all_blue = [];
    y_all_blue = [];
    
    % Loop through each boundary
    for k = 1:length(B)
        boundary = B{k};
    
        % Check if the boundary touches the edge of the image
        if any(boundary(:,1) == 1) || any(boundary(:,1) == rows) || ...
           any(boundary(:,2) == 1) || any(boundary(:,2) == cols)
            % Boundary touches the edge, add to red arrays
            x_all_red = [x_all_red; boundary(:,2); NaN];
            y_all_red = [y_all_red; boundary(:,1); NaN];
        else
            % Boundary does not touch the edge, add to blue arrays
            x_all_blue = [x_all_blue; boundary(:,2); NaN];
            y_all_blue = [y_all_blue; boundary(:,1); NaN];
        end
    end

    hold on;
    
    % Plot all red boundaries in one call
    if ~isempty(x_all_red)
        plot(x_all_red, y_all_red, 'Color', 'red', 'LineStyle', ':', 'LineWidth', 1);
    end
    
    % Plot all blue boundaries in one call
    if ~isempty(x_all_blue)
        plot(x_all_blue, y_all_blue, 'Color', 'blue', 'LineStyle', ':', 'LineWidth', 1);
    end
    
    %---------------------------------------------------------%
    % Manually selection cell ROIs
    %---------------------------------------------------------%
    continueSettingROIs = 'k'; % Initial condition to enter the loop
    line_count = 0; % Keep track of how many lines have been drawn

    cd(Save_individual_acq_dir)
    if exist('roi_corners.mat','file')
        load('roi_corners','-mat')
    else
        while continueSettingROIs == 'k'
            line_count = line_count + 1;
            [x, y] = ginput(n);
            x = round(x);
            y = round(y);

            [x, y] = IRCE_ROIcorners_stay_within_image_bounds(Ch1_corr_IRM,x,y);

            roi_corners{line_count,1} = [x,y]; 
            hold on
                rectangle('position',[min(x) min(y) max(x)-min(x) max(y)-min(y)],'EdgeColor','yellow','LineWidth',3)
                text(mean(x),mean(y),[num2str(line_count,'%03.f')],'Color','cyan','FontSize',18,'HorizontalAlignment','center','FontName','Trebuchet MS')
            hold off
            
            disp('Add another ROI? (''k'' - yes, ''l'' - no):');
            continueSettingROIs = getkey(1);
            continueSettingROIs = lower(char(continueSettingROIs));
            
            % Loop until the input is 'k','l'
            while ~any(strcmp(continueSettingROIs, {'k', 'l'}))
                disp('Invalid input. Please enter ''k - yes'', ''l - no''.');
                continueSettingROIs = input('Add another ROI? (''k'' - yes, ''l'' - no):', 's');
            end
        end
        save('roi_corners.mat','roi_corners','-v7.3')
    end
    
    if Add_more_ROIs_flag == 1
        previous_ROI_count = size(roi_corners,1);
        line_count = size(roi_corners,1);

        % Add ROI Boxs
        for i = 1:size(roi_corners,1)
            x = roi_corners{i,1}(:,1);
            y = roi_corners{i,1}(:,2);

            hold on
                rectangle('position',[min(x) min(y) max(x)-min(x) max(y)-min(y)],'EdgeColor','yellow','LineWidth',2)
            hold off
        end

        while continueSettingROIs == 'k'
            line_count = line_count + 1;
            [x, y] = ginput(n);
            x = round(x);
            y = round(y);

            [x, y] = IRCE_ROIcorners_stay_within_image_bounds(Ch1_corr_IRM,x,y);

            roi_corners{line_count,1} = [x,y];
            hold on
                rectangle('position',[min(x) min(y) max(x)-min(x) max(y)-min(y)],'EdgeColor','yellow','LineWidth',3)
                text(mean(x),mean(y),[num2str(line_count,'%03.f')],'Color','cyan','FontSize',18,'HorizontalAlignment','center','FontName','Trebuchet MS')
            hold off
            

            disp('Add another ROI? (''k'' - yes, ''l'' - no):');
            continueSettingROIs = getkey(1);
            continueSettingROIs = lower(char(continueSettingROIs));
            
            % Loop until the input is 'k','l'
            while ~any(strcmp(continueSettingROIs, {'k', 'l'}))
                disp('Invalid input. Please enter ''k - yes'', ''l - no''.');
                continueSettingROIs = input('Add another ROI? (''k'' - yes, ''l'' - no):', 's');
            end
        end
        save('roi_corners.mat','roi_corners','-v7.3')
        
    else
        previous_ROI_count = size(roi_corners,1);
    end

    close all
    
    %---------------------------------------------------------%
    % View Min Projection of IRM
    %---------------------------------------------------------%
    img = min(Ch1_corr_IRM,[],3);
    %imshow(img,[min(img,[],'all') prctile(img,99.5,'all')])
    imshow(img, IRM_LUT);
    g = gcf;
    g.WindowState = 'maximized';  
    
    % Add ROI Boxs
    for i = 1:size(roi_corners,1)
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);

        hold on
            rectangle('position',[min(x) min(y) max(x)-min(x) max(y)-min(y)],'EdgeColor','yellow','LineWidth',2)
            text(mean(x),mean(y),[num2str(i,'%03.f')],'Color','cyan','FontSize',18,'HorizontalAlignment','center','FontName','Trebuchet MS')
        hold off
    end

    % Make the last ROI number red to easily see what is the maximum ROI number
    hold on
        text(mean(x),mean(y),num2str(i,'%03.f'),'Color','red','FontSize',18,'HorizontalAlignment','center','FontName','Trebuchet MS')
    hold off

    fig = gcf;
    exportgraphics(fig,'Manually_Selected_ROIs.png','Resolution',300)

    clear g n continueSettingROIs line_count fig i Add_more_ROIs_flag x y
end