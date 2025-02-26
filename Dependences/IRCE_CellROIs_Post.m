function [Post_roi_corners, previous_ROI_count] = IRCE_CellROIs_Post(Ch1_corr_IRM, Save_individual_acq_dir, Add_more_ROIs_flag, IRM_thres, roi_corners, Data_FOVs_nums, manual_offset)
    % Last updated:
    % 20241204: KLS - Instead of showing all previous boxs, just show the
        % current ROI number
    % 20240822: KLS - Built in check for x/y positions outside
        % image dimensions
    if nargin < 7
        manual_offset = [0 0];
    end
    n = 2; % Setpoint two corners (Top left, bottom right)

    %---------------------------------------------------------%
    % View Min Projection of IRM
    %---------------------------------------------------------%
    img = min(Ch1_corr_IRM,[],3);
    
    imshow(img,[min(img,[],'all') prctile(img,99.5,'all')])
    g = gcf;
    g.WindowState = 'maximized';  
    
    median_mask = img < IRM_thres;
    
    [rows, cols] = size(median_mask);

    % Check if the endpoint data is a different number of 512 by 512 FOVs
    % compared to the assay data. If it is off set the old ROIs
    % check dimension 1 offset
    if Data_FOVs_nums(3) > Data_FOVs_nums(1)
        dim_1_offset = ((Data_FOVs_nums(3)-Data_FOVs_nums(1))*512)/2;
    elseif Data_FOVs_nums(3) < Data_FOVs_nums(1)
        disp('Endpoint data is smaller than assay data, that is not valid')
        return;
    else 
        dim_1_offset = 0;
    end
    % check dimension 2 offset
    if Data_FOVs_nums(4) > Data_FOVs_nums(2)
        dim_2_offset = ((Data_FOVs_nums(4)-Data_FOVs_nums(2))*512)/2;
    elseif Data_FOVs_nums(4) < Data_FOVs_nums(2)
        disp('Endpoint data is smaller than assay data, that is not valid')
        return;
    else
        dim_2_offset = 0;
    end
    % Add old ROI Boxs
    for i = 1:size(roi_corners,1)
        x = roi_corners{i,1}(:,1) + dim_1_offset + manual_offset(1);
        y = roi_corners{i,1}(:,2) + dim_2_offset + manual_offset(2);
    
        hold on
            rectangle('position',[min(x) min(y) max(x)-min(x) max(y)-min(y)],'EdgeColor','green','LineWidth',2,'LineStyle','--')
            text(mean(x),mean(y),[num2str(i,'%03.f')],'Color','green','FontSize',18,'HorizontalAlignment','center','FontName','Trebuchet MS')
        hold off
    end

    %---------------------------------------------------------%
    % Manually selection cell ROIs
    %---------------------------------------------------------%
    continueSettingROIs = 'k'; % Initial condition to enter the loop
    line_count = 0; % Keep track of how many lines have been drawn

    cd(Save_individual_acq_dir)
    if exist('Post_roi_corners.mat','file')
        load('Post_roi_corners','-mat')
    else
        % Initialize graphics handles
        hRectangle = [];
        hText = [];

        while line_count < size(roi_corners,1)
            % Delete previous rectangle and text if they exist
            if ~isempty(hRectangle)
                delete(hRectangle);
            end
            if ~isempty(hText)
                delete(hText);
            end

            % If there is a previous ROI, display it as a dashed box
            if line_count >= 1
                prev_x = Post_roi_corners{line_count,1}(:,1);
                prev_y = Post_roi_corners{line_count,1}(:,2);
                hold on
                    hRectangle = rectangle('position',[min(prev_x) min(prev_y) max(prev_x)-min(prev_x) max(prev_y)-min(prev_y)],'EdgeColor','green','LineWidth',2,'LineStyle','--');
                    hText = text(mean(prev_x),mean(prev_y),[num2str(line_count,'%03.f')],'Color','green','FontSize',18,'HorizontalAlignment','center','FontName','Trebuchet MS');
                hold off
            end

            line_count = line_count + 1;
            [x, y] = ginput(n);
            x = round(x);
            y = round(y);

            [x, y] = IRCE_ROIcorners_stay_within_image_bounds(Ch1_corr_IRM,x,y);

            Post_roi_corners{line_count,1} = [x,y]; 
            hold on
                rectangle('position',[min(x) min(y) max(x)-min(x) max(y)-min(y)],'EdgeColor','yellow','LineWidth',3)
                text(mean(x),mean(y),[num2str(line_count,'%03.f')],'Color','cyan','FontSize',18,'HorizontalAlignment','center','FontName','Trebuchet MS')
            hold off
        end
        save('Post_roi_corners.mat','Post_roi_corners','-v7.3')
    end
    
    if Add_more_ROIs_flag == 1
        previous_ROI_count = size(Post_roi_corners,1);
        line_count = size(Post_roi_corners,1);

        % Initialize graphics handles
        hRectangle = [];
        hText = [];

        while continueSettingROIs == 'k'
            % Delete previous rectangle and text if they exist
            if ~isempty(hRectangle)
                delete(hRectangle);
            end
            if ~isempty(hText)
                delete(hText);
            end

            % If there is a previous ROI, display it as a dashed box
            if line_count >= 1
                prev_x = Post_roi_corners{line_count,1}(:,1);
                prev_y = Post_roi_corners{line_count,1}(:,2);
                hold on
                    hRectangle = rectangle('position',[min(prev_x) min(prev_y) max(prev_x)-min(prev_x) max(prev_y)-min(prev_y)],'EdgeColor','green','LineWidth',2,'LineStyle','--');
                    hText = text(mean(prev_x),mean(prev_y),[num2str(line_count,'%03.f')],'Color','green','FontSize',18,'HorizontalAlignment','center','FontName','Trebuchet MS');
                hold off
            end

            line_count = line_count + 1;
            [x, y] = ginput(n);
            x = round(x);
            y = round(y);

            [x, y] = IRCE_ROIcorners_stay_within_image_bounds(Ch1_corr_IRM,x,y);

            Post_roi_corners{line_count,1} = [x,y];
            hold on
                rectangle('position',[min(x) min(y) max(x)-min(x) max(y)-min(y)],'EdgeColor','yellow','LineWidth',3)
                text(mean(x),mean(y),[num2str(line_count,'%03.f')],'Color','cyan','FontSize',18,'HorizontalAlignment','center','FontName','Trebuchet MS')
            hold off

            % Ask to continue
            disp('Add another ROI? (''k'' - yes, ''l'' - no):');
            continueSettingROIs = getkey(1);
            continueSettingROIs = lower(char(continueSettingROIs));

            % Loop until the input is 'k','l'
            while ~any(strcmp(continueSettingROIs, {'k', 'l'}))
                disp('Invalid input. Please enter ''k - yes'', ''l'' - no''.');
                continueSettingROIs = input('Add another ROI? (''k'' - yes, ''l'' - no):', 's');
            end
        end
        save('Post_roi_corners.mat','Post_roi_corners','-v7.3')
        
    else
        previous_ROI_count = size(Post_roi_corners,1);
    end

    close all
    
    %---------------------------------------------------------%
    % View Min Projection of IRM
    %---------------------------------------------------------%
    img = min(Ch1_corr_IRM,[],3);
    imshow(img,[min(img,[],'all') prctile(img,99.5,'all')])
    g = gcf;
    g.WindowState = 'maximized';  
    
    % Add ROI Boxes
    for i = 1:size(Post_roi_corners,1)
        x = Post_roi_corners{i,1}(:,1);
        y = Post_roi_corners{i,1}(:,2);

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
    exportgraphics(fig,'Post_Manually_Selected_ROIs.png','Resolution',300)

    clear g n continueSettingROIs line_count fig i Add_more_ROIs_flag x y hRectangle hText
end
