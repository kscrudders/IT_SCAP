function [roi_corners, previous_ROI_count] = IRCE_SingeFrame_CellROIs(img_stack, Save_individual_acq_dir, Add_more_ROIs_flag)
    % Updated 20240822: KLS - Built in check for x/y positions outside
        % image dimensions
    n = 2; % Setpoint two corners (Top left, bottom right)

    %---------------------------------------------------------%
    % View Min Projection of IRM
    %---------------------------------------------------------%
    img = img_stack(:,:,1);
    
    imshow(img,[min(img,[],'all') prctile(img,99.5,'all')])
    g = gcf;
    g.WindowState = 'maximized';  
    
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

            [x, y] = IRCE_ROIcorners_stay_within_image_bounds(img_stack,x,y);

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

            [x, y] = IRCE_ROIcorners_stay_within_image_bounds(img_stack,x,y);

            roi_corners{line_count,1} = [x,y]; %#ok<SAGROW>
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
    img = img_stack(:,:,1);
    imshow(img,[min(img,[],'all') prctile(img,99.5,'all')])
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
    sound(exp(sin(1:1500))) 
end