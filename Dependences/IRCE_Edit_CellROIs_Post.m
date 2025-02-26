function Post_roi_corners = IRCE_Edit_CellROIs_Post(Ch1_corr_IRM, Save_individual_acq_dir, IRM_thres, Post_roi_corners, Specific_ROIs)
    cd(Save_individual_acq_dir)
    if exist('roi_corners.mat','file')
        load('roi_corners','-mat')
    else
        disp('You need to generate ROIs corners before you edit them.')
        return;
    end
    if bitor(max(Specific_ROIs) > length(Post_roi_corners), min(Specific_ROIs) < 1)
        disp(['Specific_ROIs to be edited must be between 1 and ' num2str(length(Post_roi_corners))])
        return;
    end

    %---------------------------------------------------------%
    % View Min Projection of IRM
    %---------------------------------------------------------%
    img = min(Ch1_corr_IRM,[],3);
    
    imshow(img,[min(img,[],'all') prctile(img,99.5,'all')])
    g = gcf;
    g.WindowState = 'maximized';  
  
    %---------------------------------------------------------%
    % Manually edit cell ROI selection(s) 
    %---------------------------------------------------------%
    for i = Specific_ROIs
        
        x = Post_roi_corners{i,1}(:,1);
        y = Post_roi_corners{i,1}(:,2);        
        
        hold on
            rectangle('position',[min(x) min(y) max(x)-min(x) max(y)-min(y)],'EdgeColor','blue','LineWidth',2,'LineStyle','--')
            text(mean(x),mean(y),num2str(i,'%03.f'),'Color','blue','FontSize',18,'HorizontalAlignment','center','FontName','Trebuchet MS')
        hold off

        [x, y] = ginput(2);
        x = round(x);
        y = round(y);

        Post_roi_corners{i,1} = [x,y]; 
        hold on
            rectangle('position',[min(x) min(y) max(x)-min(x) max(y)-min(y)],'EdgeColor','yellow','LineWidth',3)
        hold off
    end
    
    save('Post_roi_corners.mat','Post_roi_corners','-v7.3')

    close all
    
    %---------------------------------------------------------%
    % View Min Projection of IRM
    %---------------------------------------------------------%
    img = min(Ch1_corr_IRM,[],3);
    imshow(img,[min(img,[],'all') prctile(img,99.5,'all')])
    g = gcf;
    g.WindowState = 'maximized';  
    
    % Add ROI Boxs
    for i = 1:size(Post_roi_corners,1)
        x = Post_roi_corners{i,1}(:,1);
        y = Post_roi_corners{i,1}(:,2);

        hold on
            rectangle('position',[min(x) min(y) max(x)-min(x) max(y)-min(y)],'EdgeColor','yellow','LineWidth',2)
        hold off
    end
    
    % Add ROI Numbers as Centered Text
    for i = 1:size(Post_roi_corners,1)
        x = Post_roi_corners{i,1}(:,1);
        y = Post_roi_corners{i,1}(:,2);
        hold on
            text(mean(x),mean(y),[num2str(i,'%03.f')],'Color','cyan','FontSize',18,'HorizontalAlignment','center','FontName','Trebuchet MS')
        hold off
    end
    
    % Make the last ROI number red to easily see what is the maximum ROI number
    hold on
        text(mean(x),mean(y),num2str(i,'%03.f'),'Color','red','FontSize',18,'HorizontalAlignment','center','FontName','Trebuchet MS')
    hold off

    fig = gcf;
    exportgraphics(fig,'Post_Manually_Selected_ROIs.png','Resolution',300)

    clear g n continueSettingROIs line_count fig i Add_more_ROIs_flag x y
    sound(exp(sin(1:1500))) 
end