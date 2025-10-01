function roi_corners = IRCE_Remove_CellROIs(Ch1_corr_IRM, Save_individual_acq_dir, roi_corners, Specific_ROIs)
    cd(Save_individual_acq_dir)
    if exist('roi_corners.mat','file')
        load('roi_corners','-mat')
    else
        disp('You need to generate ROIs corners before you edit them.')
        return;
    end
    if bitor(max(Specific_ROIs) > length(roi_corners), min(Specific_ROIs) < 1)
        disp(['Specific_ROIs to be edited must be between 1 and ' num2str(length(roi_corners))])
        return;
    end
    
    Rois_to_keep = setdiff(1:length(roi_corners),Specific_ROIs);
    
    %---------------------------------------------------------%
    % Manually remove cell ROI selection(s) 
    %---------------------------------------------------------%
    old_roi_corners = roi_corners;
    roi_corners = cell([length(Rois_to_keep) 1]);
    for n = 1:length(Rois_to_keep)
        i = Rois_to_keep(n);
        
        roi_corners{n,1}(:,1) = old_roi_corners{i,1}(:,1);
        roi_corners{n,1}(:,2) = old_roi_corners{i,1}(:,2);        
    end
    
    save('roi_corners.mat','roi_corners','-v7.3')

    close all
    
    %---------------------------------------------------------%
    % View Min Projection of IRM
    %---------------------------------------------------------%
    img = min(Ch1_corr_IRM,[],3);
    imshow(img,[min(img,[],'all') prctile(img,99.5,'all')])
    g = gcf;
    g.WindowState = 'maximized';  
    
    % Add ROI Boxs
    for i = 1:size(roi_corners,1)
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);

        hold on
            rectangle('position',[min(x) min(y) max(x)-min(x) max(y)-min(y)],'EdgeColor','yellow','LineWidth',2)
        hold off
    end
    
    % Add ROI Numbers as Centered Text
    for i = 1:size(roi_corners,1)
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);
        hold on
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