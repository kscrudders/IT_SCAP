function [Img_ROIs_out, roi_corners] = IRCE_MaskAndROICorrnersCrop(Dilated_label_ROIs, Channel_stack, roi_corners)
    Img_ROIs_out = cell([size(Dilated_label_ROIs,1) 1]);    
    
    pix_pad = 11;

    max_ROI = size(Dilated_label_ROIs,1);
    ROI_n = 1;
    while ROI_n <= max_ROI
        x = roi_corners{ROI_n,1}(:,1);
        y = roi_corners{ROI_n,1}(:,2);        
        Img_ROIs_out{ROI_n,1} = Channel_stack{ROI_n,1};
        
        % Get bounding boxes for each labeled region
        temp_label = Dilated_label_ROIs{ROI_n,1};
        temp_label = double(temp_label == ROI_n);
        temp_label(temp_label == 0) = nan; % ignore values not == current ROI#
        min_label = min(temp_label,[],3); % Project the mask to get max contact area
        min_label(isnan(min_label)) = 0; % Make nan values zeros for the 'regionprops function'
        stats = regionprops(min_label, 'BoundingBox','Centroid');

        % Inflate the bounding box to include some area outside the cell
        x_lims = round(stats(1).BoundingBox(1))-pix_pad:(round(stats(1).BoundingBox(1))+round(stats(1).BoundingBox(3)))+pix_pad;
        y_lims = round(stats(1).BoundingBox(2))-pix_pad:(round(stats(1).BoundingBox(2))+round(stats(1).BoundingBox(4)))+pix_pad;
        % Make sure the inflation is not too big
        [x_lims, y_lims] = KLS_stay_within_image_bounds(min_label,x_lims,y_lims);
    
        % Update the masks for the new trucation
        Img_ROIs_out{ROI_n,1} = Img_ROIs_out{ROI_n,1}(y_lims,x_lims,:);
        
        % Update the x values in roi_corners
        roi_corners{ROI_n,1}(1,1) = roi_corners{ROI_n,1}(1,1) + x_lims(1)-1;
        roi_corners{ROI_n,1}(2,1) = roi_corners{ROI_n,1}(2,1) - (size(x(1):x(2),2) - x_lims(end));
        
        % Update the y values in roi_corners
        roi_corners{ROI_n,1}(1,2) = roi_corners{ROI_n,1}(1,2) + y_lims(1)-1;
        roi_corners{ROI_n,1}(2,2) = roi_corners{ROI_n,1}(2,2) - (size(y(1):y(2),2) - y_lims(end)); 
        ROI_n = ROI_n+1;
    end
    clear x y i
end
