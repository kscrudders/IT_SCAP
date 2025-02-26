function Img_ROIs_out = IRCE_ROIchannelcrops(Dilated_label_ROIs, Channel_stack, roi_corners)
    Img_ROIs_out = cell([size(Dilated_label_ROIs,1) 1]);    
    
    for i = 1:size(Dilated_label_ROIs,1)
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);

        Img_ROIs_out{i,1} = Channel_stack(min(y):max(y),min(x):max(x),:);
    end
    
    pix_pad = 11;

    max_ROI = size(Dilated_label_ROIs,1);
    ROI_n = 1;
    while ROI_n <= max_ROI
        % Get bounding boxes for each labeled region
        temp_label = Dilated_label_ROIs{ROI_n,1};
        temp_label = double(temp_label == ROI_n);
        temp_label(temp_label == 0) = nan; % ignore values not == current ROI#
        min_label = min(temp_label,[],3); % Project the mask to get max contact area
        min_label(isnan(min_label)) = 0; % Make nan values zeros for the 'regionprops function'

        if unique(min_label) == 0
            error(['ROI #' num2str(ROI_n) ' is an empty mask, redo that ROI mask'])
        end

        stats = regionprops(min_label, 'BoundingBox','Centroid');

        % Inflate the bounding box to include some area outside the cell
        x_lims = round(stats(1).BoundingBox(1))-pix_pad:(round(stats(1).BoundingBox(1))+round(stats(1).BoundingBox(3)))+pix_pad;
        y_lims = round(stats(1).BoundingBox(2))-pix_pad:(round(stats(1).BoundingBox(2))+round(stats(1).BoundingBox(4)))+pix_pad;
        % Make sure the inflation is not too big
        [x_lims, y_lims] = KLS_stay_within_image_bounds(min_label,x_lims,y_lims);
    
        Img_ROIs_out{ROI_n,1} = Img_ROIs_out{ROI_n,1}(y_lims,x_lims,:);
        ROI_n = ROI_n+1;
    end
    
    clear x y i
end
