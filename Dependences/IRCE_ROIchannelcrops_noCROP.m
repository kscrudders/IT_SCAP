function Img_ROIs_out = IRCE_ROIchannelcrops_noCROP(Channel_stack, roi_corners)
    Img_ROIs_out = cell([size(roi_corners,1) 1]);    
    
    for i = 1:size(roi_corners,1)
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);

        Img_ROIs_out{i,1} = Channel_stack(min(y):max(y),min(x):max(x),:);
    end
    
    clear x y i
end
