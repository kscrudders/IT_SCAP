function [x_lims, y_lims] = KLS_stay_within_image_bounds(img,x_lims,y_lims)
    min_x = 1;
    min_y = 1;
    max_x = size(img,2);
    max_y = size(img,1);

    x_idx_min = x_lims <= min_x;    
    if sum(x_idx_min) > 0 
        k = find(x_idx_min,1,'last');
        x_lims = x_lims(k:end);
    end
    x_idx_max = x_lims >= max_x;    
    if sum(x_idx_max) > 0 
        k = find(x_idx_max,1,'first');
        x_lims = x_lims(1:k);
    end
    
    y_idx_min = y_lims <= min_y;    
    if sum(y_idx_min) > 0 
        k = find(y_idx_min,1,'last');
        y_lims = y_lims(k:end);
    end
    y_idx_max = y_lims >= max_y;    
    if sum(y_idx_max) > 0 
        k = find(y_idx_max,1,'first');
        y_lims = y_lims(1:k);
    end
end