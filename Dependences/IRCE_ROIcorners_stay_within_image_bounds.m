function [x, y] = IRCE_ROIcorners_stay_within_image_bounds(img,x,y)
    min_x = 1;
    min_y = 1;
    max_x = size(img,2); % Remember columns are x
    max_y = size(img,1); % Remember rows are y

    x_idx_min = x <= min_x;
    if ~isempty(find(x_idx_min))
        for i = find(x_idx_min)
            x(i) = min_x;
        end
    end
    x_idx_max = x >= max_x;
    if ~isempty(find(x_idx_max))
        for i = find(x_idx_max)
            x(i) = max_x;
        end
    end
    y_idx_min = y <= min_y;
    if ~isempty(find(y_idx_min))
        for i = find(y_idx_min)
            y(i) = min_y;
        end
    end
    y_idx_max = y >= max_y;
    if ~isempty(find(y_idx_max))
        for i = find(y_idx_max)
            y(i) = max_y;
        end
    end
end