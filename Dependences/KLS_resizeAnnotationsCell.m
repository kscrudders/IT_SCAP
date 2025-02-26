function new_roi_IRM_interface = KLS_resizeAnnotationsCell(roi_IRM_interface, maxT)
    t = size(roi_IRM_interface, 2);
    repeatFactor = ceil(maxT / t); % Calculate how many times each slice should be repeated
    
    new_roi_IRM_interface = cell(size(roi_IRM_interface, 1), t * repeatFactor, size(roi_IRM_interface, 3));  % Pre-allocate cellmask to match the size of Ch3 in the third dimension
    
    for i = 1:t
        for j = 1:size(roi_IRM_interface, 1)
            for k = 1:size(roi_IRM_interface, 3)
                for n = 1:repeatFactor
                    new_roi_IRM_interface{j, (repeatFactor * (i - 1)) + n, k} = roi_IRM_interface{j, i, k};
                end
            end
        end
    end
end