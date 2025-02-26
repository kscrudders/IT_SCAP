function size_vector = KLS_label2ROIsize(label_img, pixel_size)
    % INPUTS:
    % 1st input: (label_img) labeled image (no label == 0)
    %
    % OUTPUTS:
    % 1st output: (size_vector) 
        % vector of size (time by ROI_n)

    %---------------------------------------------------------%
    % Fit all pixel intensity to a guassian to find threshold
    %---------------------------------------------------------%

    max_ROI = max(label_img,[],'all');

    t = 1;
    max_t = size(label_img,3);

    size_vector = nan([max_t max_ROI]);
    while t <= max_t
        i = 1;
        curr_img = label_img(:,:,t);
        while i <= max_ROI
            ROI_idx = curr_img == i; % id pixels associated with a give ROI as a binary image
            size_vector(t,i) = sum(ROI_idx,'all'); % sum that binary image denoting the size of the ROI in square pixels
            i = i+1;
        end
        t = t+1;
    end
    zero_idx = size_vector == 0;
    
    size_vector = size_vector.*(pixel_size^2); % convert square pixel to square micron
    size_vector(zero_idx) = nan;
end