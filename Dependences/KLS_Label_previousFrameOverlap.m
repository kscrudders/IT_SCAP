function label_out = KLS_Label_previousFrameOverlap(label_in)
    % Make a bw_img
    label_in = label_in > 0;

    % Make an empty label image
    label_out = zeros(size(label_in));

    % Start labeling the bw_img using 2D connectivity
    label_out(:,:,1) = bwlabeln(label_in(:,:,1)>0,4);
    
    % Initialize running max label
    running_max_ROI = max(label_out(:,:,1),[],'all');

    % Loop through the remaining frames
    for i = 2:size(label_in,3)
        % Label the next bw_img frame using 2D connectivity
        label_out(:,:,i) = bwlabeln(label_in(:,:,i)>0,4);

        % Get the number of ROIs in the current and previous frames
        num_current_ROIs = max(label_out(:,:,i),[],'all');
        num_previous_ROIs = max(label_out(:,:,i-1),[],'all');

        % Setup cell arrays to store masks for each ROI
        mask_current_ROI = cell(num_current_ROIs, 1);
        mask_previous_ROI = cell(num_previous_ROIs, 1);

        % Store the masks for current frame ROIs
        for ii = 1:num_current_ROIs
            mask_current_ROI{ii} = label_out(:,:,i) == ii;
        end
        
        % Store the masks for previous frame ROIs
        for jj = 1:num_previous_ROIs
            mask_previous_ROI{jj} = label_out(:,:,i-1) == jj;
        end

        % Loop through all the current ROIs
        for ii = 1:num_current_ROIs
            % Calculate the overlap for the current ROI with all previous ROIs
            overlap = zeros(num_previous_ROIs, 1);
            for jj = 1:num_previous_ROIs
                overlap(jj) = sum(mask_current_ROI{ii} & mask_previous_ROI{jj}, 'all');
            end
            
            % Find the previous ROI with the most overlap
            [max_overlap, max_idx] = max(overlap);
            total_pixels_current = sum(mask_current_ROI{ii}, 'all');
            
            % Check if the maximum overlap is greater than 50%
            if (max_overlap / total_pixels_current) > 0.5
                % Use the label of the previous ROI with the most overlap
                temp_label_img = label_out(:,:,i);
                temp_label_img(mask_current_ROI{ii}) = max_idx;
                label_out(:,:,i) = temp_label_img;
            elseif max_overlap > 0
                % Use the label of the previous ROI with the most overlap, even if it's less than 50%
                temp_label_img = label_out(:,:,i);
                temp_label_img(mask_current_ROI{ii}) = max_idx;
                label_out(:,:,i) = temp_label_img;  % Assign it back after updating
            else
                % If no sufficient overlap, use a new label
                temp_label_img = label_out(:,:,i);
                temp_label_img(mask_current_ROI{ii}) = running_max_ROI + 1;
                label_out(:,:,i) = temp_label_img;
                running_max_ROI = running_max_ROI + 1;
            end
        end
    end
end
