function Base_label_ROIs = IRCE_RemoveSmallCellTraces(Base_label_ROIs, Save_individual_acq_dir)
% Last updated: 20241203 KLS
%    - Now accounts for the possibility that a ROI has no label (likely a
%    mistake in the initial mask process. Corrections are later in the
%    code)

    for i = 1:length(Base_label_ROIs)
        label_img = Base_label_ROIs{i}; % x by y by t matrix
        labels = unique(label_img);
        labels(labels == 0) = []; % Exclude the background label (0)
        
        for lbl = labels'
            area_reaches_threshold = false;
            
            if isempty(lbl)
                area_reaches_threshold = true;
            else
                for t = size(label_img, 3):-1:1
                    % searching backward in time gets the 
    
                    % Extract the current slice
                    current_slice = label_img(:, :, t);
                    
                    % Create a binary mask for the current label
                    label_mask = current_slice == lbl;
                    
                    % Calculate the area (number of pixels)
                    area = nnz(label_mask); % number of non-zero elements
                    
                    % Check if the area reaches the minimum threshold
                    if area >= 200 % 200 pixel is rough 5 µm^2 (0.157 µm/px)
                        area_reaches_threshold = true;
                        break; % No need to check further if the threshold is met
                    end
                end
            end

            % If the area never reaches the threshold, set the label to 0
            if ~area_reaches_threshold
                label_img(label_img == lbl) = 0;
            end
        end
        
        % Update the cell array with the modified labeled image
        Base_label_ROIs{i} = label_img;
    end
    cd(Save_individual_acq_dir)
    save('Base_label_ROIs.mat', 'Base_label_ROIs', '-v7.3')
end
