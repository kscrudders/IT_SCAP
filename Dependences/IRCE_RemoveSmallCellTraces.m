function Base_label_ROIs = IRCE_RemoveSmallCellTraces(Base_label_ROIs, Save_individual_acq_dir)
% Updated: 20241203 KLS
%    - Now accounts for the possibility that a ROI has no label (likely a
%    mistake in the initial mask process. Corrections are later in the
%    code)
% Updated: 20250923 KLS
%    - Added an additional fillholes to correct for importing prior
%    correction, if the data has previously been segemented and hand
%    corrected

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
        
        label_img = KLS_fillLabelHoles(label_img); % Fill any holes on a per slice basis

        % Update the cell array with the modified labeled image
        Base_label_ROIs{i} = label_img;
    end
    cd(Save_individual_acq_dir)
    save('Base_label_ROIs.mat', 'Base_label_ROIs', '-v7.3')
end


function filledStack = KLS_fillLabelHoles(labelStack)
%KLS_FILLLABELHOLES Fill holes inside labeled regions in a x-y-t stack
%
% Input:
%   labelStack - integer label image stack, size HxWxT
%                0 = background, 1..n = labels
%
% Output:
%   filledStack - same size, with interior holes filled for each label

    %------------------- Validate inputs -------------------%
    p = inputParser;
    addRequired(p, 'labelStack', ...
        @(x) isnumeric(x) && ndims(x)==3 && all(x(:) >= 0));
    parse(p, labelStack);

    [H, W, T] = size(labelStack); % input dimensions
    filledStack = zeros(H, W, T, 'like', labelStack); % What ever data type the input was

    %------------------- Process each frame -------------------%
    for t = 1:T
        L = labelStack(:,:,t);
        filledFrame = zeros(H, W, 'like', L);

        % Process each label individually
        labels = setdiff(unique(L), 0); % exclude background (0)
        for k = 1:numel(labels)
            mask = (L == labels(k)); % current label #

            % Fill interior holes
            maskFilled = imfill(mask, 'holes');

            % Write back into frame using the current label #
            filledFrame(maskFilled) = labels(k);
        end

        filledStack(:,:,t) = filledFrame;
    end
end