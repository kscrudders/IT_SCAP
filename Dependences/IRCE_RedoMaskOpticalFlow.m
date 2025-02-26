function [roi_mask_Manual_Crops, Base_label_ROIs] = IRCE_RedoMaskOpticalFlow(ROIs_to_Redo_Entirely, roi_corners, Ch1_corr_IRM, roi_mask_Manual_Crops, IRM_LUT, Base_label_ROIs)
    for n = 1:length(ROIs_to_Redo_Entirely) % Loop over manually selected ROIs
        i = ROIs_to_Redo_Entirely(n);

        disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);
        img = Ch1_corr_IRM(min(y):max(y),min(x):max(x),:);

        % Basic Segment of that cell ROI in time
        cellMask = zeros(size(img));
        
        Th = 0.5;          % Threshold for optical flow magnitude
        diffThreshold = 1;
        sizeFilter = 500;    % Minimum object size
        smoothingDisk = 5;  % Radius for morphological smoothing

        inverted_stack = imcomplement(img);

        cellMask = KLS_segmentCellOpticalFlow(inverted_stack, Th, diffThreshold, sizeFilter, smoothingDisk);

        frameNumber = 10;
        subplot(1,2,1);
        imshow(img(:,:,frameNumber), []);
        title(['Original Frame ', num2str(frameNumber)]);
        
        subplot(1,2,2);
        imshow(cellMask(:,:,frameNumber));
        title(['Segmentation Mask Frame ', num2str(frameNumber)]);
        
        
        %roi_mask_Manual_Crops{i,ii,1} = comb_mask;
            
         
        %save('roi_mask_Manual_Crops.mat','roi_mask_Manual_Crops','-v7.3')

    end
end

function continueDrawing = LF_getKey()
    % Initiate the first ask for mask seperating line(s)
    disp('Add line to seperate masks? (k = yes, l = no, ; = go back 1 frames): ');
    continueDrawing = getkey(1);
    continueDrawing = lower(char(continueDrawing));

    % Loop until the input is 'k','l',';'
    while ~any(strcmp(continueDrawing, {'k', 'l', ';'}))
        disp('Invalid input. Please enter ''k - yes'', ''l - no'', '';'', or ''a''.');

        disp('Add another line? (k = yes, l = no, ; = go back 1 frame): ');
        continueDrawing = getkey(1);
        continueDrawing = lower(char(continueDrawing));
    end
end

function continueDrawing = LF_getKey_continue()
    disp('Add another line? (k = yes, l = no, ; = go back 1 frame): ');
    continueDrawing = getkey(1);
    continueDrawing = lower(char(continueDrawing));

    % Loop until the input is 'k','l',';'
    while ~any(strcmp(continueDrawing, {'k', 'l', ';'}))
        disp('Add another line? (k = yes, l = no, ; = go back 1 frame): ');
        continueDrawing = getkey(1);
        continueDrawing = lower(char(continueDrawing));
    end
end

function [x_all_red, y_all_red, x_all_blue, y_all_blue] = LF_prep_boundary_color(B, rows, cols)
    % Initialize arrays to hold all x and y coordinates for red and blue boundaries
    x_all_red = []; y_all_red = []; x_all_blue = []; y_all_blue = [];
    
    % Loop through each boundary
    for k = 1:length(B)
        boundary = B{k};
    
        % Check if the boundary touches the edge of the image
        if any(boundary(:,1) == 1) || any(boundary(:,1) == rows) || ...
           any(boundary(:,2) == 1) || any(boundary(:,2) == cols)
            % Boundary touches the edge, add to red arrays
            x_all_red = [x_all_red; boundary(:,2); NaN];
            y_all_red = [y_all_red; boundary(:,1); NaN];
        else
            % Boundary does not touch the edge, add to blue arrays
            x_all_blue = [x_all_blue; boundary(:,2); NaN];
            y_all_blue = [y_all_blue; boundary(:,1); NaN];
        end
    end
end

function LF_overlay_colored_boundaries(x_all_red, y_all_red, x_all_blue, y_all_blue)
    hold on
        % Plot all red boundaries in one call
        if ~isempty(x_all_red)
            plot(x_all_red, y_all_red, 'Color', 'red', 'LineStyle', ':', 'LineWidth', 2);
        end
        
        % Plot all blue boundaries in one call
        if ~isempty(x_all_blue)
            plot(x_all_blue, y_all_blue, 'Color', 'blue', 'LineStyle', '-', 'LineWidth', 2);
        end
    hold off
end