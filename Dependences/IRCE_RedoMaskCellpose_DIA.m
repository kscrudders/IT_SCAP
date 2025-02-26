function Base_label_ROIs = IRCE_RedoMaskCellpose_DIA(ROIs_to_Redo_Entirely, roi_corners, DIA_Stack, DIA_LUT, Save_individual_acq_dir)
    
    %---------------------------------------------------------%
    % Setup the mask for separating cell ROIs that are incorrectly connected
    %---------------------------------------------------------%
    cd(Save_individual_acq_dir)
    if exist('Base_label_ROIs.mat', 'file')
        load('Base_label_ROIs.mat', 'Base_label_ROIs')
    else
        Base_label_ROIs = cell([size(roi_corners,1), 1]);
    end

    for n = 1:length(ROIs_to_Redo_Entirely) % Loop over manually selected ROIs
        i = ROIs_to_Redo_Entirely(n);

        disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);
        img = DIA_Stack(min(y):max(y),min(x):max(x),:);

        % Basic Segment of that cell ROI in time
        Base_label = zeros(size(img));

        ii = 1;
        while ii <= size(img,3)
            % Loop through time
            ii = ii+1;
        end

        t = 1;
        [rows, cols] = size(img(:,:,t));

        % Plot the boundaries
        L = Base_label(:,:,t);
        [B,~] = bwboundaries(L > 0, 'noholes');
        [x_all_red, y_all_red, x_all_blue, y_all_blue] = LF_prep_boundary_color(B, rows, cols);
        
        % Display the image
        LF_imshow(img, t, DIA_LUT);
        LF_overlay_colored_boundaries(x_all_red, y_all_red, x_all_blue, y_all_blue)
        
        Base_label_ROIs{i,1} = Base_label;
        Base_label_ROIs{i,1} = KLS_Label_previousFrameOverlap(Base_label_ROIs{i,1});
        % Save the updated masks
        save('Base_label_ROIs.mat', 'Base_label_ROIs', '-v7.3')
    end
    close all
end

function labeledImage = LF_remove_edge_labels(labeledImage)
    % Find labels that touch the edge
    % Get the labels at the edges of the image
    edgeLabels = unique([labeledImage(1,:), labeledImage(end,:), labeledImage(:,1)', labeledImage(:,end)']);
    
    % Remove background label (usually label 0, if present)
    edgeLabels(edgeLabels == 0) = [];
    
    % Remove edge labels from the labeled image
    for k = 1:length(edgeLabels)
        labeledImage(labeledImage == edgeLabels(k)) = 0;
    end
end

function continueDrawing = LF_getKey()
    % Initiate the first ask for mask editing
    disp('Add to mask (j), subtract from mask (k), approve (l), go back (;), approve remaining frames (a): ');
    continueDrawing = getkey(1);
    continueDrawing = lower(char(continueDrawing));
    
    % Loop until the input is valid
    while ~any(strcmp(continueDrawing, {'j', 'k', 'l', ';', 'a'}))
        disp('Invalid input. Please enter ''j'', ''k'', ''l'', '';'', or ''a''.');
        continueDrawing = getkey(1);
        continueDrawing = lower(char(continueDrawing));
    end
end

function continueDrawing = LF_getKey_continue()
    disp('Add to mask (j), subtract from mask (k), approve (l), go back (;), approve remaining frames (a): ');
    continueDrawing = getkey(1);
    continueDrawing = lower(char(continueDrawing));
    
    % Loop until the input is valid
    while ~any(strcmp(continueDrawing, {'j', 'k', 'l', ';', 'a'}))
        disp('Invalid input. Please enter ''j'', ''k'', ''l'', '';'', or ''a''.');
        continueDrawing = getkey(1);
        continueDrawing = lower(char(continueDrawing));
    end
end

function ImgH = LF_imshow(img, ii, IRM_LUT)
    % Get the screen size
    screen_size = get(0, 'ScreenSize'); % Returns [left bottom width height]
    half_screen_width = screen_size(3) / 2; % Use half the screen width
    
    % Get the size of the current image
    [~, image_width] = size(img(:, :, ii));
    
    % Calculate the magnification factor
    magnification = (half_screen_width / image_width) * 100; % Magnification in percentage
    
    % Display the image with the calculated magnification
    ImgH = imshow(img(:, :, ii), IRM_LUT, 'InitialMagnification', magnification, 'Border', 'tight');
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
