function Base_label_ROIs = IRCE_RedoMaskCellpose(ROIs_to_Redo_Entirely, roi_corners, Ch1_corr_IRM, cp, Base_label_ROIs, IRM_LUT, Save_individual_acq_dir)
    averageCellDiameter = 75;

    cd(Save_individual_acq_dir)
    for n = 1:length(ROIs_to_Redo_Entirely) % Loop over manually selected ROIs
        i = ROIs_to_Redo_Entirely(n);

        disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);
        img = Ch1_corr_IRM(min(y):max(y),min(x):max(x),:);

        % Basic Segment of that cell ROI in time
        Base_label = zeros(size(img));

        ii = 1;
        while ii <= size(img,3)
            labels = segmentCells2D(cp,img(:,:,ii), ImageCellDiameter=averageCellDiameter);
            labels = LF_remove_edge_labels(labels);
            Base_label(:,:,ii) = labels;
            ii = ii+1;
        end
        
        %---------------------------------------------------------%
        % Manually Edit Masks
        %---------------------------------------------------------%
        close all
        figure()
        
        % Initialize variables
        approve_remaining_flag = 0;
        ii = 1;
        while ii <= size(Ch1_corr_IRM,3) % Loop over time frames
            
            % Skip the current frame if there are no ROI present
            if isempty(find(Base_label(:,:,ii) > 0, 1))
                ii = ii + 1;
                continue;
            end
            
            % Repeat the mask correction if the mask has not changed
            if ii > 1 && all(Base_label(:,:,ii-1) == Base_label(:,:,ii), 'all')
                ii = ii + 1;
                continue;
            end
            
            if approve_remaining_flag == 0 % Skip correcting masks if user wants to auto-approve remaining frames
                [rows, cols] = size(img(:,:,ii));
                
                % Plot the boundaries
                L = Base_label(:,:,ii);
                [B,~] = bwboundaries(L > 0, 'noholes');
                [x_all_red, y_all_red, x_all_blue, y_all_blue] = LF_prep_boundary_color(B, rows, cols);
                
                % Display the image
                ImgH = LF_imshow(img, ii, IRM_LUT);
                LF_overlay_colored_boundaries(x_all_red, y_all_red, x_all_blue, y_all_blue)
                
                % Initialize new boundaries
                new_boundaries = L > 0;
                
                % Get user input
                continueDrawing = LF_getKey();
                
                while any(strcmp(continueDrawing, {'j', 'k'}))
                    if continueDrawing == 'j'
                        % Add to the mask
                        roi_add = drawfreehand('Color', 'g'); % Draw area to add
                        mask_to_add = createMask(roi_add, ImgH);
                        new_boundaries = new_boundaries | mask_to_add;
                        delete(roi_add);
                    elseif continueDrawing == 'k'
                        % Subtract from the mask
                        roi_subtract = drawfreehand('Color', 'r'); % Draw area to subtract
                        mask_to_subtract = createMask(roi_subtract, ImgH);
                        new_boundaries = new_boundaries & ~mask_to_subtract;
                        delete(roi_subtract);
                    end
                    
                    % Update the displayed mask boundaries
                    [B,~] = bwboundaries(new_boundaries, 'noholes');
                    [x_all_red, y_all_red, x_all_blue, y_all_blue] = LF_prep_boundary_color(B, rows, cols);
                    LF_overlay_colored_boundaries(x_all_red, y_all_red, x_all_blue, y_all_blue)
                    
                    % Ask the user if they want to continue drawing
                    continueDrawing = LF_getKey_continue();
                end
                
                % Handle 'l', ';', 'a' cases
                if continueDrawing == 'l'
                    % User approves the current mask
                    Base_label(:,:,ii) = new_boundaries;
                elseif continueDrawing == ';'
                    % User wants to go back one frame
                    ii = ii - 2;
                    if ii < 1
                        ii = 1;
                    end
                    continue;
                elseif continueDrawing == 'a'
                    % User wants to approve remaining frames
                    approve_remaining_flag = 1;
                    Base_label(:,:,ii) = new_boundaries;
                end
                
                % Update the mask
                Base_label(:,:,ii) = new_boundaries;
            end
            
            ii = ii + 1;
        end

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
