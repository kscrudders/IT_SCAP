function [Base_label_ROIs] = IRCE_MaskSeperation(roi_corners, Ch1_corr_IRM, IRM_thres, IRM_LUT, Save_individual_acq_dir)
    % Update 2024101 KLS
    % --- Simplified if-then checks and updated mask editing interface
    % Update 20250318 KLS
    % --- Added logic to check if ROI_corners is bigger than the prior
    %     saved data
    % Update 20250510 KLS
    % --- Added an adaptive filter threshold to better segment areas where
    %     scattered light elevates the intensity
    % Update 20250514 KLS
    % --- Options commented line that combines destrcutive interference
    %     with constructive interference. Could automate that threshold
    %     using the same workflow as IRM_thres
    % Update 20250602 KLS
    % --- Added logic to log manual mask annotations such that overlaping
    %     ROIs receive the same annotations, thus saving time when the same
    %     annotations are needed for different ROIs
    %---------------------------------------------------------------------%
    % Manually Edit Masks -- Adjust ROIs with Paint Brush Interface
    %---------------------------------------------------------------------%
    if exist(fullfile(Save_individual_acq_dir,'globalAnnotMask.mat'), 'file')
        load(fullfile(Save_individual_acq_dir,'globalAnnotMask.mat'), 'globalAnnotMask')
    else
        globalAnnotMask = true(size(Ch1_corr_IRM));
    end    

    %---------------------------------------------------------%
    % Adaptive threshold setting
    %---------------------------------------------------------%
    Base_label_ROIs = cell([size(roi_corners,1), 1]);

    gMed = median(Ch1_corr_IRM,'all'); % global median intensity

    BW_adaptive_thres = nan(size(Ch1_corr_IRM)); % Preallocate memory

    winsz     = [129 129];       % sliding window size (should be the 1.5-2.5x the size of a cell (10 µm == 64 px)
    h = ones(winsz)/prod(winsz);       % 129×129 averaging kernel, example [1 1 1; 1 1 1; 1 1 1;] -> [1/9 1/9 1/9]    
    for t = 1:size(Ch1_corr_IRM,3)
        I = Ch1_corr_IRM(:,:,t);

        %––– rolling average –––
        Iavg = imfilter(I, h, 'symmetric');  % or conv2(I, h, 'same')
        
        %––– build a threshold map –––
        %  increase above IRM_thres by (Iavg − gMed), but never drop below IRM_thres
        thrMap = IRM_thres + max( (Iavg - gMed), 0);
        
        %––– 4) binarize –––
        BW_adaptive_thres(:,:,t) = I < thrMap; % cells
        
        % 20250514 KLS
        % If constructive interference is self contained under cells:
        % Hard code threshold
        %BW_adaptive_thres(:,:,t) = I < thrMap | I > 2700; % cells
    end
    
    %---------------------------------------------------------%
    % Setup the mask for separating cell ROIs that are incorrectly connected
    %---------------------------------------------------------%
    %cd(Save_individual_acq_dir) % edited to avoid moving the directory KLS 20250602

    if exist(fullfile(Save_individual_acq_dir,'Base_label_ROIs.mat'), 'file') % edited to avoid moving the directory KLS 20250602
        load(fullfile(Save_individual_acq_dir,'Base_label_ROIs.mat'), 'Base_label_ROIs') % edited to avoid moving the directory KLS 20250602
    end
    
    if exist('Base_label_ROIs','var')
        if size(Base_label_ROIs,1) < size(roi_corners,1)
            need_2_mask_flag = 1;
            i_start = size(roi_corners,1)+1;
        else
            need_2_mask_flag = 0;
        end
    else
        need_2_mask_flag = 1;
        i_start = 1;
    end

    any_empty_flag = 0;
    for i = 1:size(Base_label_ROIs,1)
        if isempty (Base_label_ROIs{i,1})
            any_empty_flag = 1;
            i_start = 1;
        end
    end

    if need_2_mask_flag == 1 | any_empty_flag == 1
        %---------------------------------------------------------%
        % Populate empty ROI masks or add more ROIs
        %---------------------------------------------------------%
        % Loop over manually selected ROIs
        for i = i_start:size(roi_corners,1)
            if ~isempty(Base_label_ROIs{i})
                continue;
            end
            disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])
            
            x = roi_corners{i,1}(:,1);
            y = roi_corners{i,1}(:,2);
            img = Ch1_corr_IRM(min(y):max(y),min(x):max(x),:);
            
            %---------------------------------------------------------%
            % Adaptive Threshold
            %---------------------------------------------------------%
            % Adaptive Threshold Segment of that cell ROI in time
            Base_label = BW_adaptive_thres(min(y):max(y),min(x):max(x),:);
            

            %---------------------------------------------------------%
            % Clean up that basic Segmentation
            %---------------------------------------------------------%
            
            % Remove Small ROIs
            P = 202; % at least 5 µm^2 in size (0.157 µm/px)
            small_obj_removed = zeros(size(img));
            for t = 1:size(img,3)
                CC = bwconncomp(Base_label(:,:,t), 8); % 8-way connectivity
                S = regionprops(CC, 'Area');
                L = labelmatrix(CC);
                small_obj_removed(:,:,t) = ismember(L, find([S.Area] >= P));
            end
            small_obj_removed = small_obj_removed > 0;
            
            SE = strel('disk', 6); % Structuring Element for dilation
            
            % Dilate Mask
            for t = 1:size(img,3)
                Base_label(:,:,t) = imdilate(small_obj_removed(:,:,t), SE);
            end
            
            % Fill holes (gaps) in mask
            for t = 1:size(img,3)
                Base_label(:,:,t) = imfill(Base_label(:,:,t), 'holes');
            end
            
            % Erode mask to remove the effect of dilation
            for t = 1:size(img,3)
                Base_label(:,:,t) = imerode(Base_label(:,:,t), SE);
            end
    
            Base_label = Base_label > 0;
            Base_label = Base_label .* globalAnnotMask(min(y):max(y),min(x):max(x),:); % apply prior annotations by remove label pixel values where globalAnnotMask == 0
            %---------------------------------------------------------%
            % Manually Edit Masks
            %---------------------------------------------------------%
            close all
            figure()
            
            % Initialize variables
            approve_remaining_flag = 0; % If the user wants to approve remaining frame without edits
            initial_valid_mask_flag = 0; % Only start prompting user for edits if there has been one valid mask first
            go_quick_flag = 1; % At first, auto-prompt the user for substractive edits to speed up editing
            go_back_flag = 0;
            ii = 1;
            while ii <= size(Ch1_corr_IRM,3) % Loop over time frames
                
                % Skip the current frame if there are no ROI present
                if isempty(find(Base_label(:,:,ii) > 0, 1))
                    if go_back_flag == 1 % move back in frames logic
                        ii = ii - 2;
                        if ii <= 0
                            ii = 0;
                            go_back_flag = 0;
                        end
                    end
                    ii = ii + 1;
                    continue;
                end
                
                go_back_flag = 0; % reset flag
                
                % Repeat the mask annotation if the mask has not changed
                if ii > 1 && all(Base_label(:,:,ii-1) == Base_label(:,:,ii), 'all')
                    ii = ii + 1;
                    continue;
                end
                
                if approve_remaining_flag == 0 % Skip correcting masks if user wants to auto-approve remaining frames
                    [rows, cols] = size(img(:,:,ii));
                    
                    % Plot the boundaries
                    L = Base_label(:,:,ii);
                    L = L .* globalAnnotMask(min(y):max(y),min(x):max(x),ii); % apply prior annotations by remove label pixel values where globalAnnotMask == 0

                    [B,~] = bwboundaries(L > 0, 'noholes');
                    [x_all_red, y_all_red, x_all_blue, y_all_blue] = LF_prep_boundary_color(B, rows, cols);
                    
                    % If there is 1+ valid ROI, skip prompting for edits
                    if ~isempty(x_all_blue) >= 1
                        initial_valid_mask_flag = 1;
                        ii = ii+1;
                        continue
                    end
    
                    % If there hasn't yet been a valid ROI and there are
                    % only invalid ROIs in the current frame, skip
                    % prompting for edits
                    if initial_valid_mask_flag == 0 && isempty(x_all_blue)
                        ii = ii+1;
                        continue
                    end

                    % Display the image
                    ImgH = LF_imshow(img, ii, IRM_LUT);
                    LF_overlay_colored_boundaries(x_all_red, y_all_red, x_all_blue, y_all_blue)
                    
                    % Initialize new boundaries
                    new_boundaries = L > 0;
                    
                    if go_quick_flag == 1
                        % For the go quick setting, only subtractions are
                        % going to be asked for.
                        % If for some reason subtraction is not needed,
                        % then after 3 subtraction edits where no valid ROI
                        % is generated, then the options for edits will
                        % switch to the full set of options
                        continueDrawing = 'k';
                        
                        SE_line_dilate = strel('disk', 3); % Structuring Element for dilation
                        
                        frustration_counter = 0;
                        while any(strcmp(continueDrawing, {'j', 'k'}))
                            if continueDrawing == 'k'
                                % Subtract from the mask
                                roi_subtract = drawfreehand('Color', 'r','LineWidth',5,'Closed',false,'DrawingArea','auto','Multiclick',true); % Draw area to subtract
                                mask_to_subtract = createMask(roi_subtract, ImgH);
                                mask_to_subtract= imdilate(mask_to_subtract, SE_line_dilate);
        
                                globalAnnotMask(min(y):max(y),min(x):max(x),ii) = globalAnnotMask(min(y):max(y),min(x):max(x),ii) & ~mask_to_subtract; % update global annotation mask
                                new_boundaries = new_boundaries & ~mask_to_subtract;
                                delete(roi_subtract);
                            end
        
                            ImgH = LF_imshow(img, ii, IRM_LUT);
                            % Update the displayed mask boundaries
                            [B,~] = bwboundaries(new_boundaries, 'noholes');
                            [x_all_red, y_all_red, x_all_blue, y_all_blue] = LF_prep_boundary_color(B, rows, cols);
                            LF_overlay_colored_boundaries(x_all_red, y_all_red, x_all_blue, y_all_blue)
                            
                            % Ask the user if they want to continue drawing
                            %continueDrawing = LF_getKey_continue();
        
                            % If there is 1+ valid ROI, 
                            if ~isempty(x_all_blue) >= 1
                                continueDrawing = 'l';
                                continue
                            else
                                continueDrawing = 'k';
                                frustration_counter = frustration_counter+1;
                            end

                            if frustration_counter == 3
                                continueDrawing = ';';
                                go_quick_flag = 0;
                            end
                        end
                        
                        % Handle 'l', ';', 'a' cases
                        if continueDrawing == 'l'
                            % User approves the current mask
                            %Base_label(:,:,ii) = new_boundaries; % This line
                            % occurs after this if statememt
                        elseif continueDrawing == ';'
                            % User wants to go back one frame
                            ii = ii - 1;
                            if ii < 1
                                ii = 1;
                            end
                            go_back_flag = 1;
                            continue;
                        elseif continueDrawing == 'a'
                            % User wants to approve remaining frames
                            approve_remaining_flag = 1;
                            %Base_label(:,:,ii) = new_boundaries; % This line
                            % occurs after this if stateme
                        end
                    else
                        % Get user input
                        continueDrawing = LF_getKey();
                        
                        SE_line_dilate = strel('disk', 3); % Structuring Element for dilation
                        while any(strcmp(continueDrawing, {'j', 'k'}))
                            if continueDrawing == 'j'
                                % Add to the mask
                                roi_add = drawfreehand('Color', 'g','LineWidth',5,'Closed',true,'DrawingArea','auto','Multiclick',false); % Draw area to add
                                mask_to_add = createMask(roi_add, ImgH);
                                mask_to_add= imdilate(mask_to_add, SE_line_dilate);

                                globalAnnotMask(min(y):max(y),min(x):max(x),ii) = globalAnnotMask(min(y):max(y),min(x):max(x),ii) | mask_to_add; % update global annotation mask
                                new_boundaries = new_boundaries | mask_to_add;
                                delete(roi_add);
                            elseif continueDrawing == 'k'
                                % Subtract from the mask
                                roi_subtract = drawfreehand('Color', 'r','LineWidth',5,'Closed',false,'DrawingArea','auto','Multiclick',true); % Draw area to subtract
                                mask_to_subtract = createMask(roi_subtract, ImgH);
                                mask_to_subtract= imdilate(mask_to_subtract, SE_line_dilate);
        
                                globalAnnotMask(min(y):max(y),min(x):max(x),ii) = globalAnnotMask(min(y):max(y),min(x):max(x),ii) & ~mask_to_subtract; % update global annotation mask
                                new_boundaries = new_boundaries & ~mask_to_subtract;
                                delete(roi_subtract);
                            end
        
                            ImgH = LF_imshow(img, ii, IRM_LUT);
                            % Update the displayed mask boundaries
                            [B,~] = bwboundaries(new_boundaries, 'noholes');
                            [x_all_red, y_all_red, x_all_blue, y_all_blue] = LF_prep_boundary_color(B, rows, cols);
                            LF_overlay_colored_boundaries(x_all_red, y_all_red, x_all_blue, y_all_blue)
                            
                            % Ask the user if they want to continue drawing
                            %continueDrawing = LF_getKey_continue();
        
                            % If there is 1+ valid ROI, 
                            if ~isempty(x_all_blue) >= 1
                                continueDrawing = 'l';
                                continue
                            else
                                continueDrawing = LF_getKey_continue();
                            end
                        end
                        
                        % Handle 'l', ';', 'a' cases
                        if continueDrawing == 'l'
                            % User approves the current mask
                            %Base_label(:,:,ii) = new_boundaries; % This line
                            % occurs after this if statememt
                        elseif continueDrawing == ';'
                            % User wants to go back one frame
                            ii = ii - 1;
                            if ii < 1
                                ii = 1;
                            end
                            go_back_flag = 1;
                            continue;
                        elseif continueDrawing == 'a'
                            % User wants to approve remaining frames
                            approve_remaining_flag = 1;
                            %Base_label(:,:,ii) = new_boundaries; % This line
                            % occurs after this if stateme
                        end
                    end
                    % Update the mask
                    Base_label(:,:,ii) = new_boundaries;
                end
                
                ii = ii + 1;
            end
            
            % Update the Base_label_ROIs
            %Base_label_ROIs{i,1} = KLS_Label_previousFrameOverlap(Base_label);
            %Base_label_ROIs{i,1} = LF_remove_edge_labels(Base_label_ROIs{i,1});
    
            Base_label_ROIs{i,1} = KLS_Label_fullStackTracking(Base_label);      
    
            % Save the updated masks
            save(fullfile(Save_individual_acq_dir,'Base_label_ROIs.mat'), 'Base_label_ROIs', '-v7.3')
            % Save the global mask annotation
            save(fullfile(Save_individual_acq_dir,'globalAnnotMask.mat'), 'globalAnnotMask', '-v7.3')
        end
    end
   
    clc
    close all
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
