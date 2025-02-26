function [roi_mask_Manual_Crops, Base_label_ROIs] = IRCE_RedoMaskActiveContour(ROIs_to_Redo_Entirely, roi_corners, Ch1_corr_IRM, roi_mask_Manual_Crops, IRM_LUT, Base_label_ROIs)
    for n = 1:length(ROIs_to_Redo_Entirely) % Loop over manually selected ROIs
        i = ROIs_to_Redo_Entirely(n);

        disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);
        img = Ch1_corr_IRM(min(y):max(y),min(x):max(x),:);

        % Basic Segment of that cell ROI in time
        Base_label = zeros(size(img));
        
        figure()
        continueDrawing = 'k';
        
        ii = 1;
        while ii <= size(Ch1_corr_IRM,3) % Loop over time
            %imshow(img(:,:,ii),[],'InitialMagnification',200)
            
            if continueDrawing == 'l'
                % Initiate the first ask for mask seperating line(s)
                disp('Add a ROI? (k = yes, l = no:');
                continueDrawing = getkey(1);
                continueDrawing = lower(char(continueDrawing));

                % Loop until the input is 'k','l'
                while ~any(strcmp(continueDrawing, {'k', 'l'}))
                    disp('Invalid input. Please enter ''k - yes'', or ''l - no''.');

                    disp('Add a ROI? (k = yes, l = no): ');
                    continueDrawing = getkey(1);
                    continueDrawing = lower(char(continueDrawing));
                end            
            end
            
            % Keep generating mask seperating lines, when asked.
            if continueDrawing == 'k'
                preconditioned_image = KLS_nonnegative_mixed_norm_preconditioning(img(:,:,ii));
                constructive_interference = preconditioned_image > 0;
                filled_mask = imfill(constructive_interference, 'holes');
                non_cell = constructive_interference + ~filled_mask;

                BWs = ~non_cell; % Rough cell mask
                BWs = imfill(BWs,'holes');

                P = 75; 
                CC = bwconncomp(BWs, 8); % 8-way connectivity
                S = regionprops(CC, 'Area');
                L = labelmatrix(CC);
                small_obj_removed = ismember(L, find([S.Area] >= P));

                BWs = small_obj_removed>0; 

                %RGB = img(:,:,ii);
                %L = superpixels(img(:,:,ii),500);

                %foreground_annotation = drawfreehand('color','b','Closed',1);
                %foreground = createMask(foreground_annotation);

                %bkgd_annotation = drawfreehand('color','r','Closed',1);
               % bkgd = createMask(bkgd_annotation);

                %BWs = lazysnapping(RGB,L,foreground,bkgd);

                %[~,threshold] = edge(img(:,:,ii),'sobel');
                %fudgeFactor = 1.75;
                %BWs = edge(img(:,:,ii),'sobel',threshold * fudgeFactor);
                %imshow(BWs,[])

                SE = strel('disk', 9); % Flat Structuring Element for Image dilation, disk of size 9
                BWsdil = imdilate(BWs,SE);

                P = 150; 
                CC = bwconncomp(BWsdil, 8); % 8-way connectivity
                S = regionprops(CC, 'Area');
                L = labelmatrix(CC);
                small_obj_removed = ismember(L, find([S.Area] >= P));

                BWsdil = small_obj_removed>0;

                %imshow(BWsdil,[])

                BWdfill = imfill(BWsdil,'holes');

                seD = strel('diamond',1);
                BWfinal = imerode(BWdfill,SE);
                BWfinal = imerode(BWfinal,seD);
                %imshow(BWfinal)

                Base_label(:,:,ii) = BWfinal;
            end
                
            ii = ii+1;
        end
        
        % Populate the Non-edited Label for each cell ROI
        Base_label_ROIs{i,1} = Base_label;

        close all
        figure()
        ii = 1;
        while ii <= size(Ch1_corr_IRM,3) % Loop over time
            % slip the current frame if there are no ROI present
            if isempty(find(Base_label(:,:,ii) > 0,1))
                roi_mask_Manual_Crops{i,ii,1} = zeros(size(Base_label,[1 2]));
                ii = ii+1;
                continue
            end

            imshow(img(:,:,ii),IRM_LUT)

            [rows, cols] = size(img(:,:,ii));
            
            L = Base_label(:,:,ii);
            hold on 
            [B,L] = bwboundaries(L>0,'noholes');
            hold on
            for k = 1:length(B)
                boundary = B{k};
                if any(boundary(:,1) == 1) || any(boundary(:,1) == rows) || any(boundary(:,2) == 1) || any(boundary(:,2) == cols)
                    % Boundary touches the edge, color red
                    plot(boundary(:,2), boundary(:,1), 'Color', 'red', 'LineStyle', ':', 'LineWidth', 2)
                else
                    % Boundary does not touch the edge, color yellow
                    plot(boundary(:,2), boundary(:,1), 'Color', 'yellow', 'LineStyle', ':', 'LineWidth', 2)
                end
            end
            hold off

            g = gcf;
            g.WindowState = 'maximized';

            line_count = 0; % Keep track of how many lines have been drawn
            
            % Initiate the first ask for mask seperating line(s)
            continueDrawing = LF_getKey();

            % Keep generating mask seperating lines, when asked.
            new_boundaries = L>0;
            while continueDrawing == 'k'
                line_count = line_count + 1;
                roi_interface{i,ii,line_count} = drawfreehand('color','w','Closed',0);

                SE = strel('disk', 2); % Flat Structuring Element for Image dilation, disk of size 2
                new_boundaries = bitand(new_boundaries, ~imdilate(createMask(roi_interface{i,ii,line_count},L),SE));
                [B,~] = bwboundaries(new_boundaries,'noholes');
                [x_all_red, y_all_red, x_all_blue, y_all_blue] = LF_prep_boundary_color(B , rows, cols);
                LF_overlay_colored_boundaries(x_all_red, y_all_red, x_all_blue, y_all_blue)

                continueDrawing = LF_getKey_continue();
            end
            
            comb_mask = zeros(size(Base_label,[1 2]));
            iii = 1;
            while iii <= line_count
                comb_mask = bitor(comb_mask, createMask(roi_interface{i,ii,iii}));
                iii = iii+1;
            end
            roi_mask_Manual_Crops{i,ii,1} = comb_mask;
            
            if continueDrawing == ';'
               ii = ii - 2; % Go back one frame 
               if ii < 0
                   ii = 0;
               end
            end
            ii = ii+1;
        end
        
        close all
        save('roi_mask_Manual_Crops.mat','roi_mask_Manual_Crops','-v7.3')
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