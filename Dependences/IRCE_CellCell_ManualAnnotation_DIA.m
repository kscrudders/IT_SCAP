function ROIs_Annotations = IRCE_CellCell_ManualAnnotation_DIA(Specific_ROIs, roi_corners, DIA_Stack, DIA_LUT, Save_individual_acq_dir)
    
    %---------------------------------------------------------%
    % Setup the mask for separating cell ROIs that are incorrectly connected
    %---------------------------------------------------------%
    cd(Save_individual_acq_dir)
    if exist('ROIs_Annotations.mat', 'file')
        load('ROIs_Annotations.mat', 'ROIs_Annotations')
    else
        ROIs_Annotations = cell([size(roi_corners,1), size(DIA_Stack,3)]);
    end

    % Loop through each annotation task in a separate phase
    for annotation_phase = 1:2
        % annotation_phase 1: Target cell, 2: Effector cells
        for n = 1:length(Specific_ROIs) % Loop over manually selected ROIs
            i = Specific_ROIs(n);
            
            % Skip if there are already annotations for this ROI
            if all(~cellfun(@isempty, ROIs_Annotations(i, :)))
                disp(['***** Skipping ROI = ' num2str(i,'%03.f') ' (already annotated) *****'])
                continue;
            end
            
            disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])
            x = roi_corners{i,1}(:,1);
            y = roi_corners{i,1}(:,2);
            img_stack = DIA_Stack(min(y):max(y),min(x):max(x),:);

            % Load previous annotations or initialize variables
            if exist('ROIs_Annotations', 'var') && size(ROIs_Annotations, 1) >= i && size(ROIs_Annotations, 2) >= size(img_stack, 3)
                target_annotations = cell(size(img_stack, 3), 1);
                for frame_idx = 1:size(img_stack, 3)
                    if ~isempty(ROIs_Annotations) && size(ROIs_Annotations, 1) >= i && size(ROIs_Annotations, 2) >= frame_idx
                        target_annotations{frame_idx} = ROIs_Annotations{i, frame_idx}.Target;
                    end
                end
                effector_annotations = cell(size(img_stack, 3), 1);
                for frame_idx = 1:size(img_stack, 3)
                    if ~isempty(ROIs_Annotations) && size(ROIs_Annotations, 1) >= i && size(ROIs_Annotations, 2) >= frame_idx
                        effector_annotations{frame_idx} = ROIs_Annotations{i, frame_idx}.Effectors;
                    end
                end
            else
                target_annotations = cell(size(img_stack, 3), 1);
                effector_annotations = cell(size(img_stack, 3), 1);
            end

            switch annotation_phase
                case 1  % Annotate Target Cell
                    disp('Draw annotation for the target cell (Alive/Dead)');
                    disp('Left click to annotate the target cell as Alive, right click to annotate as Dead.');
                case 2  % Annotate Effector Cells
                    clc
                    disp('Annotate effector cells');
                    disp('Left click to annotate as Suspension effector cell, right click to annotate as Adherent effector cell. Press Enter when done.');
            end
            
            % Process each frame for the specific annotation type
            for ii = 1:size(img_stack, 3)
                img = img_stack(:,:,ii);

                % Skip further annotation once the target cell is marked dead
                if annotation_phase == 1 && ii > 1 && ~isempty(target_annotations{ii - 1}) && strcmp(target_annotations{ii - 1}.Status, 'Dead')
                    break;
                end

                switch annotation_phase
                    case 1  % Annotate Target Cell
                        % Draw or edit the mask annotation for the target cell
                        target_annotations{ii} = annotate_target_cell(img, DIA_LUT);
                    case 2  % Annotate Effector Cells
                        % Only proceed if the target cell is still alive
                        if ~isempty(target_annotations{ii}) && strcmp(target_annotations{ii}.Status, 'Alive')
                            effector_annotations{ii} = annotate_effector_cells(img, DIA_LUT);
                        end
                end

                % Save annotations for the current frame
                ROIs_Annotations{i, ii} = struct('Target', target_annotations{ii}, 'Effectors', effector_annotations{ii});
            end
            % Save the updated annotations to disk
            save('ROIs_Annotations.mat', 'ROIs_Annotations', '-v7.3');
        end
    end
    close all;
end

function target_annotation = annotate_target_cell(img, DIA_LUT)
    LF_imshow(img, DIA_LUT);

    [x_click, y_click, button] = ginput(1);
    radius = (30 / 0.157) / 2; % Calculate radius in pixels (30 microns diameter, 0.157 micron/pixel)
    viscircles([x_click, y_click], radius, 'Color', 'g');
    if button == 1  % Left click
        status = 'Alive';
    elseif button == 3  % Right click
        status = 'Dead';
    else
        status = 'Alive'; % Default to Alive if unexpected input
    end
    target_annotation = struct('Center', [x_click, y_click], 'Radius', radius, 'Status', status);
end

function effector_annotation = annotate_effector_cells(img, DIA_LUT)
    LF_imshow(img, DIA_LUT);
    
    effector_cells = [];
    while true
        [x_click, y_click, button] = ginput(1);
        if isempty(x_click)
            break;
        end
        radius = (10 / 0.157) / 2; % Calculate radius in pixels (30 microns diameter, 0.157 micron/pixel)
        
        if button == 1  % Left click
            effector_type = 'Suspension';
            viscircles([x_click, y_click], radius, 'Color', 'r');
        elseif button == 3  % Right click
            effector_type = 'Adherent';
            viscircles([x_click, y_click], radius, 'Color', 'b');
        else
            effector_type = 'Suspension'; % Default to Suspension if unexpected input
            viscircles([x_click, y_click], radius, 'Color', 'r');
        end
        effector_cells = [effector_cells, struct('Center', [x_click, y_click], 'Radius', radius, 'Type', effector_type)]; %#ok<AGROW>
    end
    effector_annotation = effector_cells;
end

function LF_imshow(img, img_LUT)
    % Get the screen size
    screen_size = get(0, 'ScreenSize'); % Returns [left bottom width height]
    half_screen_width = screen_size(3) / 2; % Use half the screen width
    
    % Get the size of the current image
    [~, image_width] = size(img);
    
    % Calculate the magnification factor
    magnification = (half_screen_width / image_width) * 100; % Magnification in percentage
    
    % Display the image with the calculated magnification
    imshow(img, img_LUT, 'InitialMagnification', magnification, 'Border', 'tight');
end
