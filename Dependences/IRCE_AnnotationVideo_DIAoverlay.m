function IRCE_AnnotationVideo_DIAoverlay(Save_individual_acq_dir, Specific_ROIs, roi_corners, DIA_Stack, DIA_LUT)
    % Generates a video showing annotations overlaid on DIA images for each ROI.

    % Load the annotations
    cd(Save_individual_acq_dir);
    if exist('ROIs_Annotations.mat', 'file')
        load('ROIs_Annotations.mat', 'ROIs_Annotations');
    else
        error('ROIs_Annotations.mat not found in the specified directory.');
    end

    % Create a directory to save videos
    video_folder = fullfile(Save_individual_acq_dir, 'Annotation_Videos');
    if ~exist(video_folder, 'dir')
        mkdir(video_folder);
    end

    % Ensure Specific_ROIs is properly defined
    if isempty(Specific_ROIs)
        Specific_ROIs = 1:size(roi_corners,1);
    else
        Specific_ROIs = Specific_ROIs(~isnan(Specific_ROIs));
    end

    % Loop over ROIs
    for n = 1:length(Specific_ROIs)
        i = Specific_ROIs(n);
        disp(['Processing ROI ' num2str(i,'%03.f')]);

        % Get the annotations for this ROI
        annotations = ROIs_Annotations(i, :);

        % Get the image stack for this ROI
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);
        img_stack = DIA_Stack(min(y):max(y), min(x):max(x), :);

        % Create a video writer
        video_name = fullfile(video_folder, ['Annotation_Video_ROI_' num2str(i,'%03.f') '.mp4']);
        vidObj = VideoWriter(video_name, 'MPEG-4');
        vidObj.FrameRate = round(min(size(DIA_Stack,3) / 10, 120)); % dynamically set the frame rate to have a video length of 10s maximum
        open(vidObj);

        num_frames = size(img_stack, 3);

        % Loop over frames
        for frame_idx = 1:num_frames
            img = img_stack(:,:,frame_idx);

            % Get annotations for this frame
            frame_annotation = annotations{frame_idx};
            if isempty(frame_annotation)
                frame_annotation = struct('Target', [], 'Effectors', []);
            end

            % Display image
            h_fig = figure('Visible', 'off');
            LF_imshow(img, DIA_LUT);

            hold on;

            % Overlay target annotation
            if ~isempty(frame_annotation.Target)
                target = frame_annotation.Target;
                center = target.Center;
                radius = target.Radius;

                if strcmp(target.Status, 'Alive')
                    color = 'g';
                else
                    color = 'r';
                end
                viscircles(center, radius, 'Color', color, 'LineWidth', 2);
                % Add text 'Target' near the circle
                text(center(1), center(2) - radius - 10, 'Target', 'Color', color, ...
                    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
            end

            % Overlay effector annotations
            if ~isempty(frame_annotation.Effectors)
                effectors = frame_annotation.Effectors;
                for eff_idx = 1:length(effectors)
                    effector = effectors(eff_idx);
                    center = effector.Center;
                    radius = effector.Radius;
                    if strcmp(effector.Type, 'Suspension')
                        color = 'r';
                    else
                        color = 'b';
                    end
                    viscircles(center, radius, 'Color', color, 'LineWidth', 2);
                    % Add text 'Effector' near the circle
                    text(center(1), center(2) - radius - 10, 'Effector', 'Color', color, ...
                        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
                end
            end

            hold off;

            % Adjust figure properties
            set(gcf, 'Color', [0 0 0]);
            set(gca, 'Visible', 'off');
            set(gcf, 'NextPlot', 'add');
            set(gcf, 'InvertHardCopy', 'off');
            set(gcf, 'PaperPositionMode', 'auto');

            % Capture frame
            currFrame = getframe(gca);
            writeVideo(vidObj, currFrame);

            clf
        end

        close(vidObj);
    end

    close all;
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