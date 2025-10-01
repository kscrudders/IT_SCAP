function IRCE_FinalMaskVideo_IRMoverlay(Save_individual_acq_dir, Specific_ROIs, label_ROIs, filtered_IRM_data, IRM_LUT, roi_corners)
    cd(Save_individual_acq_dir)
    clc

    if isempty(Specific_ROIs)
        Specific_ROIs = 1:size(label_ROIs,1);
    else
        Specific_ROIs = Specific_ROIs(~isnan(Specific_ROIs));
    end

    %---------------------------------------------------------%
    % Generate Label Video to Review
    %---------------------------------------------------------%
    cd(Save_individual_acq_dir)
    folderName = 'Label_Videos_PostManual_Edits'; % Generate a folder name for this data set
    if ~exist(fullfile(cd, folderName), 'dir')
        % If folder doesn't exist, create it
        mkdir(folderName);
    end
    % Move to the folder
    cd(folderName);

    % #001 --- Generate label video
    n = 1;
    while  n <= width(Specific_ROIs)
        i = Specific_ROIs(n);

        % If the only values in the mask are
        current_valid_labels = [0 i 999];
        current_labels = unique(label_ROIs{i,1}(:));
        z = length(current_labels);
        switch z
            case 2 % Check if the current ROI has only labels 0 and i
                if sum(ismember(current_labels, current_valid_labels(1:2))) == 2
                    n = i+1;
                    continue
                end
            case 3 % Check if the current ROI has only labels 0, i and 999
                if sum(ismember(current_labels, current_valid_labels)) == 3
                    n = i+1;
                    continue
                end
        end

        disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])

        ch1_label = label_ROIs{i,1};

        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);
        ch1_IRM = filtered_IRM_data(min(y):max(y),min(x):max(x),:);

        max_LUT = max(ch1_label,[],'all');

        if max_LUT == 0
           n = n+1;
           continue;
        end

        % Generate video object
        video_name = ['Label_Video_PostManual_Edits_ROI_' num2str(i,'%03.f')];
        vidObj = VideoWriter([video_name '.mp4'], 'MPEG-4');
        vidObj.FrameRate = round(min(ceil(size(ch1_label,3) / 10), 120)); % dynamically set the frame rate to have a video length of 10s maximum
        vidObj.Quality = 100; % Optional: Adjust video quality if needed
        
        close(vidObj);
        open(vidObj);

        numframes = size(ch1_IRM,3);

        if exist('h','var')
            close(h)
        end
        for t = 1:numframes
            L = ch1_label(:,:,t);
            img = ch1_IRM(:,:,t);

            % convert label image to rgb
            %RGB2 = label2rgb(L, prism(max_LUT),'k');

            % show that color coded labeled image
            imshow(img, IRM_LUT, 'Border', 'tight','InitialMagnification',200)

            hold on 
            s = regionprops(L, 'centroid');
            for ii = 1:numel(s)
                c = s(ii).Centroid; % Get the centroid coordinates
                
                % 20240506 Attempt to speed up labeling when there is a
                % label with a large number, i.e. 999
                    % 20240519 KLS - this works. There is propably a better
                    % method ...
                if isnan(c(1))
                    continue
                end
                text(c(1), c(2), sprintf('%d', ii), 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', 'Color', 'white','FontSize',14);
            end
            
            [B,~] = bwboundaries(L,'noholes');
            [rows, cols] = size(L);
            [x_all_red, y_all_red, x_all_blue, y_all_blue] = LF_prep_boundary_color(B, rows, cols);
            LF_overlay_colored_boundaries(x_all_red, y_all_red, x_all_blue, y_all_blue);

            % Various figure details:
            %screen_size = get(0,'ScreenSize');
            set(gcf,'Color',[0 0 0])
            %set(gcf, 'Position', [100 200 1150 screen_size(4)*0.28]);
            set(gca,'Visible','off')
            set(gcf,'NextPlot','add');
            set(gcf,'InvertHardCopy','off');
            set(gcf,'PaperPositionMode','auto')
            h = gcf;
            currFrame = getframe(h);

            writeVideo(vidObj, currFrame);
            hold off
            %close(h)

            clf
        end

        close(vidObj);

        n = n+1;
    end

    close all
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
            plot(x_all_red, y_all_red, 'Color', 'red', 'LineStyle', ':', 'LineWidth', 3);
        end
        
        % Plot all blue boundaries in one call
        if ~isempty(x_all_blue)
            plot(x_all_blue, y_all_blue, 'Color', 'blue', 'LineStyle', '-', 'LineWidth', 3);
        end
    hold off
end