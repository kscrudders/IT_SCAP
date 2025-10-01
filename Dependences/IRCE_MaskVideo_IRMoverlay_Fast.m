function [no_output_gen_video] = IRCE_FinalMaskVideo_IRMoverlay_Fast(Save_individual_acq_dir, Specific_ROIs, Dilated_label_ROIs, filtered_IRM_data, IRM_LUT, roi_corners)

    Channels = 1;
    channel_LUTs{1} = IRM_LUT;
    channel_colors{1} = '#ffffff';
    num_ch_included = 1;
    Ch_to_include = 1;
    channel_labels{1} = 'IRM';

    if isempty(Specific_ROIs)
        Specific_ROIs = 1:size(Dilated_label_ROIs,1);
    end
    
    %---------------------------------------------------------%
    % Generate Label Video to Review
    %---------------------------------------------------------%

    folderName = 'Label_Videos_PostManual_Edits'; % Generate a folder name for this data set
    if ~exist(fullfile(Save_individual_acq_dir, folderName), 'dir')
        % If folder doesn't exist, create it
        mkdir(fullfile(Save_individual_acq_dir, folderName));
    end


    % #001 --- Generate label video
    n = 1;
    while  n <= width(Specific_ROIs)
        i = Specific_ROIs(n);
        disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])

        ch1_label = Dilated_label_ROIs{i,1};
        ch1_label = imresize(ch1_label, 2, 'bilinear');   % smooth up-sample

        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);
        ch1_IRM = filtered_IRM_data(min(y):max(y),min(x):max(x),:);
        ch1_IRM = imresize(ch1_IRM, 2, 'bilinear');   % smooth up-sample

        maxT = size(ch1_IRM,3);
        
        max_LUT = max(ch1_label,[],'all');

        if max_LUT == 0
           n = n+1;
           continue;
        end
  
        % Set video object
        video_name = ['Label_Videos_PostManual_Edits_ROI_' num2str(i,'%03.f')];
        vidObj = VideoWriter(fullfile(Save_individual_acq_dir, folderName, [video_name '.mp4']), 'MPEG-4');
        vidObj.FrameRate = round(min(ceil(size(ch1_label,3) / 10), 120)); % dynamically set the frame rate to have a video length of 10s maximum
        vidObj.Quality = 100; % Optional: Adjust video quality if needed 

    
        open(vidObj);

        resized_data{1,1} = ch1_IRM;

        %----------------------------------------------------------
        % Preprocessing: Sharpen and Normalize outside loop
        %----------------------------------------------------------
        for c = 1:Channels
            % Apply sharpening to entire stack
            for z = 1:size(resized_data{c,1},3)
                resized_data{c,1}(:,:,z) = imsharpen(resized_data{c,1}(:,:,z), 'Radius',2, 'Amount',0.8);
            end

            % Normalize data using LUTs
            currData = resized_data{c,1};
            minVal = channel_LUTs{c}(1);
            rangeVal = channel_LUTs{c}(2) - minVal;
            resized_data{c,1} = (currData - minVal) / rangeVal;
            resized_data{c,1}(resized_data{c,1}<0) = 0; 
            resized_data{c,1}(resized_data{c,1}>1) = 1;
        end
    
        %----------------------------------------------------------
        % Combine channels into one 3D array for convenience
        % data dimensions: height x (width*num_ch) x time
        %----------------------------------------------------------
        baseHeight = size(resized_data{1,1},1); % row = height?
        baseWidth  = size(resized_data{1,1},2); % cols = width?
        data = zeros(baseHeight, baseWidth*num_ch_included, maxT, 'like', resized_data{1,1});
        for i = 1:num_ch_included
            chIdx = Ch_to_include(i);
            yStart = (i-1)*baseWidth+1;
            yEnd   = i*baseWidth;
            data(:, yStart:yEnd, :) = resized_data{chIdx,1};
        end
    
        %----------------------------------------------------------
        % Precompute annotations
        %----------------------------------------------------------
        % Scale bar parameters (example values)
        scaleBarLength_um = 10;
        pixelSize_um = 0.157; 
        scaleBarLength_px = round(scaleBarLength_um / pixelSize_um);
        scaleBarY = baseHeight - 20;
        scaleBarX = (baseWidth * Channels) - 20 - (scaleBarLength_um ./ pixelSize_um);
        
        channelFontSize = 12;
    
        % Vertical lines between channels
        channelDivisions = (1:num_ch_included-1)*baseWidth + 1;
    
        % Precompute label positions and colors
        channelTextPositions = zeros(num_ch_included,2);
        for i = 1:num_ch_included
            xPos = 5 + (i-1)*baseWidth; 
            yPos = baseHeight - channelFontSize*1.5; 
            channelTextPositions(i,:) = [xPos,yPos];
        end
    
        channelTextColors = cell(num_ch_included,1);
        for i = 1:num_ch_included
                rgb = sscanf(channel_colors{i}(2:end),'%2x%2x%2x')/255; % convert hex code to rgb
                channelTextColors{i} = rgb';
        end
    
        % Convert data to RGB for annotations
        % Assuming grayscale channels combined side by side, replicate to 3-ch RGB
        % If you have single-channel data, we can just replicate to RGB
        % If channels are already colored, you could map them accordingly.
        % Here we assume grayscale to white. For demonstration:
        % Each pixel is grayscale; convert to RGB by repeating the channel.
        % If you want pseudo-color, you'd need to apply a colormap per channel beforehand.
        RGBframes = zeros(baseHeight, baseWidth*num_ch_included, 3, maxT, 'double');
        for t = 1:maxT
            grayFrame = data(:,:,t);
            RGBframes(:,:,:,t) = repmat(grayFrame,1,1,3);
        end
    
        %----------------------------------------------------------
        % Main loop: Add annotations directly to image
        %----------------------------------------------------------
        for t = 1:maxT
            frameRGB = RGBframes(:,:,:,t);
            maskFrame = ch1_label(:,:,t) > 0; % logical mask for this time-point

            % -------- overlay mask outlines with red/blue logic --------------  
            [B,~]      = bwboundaries(maskFrame,'noholes');
            % Sort boundaries into red (touch edge) vs blue (inside)
            [xR, yR, xB, yB] = LF_prep_boundary_color(B, size(maskFrame,1), size(maskFrame,2));

            % Convert the NaN-separated vectors returned by LF_prep_boundary_color
            % into Nx(2·k) matrices that `insertShape` expects.
            redPolys  = polyVec2Mat(xR,yR);       % cell-safe helper below
            bluePolys = polyVec2Mat(xB,yB);
    
            % ------------ Insert scale bar ----------
            % Draw a white rectangle for scale bar
            frameRGB = insertShape(frameRGB, 'FilledRectangle', ...
                [scaleBarX scaleBarY scaleBarLength_px 4], 'Color', 'white', ...
                'Opacity', 1,'SmoothEdges',true);
            % ------------------------------------------------------------------------
    
            % ------------ Insert channel labels ----------
            for i = 1:num_ch_included
                chText = channel_labels{Ch_to_include(i)};
                frameRGB = insertText(frameRGB, channelTextPositions(i,:), chText, ...
                    'FontSize', channelFontSize, ...
                    'TextColor', channelTextColors{i}./255, 'BoxOpacity', 0);
            end
            % ------------------------------------------------------------------------
    
            % ------------ Insert vertical division lines ----------
            for divPos = channelDivisions
                frameRGB = insertShape(frameRGB, 'Line', [divPos 0 divPos baseHeight], ...
                    'Color', 'white', 'LineWidth', 2,'SmoothEdges',true);
            end
            % ------------------------------------------------------------------------
    
            % ------------ overlay mask outline (1–2 px wide, anti-aliased) ----------
            if ~isempty(redPolys)
                for nn = 1:size(bluePolys,1)
                    poly_in = redPolys(nn,:);
                    poly_in = poly_in(~isnan(poly_in));
                    frameRGB = insertShape(frameRGB,'Polygon',poly_in, ...
                                 'Color','red','LineWidth',2, ...
                                 'SmoothEdges',true);
                end

            end
            if ~isempty(bluePolys)
                for nn = 1:size(bluePolys,1)
                    poly_in = bluePolys(nn,:);
                    poly_in = poly_in(~isnan(poly_in));
                    frameRGB = insertShape(frameRGB,'Polygon',poly_in, ...
                                 'Color','blue','LineWidth',2, ...
                                 'SmoothEdges',true);
                end
            end
            % ------------------------------------------------------------------------

            % Some times the insertShape calls adjust the RGB values out of
            % the range [0 1], check for this and correct the image values.
            if max(frameRGB,[],'all') > 1 || min(frameRGB,[],'all') < 0
                for z = 1:3
                    over_idx = frameRGB(:,:,z,1) > 1;
                    temp = frameRGB(:,:,z,1);
                    temp(over_idx) = 1;
                    frameRGB(:,:,z,1) = temp;

                    under_idx = frameRGB(:,:,z,1) < 0;
                    temp = frameRGB(:,:,z,1);
                    temp(under_idx) = 0;
                    frameRGB(:,:,z,1) = temp;
                end
                
            end

            % Write to video
            writeVideo(vidObj, frameRGB);
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

function polys = polyVec2Mat(xVec,yVec)
% Converts NaN-separated x/y vectors into an M×(2k) matrix for insertShape
    % Empty input → nothing to draw
    if isempty(xVec)
        polys = [];
        return
    end

    % -- split NaN-delimited chains into cell array of polygons ------------
    nanBreak = isnan(xVec);
    idxBreak = [0; find(nanBreak); numel(xVec)+1];
    polyCells = {};
    for k = 1:numel(idxBreak)-1
        seg = (idxBreak(k)+1):(idxBreak(k+1)-1);
        if isempty(seg)     % skip empty fragments
            continue
        end
        xy = [xVec(seg) yVec(seg)];      % [x y] pairs
        polyCells{end+1} = reshape(xy.',1,[]);   %#ok<AGROW>
    end

    if isempty(polyCells)
        polys = [];         % only NaNs were present
        return
    end

    % -- pad rows so that they all have the same (even) length ------------
    rowLen  = cellfun(@numel,polyCells);
    maxLen  = max(rowLen);
    if mod(maxLen,2) ~= 0           % keep an even number of columns
        maxLen = maxLen + 1;
    end

    polys = NaN(numel(polyCells), maxLen,'like',xVec);   % NaN padding
    for k = 1:numel(polyCells)
        polys(k,1:rowLen(k)) = polyCells{k};
    end
end
