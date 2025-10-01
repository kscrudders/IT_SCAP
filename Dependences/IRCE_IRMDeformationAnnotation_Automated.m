function [roi_IRM_interface] = IRCE_IRMDeformationAnnotation_Automated(Save_individual_acq_dir, roi_corners, Ch1_corr_IRM, Base_label_ROIs, IRM_LUT, Dilated_label_ROIs,IRM_thres)

    %---------------------------------------------------------%
    % loop over each cell ROI
    %---------------------------------------------------------%
    for i = 1:size(roi_corners,1) % Parallel loop over manually selected ROIs
        clc
        disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])
        pause(0.5)
        %Grab current ROI pixels: x and y
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);
        img_stack = Ch1_corr_IRM(min(y):max(y),min(x):max(x),:);

        mask_stack = Base_label_ROIs{i,1}  == i;

        cell_IRM_px_intensities = img_stack(mask_stack);

        %-----%
        % Find pixels that are in the lower
        % percentile of IRM destructive interference
        %-----%
        [~, binEdges] = histcounts(cell_IRM_px_intensities, 'BinMethod', 'fd'); % The Freedman-Diaconis rule is less sensitive to outliers in the data, and might be more suitable for data with heavy-tailed distributions. It uses a bin width of 2*iqr(X(:))*numel(X)^(-1/3).

        [counts, edges] = histcounts(cell_IRM_px_intensities, binEdges, ...
                                     'Normalization','pdf');

        edges = movmean(edges,2); edges = edges(2:end);

        sig = [4 3]; % standard deviations from the guassian mean
        [fitresult, lower_thres, upper_thres] = LF_RICMFit_1guass(edges, counts, sig);

        %---------------------------------------------------------%
        % Loop over every frame
        %---------------------------------------------------------%
        close all    
        figure()

        ii = 1;
        while ii <= size(Ch1_corr_IRM,3) % Loop over time
            % skip the current frame if the current cell ROI is absent
            if isempty(find(Base_label_ROIs{i,1}(:,:,ii) == i,1))
                ii = ii+1;
                continue
            end

            img = img_stack(:,:,ii);

            %-----%
            % Mask for the cell and the background
            %-----%
            mask_cell = (Base_label_ROIs{i,1}(:,:,ii) == i);
            SE = strel("disk",25);

            %-----%
            % erode the cell mask to remove the edge of the cell
            %-----%
            mask_cell_erode = imerode(mask_cell,SE);
            mask_bkgd = (Base_label_ROIs{i,1}(:,:,ii) == 0);

        

            % Works with multi ROIs in one Image
            % Add all cell boundary to ch2 and ch3 (impulse and response) 
            [B,~] = bwboundaries(Dilated_label_ROIs{i,1}(:,:,ii) == i,'noholes');
            hold on
            for k = 1:length(B)
               boundary = B{k};
               plot(boundary(:,2), boundary(:,1),'Color',"#00FFFF",'LineStyle',':','LineWidth',2)
            end
            hold off


            
            ii = ii+1;
        end
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

%{
    %-----%
    % Find the image gradients
    %-----%
    [Gmag_cell,~] = imgradient(imgaussfilt(img,2) .* mask_cell);
    [Gmag_bkgd,~] = imgradient(img .* mask_bkgd);

    Gmag_cell = Gmag_cell   .* mask_cell_erode;

    %-----%
    % Laplacian of gaussian filter
    %-----%
    % assume imGrad is your gradient image
    sigma    = 3;                              % approximate spot radius in pixels
    hsize    = 3*ceil(2*sigma)+1;              % filter size
    hLoG     = fspecial('log', hsize, sigma);  % Laplacian‐of‐Gaussian
    resp     = imfilter(Gmag_cell, hLoG, 'same', 'replicate');
    
    %-----%
    % Invert image to flip dark spots to bright, IRM -> we want the
    % dark spots
    %-----%
    resp = imcomplement(resp);
    
    % find local maxima of the LoG response
    peaks    = imregionalmax(resp);
    
    % threshold to discard weak responses
    th       = mean(resp(:)) + 1.0*std(resp(:));  
    spotsBW  = peaks & (resp > th);
    
    % extract centroids
    stats    = regionprops(spotsBW, 'Centroid');
    centroids = vertcat(stats.Centroid);


    %LF_imshow(img, 1, [min(img,[],'all') max(img,[],'all')]); hold on;
    %LF_imshow(resp, 1, [min(resp,[],'all') max(resp,[],'all')]); hold on;

    %-----%
    % Find spots that are in pixels that are in the lower
    % percentile of IRM destructive interference
    %-----%
    SE = strel("disk",9);
    darkest_pixels_mask = img < lower_thres;
    darkest_pixels_mask = imdilate(darkest_pixels_mask,SE);

    LF_imshow(darkest_pixels_mask, 1, [0 1]); hold on;

    if ~isempty(centroids)
        % Convert centroid coordinates to the nearest pixel indices
        rowIdx = round(centroids(:,2));   % Y → rows
        colIdx = round(centroids(:,1));   % X → columns
        
        % Keep only indices that fall inside image bounds
        valid = rowIdx >= 1 & rowIdx <= size(darkest_pixels_mask,1) & ...
                colIdx >= 1 & colIdx <= size(darkest_pixels_mask,2);
        
        rowIdx = rowIdx(valid);
        colIdx = colIdx(valid);
        
        % Linear indices of the centroid pixels
        linIdx = sub2ind(size(darkest_pixels_mask), rowIdx, colIdx);

        % Test whether any centroid lies inside the mask
        anyInside = any(darkest_pixels_mask(linIdx));
        
        % Optional: list of centroids that are inside
        insideCentroids = centroids(valid & darkest_pixels_mask(linIdx), :);

        % visualize
        
        plot(centroids(:,1)+.5, centroids(:,2)+.5, 'o','MarkerEdgeColor',[0 0 1],'LineWidth',2,'MarkerSize',45)
        if anyInside
            plot(insideCentroids(:,1) + 0.5, insideCentroids(:,2) + 0.5, 'o','MarkerEdgeColor','#82fcfd','LineWidth',2,'MarkerSize',45)
        end
    end
    hold off
%}

function [fitresult, lower_thres, upper_thres] = LF_RICMFit_1guass(X, Y, sig)
    %CREATEFIT(X,Y)
    %  Create a fit.
    %
    %  Data for 'untitled fit 1' fit:
    %      X Input : X
    %      Y Output: Y
    %  Output:
    %      fitresult : a fit object representing the fit.
    %      gof : structure with goodness-of fit info.

    %  Auto-generated by MATLAB on 28-Jul-2022 12:19:06
    [xData, yData] = prepareCurveData( X, Y );

    % Set up fittype and options.
    ft = fittype( 'gauss1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [0 min(X,[],'all') 0];
    opts.Upper = [inf max(X,[],'all') inf];
    
    max_y = max(Y,[],'all');
    idx_max_y = Y == max_y;
    amplitude_start = max_y;
    peak_centroid_in_x = X(idx_max_y);
    gaussian_width = 250;
    
    opts.StartPoint = [amplitude_start peak_centroid_in_x gaussian_width];

    % Fit model to data.
    [fitresult, ~] = fit( xData, yData, ft, opts );

    lower_thres = round(fitresult.b1 - (sig(1)*fitresult.c1));
    upper_thres = round(fitresult.b1 + (sig(2)*fitresult.c1));

    xline(lower_thres)
    xline(upper_thres)
end

