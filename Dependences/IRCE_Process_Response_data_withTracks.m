function Stats_ROIs = IRCE_Process_Response_data_withTracks(Base_label_ROIs, Ch3_Response_BaSiC_ROIs, sigma_threshold, Stats_ROIs, Response_Tracking_ROIs)
% 

% Edit logs
    % 20240919 KLS, instead of having a per frame threshold for signal by 
    % fitting a guassian over the intensity histogram, then setting the
    % threshold as some multiple of the fit's sigma above the guassian peak
    % [mean] (assuming the majority of pixels are bkgd, not signal)
    % take the median of all those values. This finds the bkgd
    % threshold across all time for a given cell.
    
    i = 1;
    while i <= size(Base_label_ROIs,1)
        ROI_n = i;
        label_base = Base_label_ROIs{i,1};
        Response = Ch3_Response_BaSiC_ROIs{i,1};

        tracks = Response_Tracking_ROIs{ROI_n}.STLN_Tracks;
        tracks(tracks ~= 0) = tracks(tracks ~= 0) + 1; % Trackmate origin is 0,0, Matlab origin is 1,1
        tracks(tracks == 0) = nan; % remove zero, zero coords
        tracks = KLS_Fill_TrackZeros(tracks); % fill track gaps  
        track_is_in = ones([size(tracks,1) 1]); % Tracks that are outside ROI

        %----------------------------------------------------------
        % Check if Ch1, or Ch3 needs to be resized to the largest data in t
        %----------------------------------------------------------
        maxT = max([size(label_base,3), size(Response,3)]);
        label_base = KLS_resizeMatrix(label_base, maxT);
        Response = KLS_resizeMatrix(Response, maxT);
        
        cell_idx = (label_base == ROI_n);
        
        masked_response = Response .* cell_idx;

        guass_thres = zeros([size(Response,3) 1]);
        Stats_ROIs{i,1}.integrated_response_above_background_guass_withlocalization = zeros([size(Response,3) 1]);

        for t = 1:size(Response,3)
            curr_frame_response = masked_response(:,:,t);

            y = curr_frame_response(curr_frame_response ~= 0); % non-zero values

            if ~isnan(y) & length(y) > 3
                binWidth = 0.25;
                % Old way of getting histogram data
                %{
                y2 = curr_frame_response(curr_frame_response > 0);
                h = histogram(y2,'BinWidth',binWidth,'Normalization','PDF');
                x2 = h.BinEdges;
                x2 = movmean(x2,2);
                x2 = x2(2:end);
                y2 = h.Values;
                %}
                
                % Get the data without making a histogram plot
                % Should be faster
                edges = min(y):binWidth:max(y)+binWidth; % Create bin edges
                [y, edges] = histcounts(y, edges, 'Normalization', 'pdf');
                x = edges(1:end-1) + diff(edges)/2;
                
                [xData, yData] = prepareCurveData( x, y );

                % Set up fittype and options.
                ft = fittype( 'gauss1' );
                opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                opts.Display = 'Off';
                opts.Lower = [0 min(x,[],'all') 0];
                opts.Upper = [inf max(x,[],'all') inf];

                max_y = max(y,[],'all');
                idx_max_y = y == max_y;
                amplitude_start = max_y;
                peak_centroid_in_x = x(find(idx_max_y,1));
                gaussian_width = 250;
                %xlim([600 3500])
                opts.StartPoint = [amplitude_start peak_centroid_in_x gaussian_width];

                [fitresult, ~] = fit( xData, yData, ft, opts );

                guass_thres(t) = fitresult.b1+(fitresult.c1*sigma_threshold);
            end
        end

        % median threshold for signal aross all frames
        Response_Intensity_threshold = median(guass_thres(guass_thres > 0));

        for t = 1:size(Response,3)
            curr_frame_response = masked_response(:,:,t);
            
            curr_frame_signal_mask = bitand((curr_frame_response >= Response_Intensity_threshold), cell_idx(:,:,t));
    
            %---------------------------------------------------------%
            % Figure out if any localizations exist in the Response Mask
            % 
            %---------------------------------------------------------%
            x = tracks(:,t,1); % all localization in this frame
            y = tracks(:,t,2);
            x(x == 0) = nan;
            y(y == 0) = nan;
    
            in_this_frame = x > 0; % which tracks have localizations this frame?
            track_is_in = zeros(size(in_this_frame));

            % Works with multi ROIs in one Image
            [B, L] = bwboundaries(curr_frame_signal_mask,'noholes');
            hold on
            valid_boundaries = zeros([1 length(B)]);
           for k = 1:length(B)
               boundary = B{k};
               x_outline = boundary(:,2);
               y_outline = boundary(:,1);

               % -- not needed. The tracked image already had the
               % mask applied to it.
               %[in,~] = inpolygon(x,y, x_outline, y_outline);
               
               in = ~isnan(tracks(:,t,1));
               track_is_in = bitand(in_this_frame, in); % in the frame, and in the mask, set good track = true
               if sum(in) > 0
                   valid_boundaries(k) = 1;
               end
            end   

            if sum(track_is_in) > 0
                valid_boundary_idx = find(valid_boundaries);
                for k = 1:length(valid_boundary_idx)
                    ii = valid_boundary_idx(k);
                    curr_boundary_signal = sum(curr_frame_response(L == ii), 'all');
                    Stats_ROIs{i,1}.integrated_response_above_background_guass_withlocalization(t) = Stats_ROIs{i,1}.integrated_response_above_background_guass_withlocalization(t) + curr_boundary_signal;
                end
            else
                Stats_ROIs{i,1}.integrated_response_above_background_guass_withlocalization(t) = 0;
            end
        end

        i = i +1;
    end
end