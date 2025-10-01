function Stats_ROIs = IRCE_Process_data(Base_label_ROIs, Ch3_Response_BaSiC_ROIs, sigma_threshold, Stats_ROIs)
%

% Edit logs
    % 20240919 KLS, instead of having a per frame threshold for signal by 
    % fitting a guassian over the intensity histogram, then setting the
    % threshold as some multiple of the fit's sigma above the guassian peak
    % [mean] (assuming the majority of pixels are bkgd, not signal)
    % take the median of all those values. This finds the bkgd
    % threshold across all time for a given cell.
    % 20250517 KLS, updated the intital histogram bin edges to using an
    % automated binning method instead of hardcoding a bin width

    i = 1;
    while i <= size(Base_label_ROIs,1)
        ROI_n = i;
        label_base = Base_label_ROIs{i,1};
        Response = Ch3_Response_BaSiC_ROIs{i,1};

        %----------------------------------------------------------
        % Check if Ch1, or Ch3 needs to be resized to the largest data in t
        %----------------------------------------------------------
        maxT = max([size(label_base,3), size(Response,3)]);
        label_base = KLS_resizeMatrix(label_base, maxT);
        Response = KLS_resizeMatrix(Response, maxT);
        
        cell_idx = (label_base == ROI_n);
        
        %landing_idx = find(sum(cell_idx,[1 2]),1,'first');
        %last_contact_idx = find(sum(cell_idx,[1 2]),1,'last');

        masked_response = Response .* cell_idx;

        Stats_ROIs{i,1}.mask_response_signal = zeros(size(Response));
        Stats_ROIs{i,1}.mask_response_bkgd = zeros(size(Response));
        Stats_ROIs{i,1}.integrated_response_above_background_guass = zeros([size(Response,3) 1]);
        Stats_ROIs{i,1}.integrated_background_response = zeros([size(Response,3) 1]);

        Stats_ROIs{i,1}.Response_Intensity_threshold = nan(1);

        guass_thres = zeros([size(Response,3) 1]);
        for t = 1:size(Response,3)
            curr_frame_response = masked_response(:,:,t);

            y = curr_frame_response(curr_frame_response ~= 0); % non-zero values

            if ~isnan(y) & length(y) > 3
                %binWidth = 0.25;
                
                % Get the data, and bin the intensities for a guassian fit
                %edges = min(y):binWidth:max(y)+binWidth; % Create bin edges

                % The Freedman-Diaconis rule is less sensitive to outliers in the data, and might be more suitable for data with heavy-tailed distributions. It uses a bin width of 2*iqr(X(:))*numel(X)^(-1/3).
                [~, edges] = histcounts(y, 'BinMethod', 'fd'); 

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

               assert(numel(xData) > 3, 'Binning has failed to product at least 3 binds. Data likely has narrow intensity range.')

                [fitresult, ~] = fit( xData, yData, ft, opts );

                guass_thres(t) = fitresult.b1+(fitresult.c1*sigma_threshold);
            end
        end

        % What is the bkgd level of Response over all frames?
        Response_Intensity_threshold = median(guass_thres(guass_thres > 0));
        Stats_ROIs{i,1}.Response_Intensity_threshold = Response_Intensity_threshold;
        %Stats_ROIs{i,1}.Response_Intensity_threshold = guass_thres;

        for t = 1:size(Response,3)
            curr_frame_response = masked_response(:,:,t);
            
            % Reponse using fixed threshold across all frames
            %
            Stats_ROIs{i,1}.integrated_response_above_background_guass(t) = sum(curr_frame_response(bitand((curr_frame_response >= Response_Intensity_threshold),cell_idx(:,:,t))), 'all');
            Stats_ROIs{i,1}.integrated_background_response(t) = sum(curr_frame_response(bitand((curr_frame_response < Response_Intensity_threshold),cell_idx(:,:,t))), 'all');
            Stats_ROIs{i,1}.mask_response_signal(:,:,t) = bitand((curr_frame_response >= Response_Intensity_threshold),cell_idx(:,:,t)); % mask of px above bkgd
            Stats_ROIs{i,1}.mask_response_bkgd(:,:,t) = bitand((curr_frame_response < Response_Intensity_threshold),cell_idx(:,:,t)); % mask of px below bkgd
            %}

            % Reponse using threshold set for each frame
            %{
            Stats_ROIs{i,1}.integrated_response_above_background_guass(t) = sum(curr_frame_response(bitand((curr_frame_response >= guass_thres(t)),cell_idx(:,:,t))), 'all');
            Stats_ROIs{i,1}.integrated_background_response(t) = sum(curr_frame_response(bitand((curr_frame_response < guass_thres(t)),cell_idx(:,:,t))), 'all');
            Stats_ROIs{i,1}.mask_response_signal(:,:,t) = bitand((curr_frame_response >= guass_thres(t)),cell_idx(:,:,t)); % mask of px above bkgd
            Stats_ROIs{i,1}.mask_response_bkgd(:,:,t) = bitand((curr_frame_response < guass_thres(t)),cell_idx(:,:,t)); % mask of px below bkgd
            %}
        end

        i = i +1;
    end
end