function Stats_ROIs = IRCE_Process_data_Impulse(Base_label_ROIs, Impulse_BaSiC_ROIs, sigma_threshold, Stats_ROIs)
% 20250517 KLS, updated the intital histogram bin edges to using an
% automated binning method instead of hardcoding a bin width

    i = 1;
    while i <= size(Base_label_ROIs,1)
        ROI_n = i;
        label_base = Base_label_ROIs{i,1};
        Impulse = Impulse_BaSiC_ROIs{i,1};

        %----------------------------------------------------------
        % Check if Ch1, or Ch3 needs to be resized to the largest data in t
        %----------------------------------------------------------
        maxT = max([size(label_base,3), size(Impulse,3)]);
        label_base = KLS_resizeMatrix(label_base, maxT);
        Impulse = KLS_resizeMatrix(Impulse, maxT);
        
        cell_idx = (label_base == ROI_n);

        SE = strel('disk', 9); % Flat Structuring Element for Image dilation, disk of size 6
        off_cell_idx = zeros(size(cell_idx));
        for t = 1:size(Impulse,3)
            off_cell_idx(:,:,t) = imdilate(cell_idx(:,:,t),SE);
            off_cell_idx(:,:,t) = off_cell_idx(:,:,t) .* ~cell_idx(:,:,t);
        end
        
        masked_impulse = Impulse .* cell_idx; % mask the response for off the cell (bkgd level of receptors)
        masked_bkgd = Impulse .* off_cell_idx; % mask the response for off the cell (bkgd level of receptors)

        Stats_ROIs{i,1}.mask_impulse_signal = zeros(size(Impulse));
        Stats_ROIs{i,1}.mask_impulse_bkgd = zeros(size(Impulse));

        guass_thres = zeros([size(Impulse,3) 1]);
        Stats_ROIs{i,1}.integrated_impulse_above_background = zeros([size(Impulse,3) 1]);
        Stats_ROIs{i,1}.integrated_impulse_background = zeros([size(Impulse,3) 1]);

        for t = 1:size(Impulse,3)
            curr_frame_response = masked_impulse(:,:,t);
            curr_frame_bkgd = masked_bkgd(:,:,t);

            y = curr_frame_bkgd(curr_frame_bkgd ~= 0); % non-zero values

            if ~isnan(y) & length(y) > 3
                %guass_thres(t) = median(y);

                %binWidth = 0.25;
                
                % Get the data, and bin the intensities for a guassian fit
                %edges = min(y):binWidth:max(y)+binWidth; % Create bin edges using a simple redefined binning

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

            Stats_ROIs{i,1}.integrated_impulse_above_background(t) = sum(curr_frame_response(bitand((curr_frame_response >= guass_thres(t)),cell_idx(:,:,t))), 'all');
            Stats_ROIs{i,1}.integrated_impulse_background(t) = sum(curr_frame_response(bitand((curr_frame_response < guass_thres(t)),cell_idx(:,:,t))), 'all');
            Stats_ROIs{i,1}.mask_impulse_signal(:,:,t) = bitand((curr_frame_response >= guass_thres(t)),cell_idx(:,:,t)); % mask of px above bkgd
            Stats_ROIs{i,1}.mask_impulse_bkgd(:,:,t) = bitand((curr_frame_response < guass_thres(t)),cell_idx(:,:,t)); % mask of px below bkgd
        end

        i = i +1;
    end
end