function Stats_ROIs = IRCE_Process_Impulse_data_withTracks(Base_label_ROIs, Impulse_BaSiC_ROIs, sigma_threshold, Stats_ROIs, Impulse_Tracking_ROIs)
    i = 1;
    while i <= size(Base_label_ROIs,1)
        ROI_n = i;
        label_base = Base_label_ROIs{i,1};
        Impulse = Impulse_BaSiC_ROIs{i,1};

        tracks = Impulse_Tracking_ROIs{ROI_n}.STLN_Tracks;
        tracks(tracks ~= 0) = tracks(tracks ~= 0) + 1; % Trackmate origin is 0,0, Matlab origin is 1,1
        tracks(tracks == 0) = nan; % remove zero, zero coords
        tracks = KLS_Fill_TrackZeros(tracks); % fill track gaps  
        track_is_in = ones([size(tracks,1) 1]); % Tracks that are outside ROI

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
        
        masked_response = Impulse .* cell_idx; % mask the response for off the cell (bkgd level of receptors)
        masked_bkgd = Impulse .* off_cell_idx; % mask the response for off the cell (bkgd level of receptors)

        guass_thres = zeros([size(Impulse,3) 1]);
        Stats_ROIs{i,1}.integrated_impulse_above_background_withlocalization = zeros([size(Impulse,3) 1]);

        for t = 1:size(Impulse,3)
            curr_frame_response = masked_response(:,:,t);
            curr_frame_bkgd = masked_bkgd(:,:,t);

            y = curr_frame_bkgd(curr_frame_bkgd ~= 0); % non-zero values

            if ~isnan(y)
                guass_thres(t) = median(y);
            end

            curr_frame_signal_mask = bitand((curr_frame_response >= guass_thres(t)), cell_idx(:,:,t));
    
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

               [in,~] = inpolygon(x,y, x_outline, y_outline);

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
                    Stats_ROIs{i,1}.integrated_impulse_above_background_withlocalization(t) = Stats_ROIs{i,1}.integrated_impulse_above_background_withlocalization(t) + curr_boundary_signal;
                end
            else
                Stats_ROIs{i,1}.integrated_impulse_above_background_withlocalization(t) = 0;
            end
        end

        i = i +1;
    end
end