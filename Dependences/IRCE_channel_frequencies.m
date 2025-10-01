function channel_freq = IRCE_channel_frequencies(Raw_data, num_ch)
    % Small helper funcation to find the unique frames in each channel
    channel_freq = cell([num_ch, 1]);

    for i = 1:num_ch
        % Pull out each channel's image data. Expectation is that the
        % channel data are interveved, ie Ch1_frame01, Ch2_frame01, Ch1_frame02, Ch2_frame02
        current_ch = Raw_data(:,:,i:num_ch:end);
        
        T = size(current_ch, 3); % number of frames
        frames_2d = reshape(current_ch, [], T); % Instead of dimenaions HxWxT rehaped to 2D vector (H*W) x T
        
        % Compare adjacent columns (frames)
        diff_frames = any(frames_2d(:,2:end) ~= frames_2d(:,1:end-1), 1); % use any so that the search stops after one value is different
        
        unique_mask = false(T,1); % Allicate some memory for the unique frames mask, start with values false

        unique_mask(1) = true; % First frame always unique, set it to true
        unique_mask(2:end) = diff_frames; % True if frame differs from previous

        not_all_zeros_mask = squeeze(any(frames_2d ~= 0,[1]))'; %
        channel_freq{i} = unique_mask & not_all_zeros_mask;
    end
end