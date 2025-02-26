function tracks = KLS_Tracks_in_Mask(binary_img, tracks)
    tracks(tracks == 0) = nan; % remove zero, zero coords
    tracks = KLS_Fill_TrackZeros(tracks); % fill track gaps    

    %---------------------------------------------------------%
    % Filter Tracks for entirely within in the cell ROI
    %---------------------------------------------------------%
    out = ones([size(tracks,1) 1]);
    for i = 1:size(tracks,2) % loop aross time
        x = tracks(:,i,1); % all localization in this frame
        y = tracks(:,i,2);
        x(x == 0) = nan;
        y(y == 0) = nan;

        in_this_frame = x >0; % which tracks have localizations this frame?

        %   boundary on top of the plot
        [row,col] = find(binary_img(:,:,i),1,'first'); % cell mask in this frame
        if ~isempty(row) % if there are localization in this frame
            boundary = bwtraceboundary(binary_img(:,:,i),[row, col],'E');
            x_outline = boundary(:,2);
            y_outline = boundary(:,1);

            [in,~] = inpolygon(x,y, x_outline, y_outline);
            %tracks(:,i,1) = x.*in; % save the localizations in this frame
            %tracks(:,i,2) = y.*in;
            out = out .* ~bitand(in_this_frame,~in); % in the frame, but not in the mask, set bad track = true
        else
            in = zeros(size(x)); % no localizations in the cell mask this frame
            out = out .* ~bitand(in_this_frame,~in); % in the frame, but not in the mask, set bad track = true
        end

        for i = 1:size(tracks,2) % loop aross time
            tracks(:,i,1) = tracks(:,i,1) .* out;
            tracks(:,i,2) = tracks(:,i,2) .* out;
        end
        tracks(tracks == 0) = nan;
        tracks = tracks+(0.5);
    end
end