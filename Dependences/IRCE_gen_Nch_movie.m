function IRCE_gen_Nch_movie(Time_stamps_address, ROI_n, ...
    Individual_ROI_data, Dilated_label_forROIs, channel_freq, channel_LUTs, ...
    channel_labels, channel_colors, roi_IRM_interface, STLN_Tracks,...
    Impulse_track_annotation_channels, Video_Name, Reorder_channels)
    
    if isempty(STLN_Tracks)
        no_Impulse_tracks_flag = 1;
    else
        no_Impulse_tracks_flag = 0;
    end

    Channels = size(Individual_ROI_data,1);

    %----------------------------------------------------------
    % Check if Ch1, Ch2 or Ch3 needs to be resized to the largest data in t
    %----------------------------------------------------------
    largest_ch_idx = find(channel_freq == min(channel_freq));
    maxT = size(Individual_ROI_data{largest_ch_idx(1),1},3);
    
    resized_data = cell(size(Individual_ROI_data));
    for i = 1:Channels
        resized_data{i,1} = KLS_resizeMatrix(Individual_ROI_data{i,1}, maxT);
    end
    cell_mask = KLS_resizeMatrix(Dilated_label_forROIs, maxT);
    
    Timestamps = IRCE_ND2_TimeStamps(Time_stamps_address);

    % Identify when the cell lands (ie when the masking starts)
    size_vector = KLS_label2ROIsize(cell_mask == ROI_n, 0.157);
    first_frame_landed = find(size_vector > 0,1,'first'); % first frame the cell landed

    % First frame showing a IRM deformation annotation
    deformation_flag = 0; % Pre-populate a deformation flag
    IRM_deformation_frame = [];
    
    for j = 1:size(roi_IRM_interface, 2)
        for k = 1:size(roi_IRM_interface, 3)
            if ~isempty(roi_IRM_interface{1, j, k})
                IRM_deformation_frame = j;
                deformation_flag = 1;
                break;
            end
        end
            % if there is a deformation stop looking for more
            if ~isnan(IRM_deformation_frame)
                break;
            end
    end

    % Generate video object
    vidObj = VideoWriter([Video_Name '.mp4'], 'MPEG-4');
    vidObj.FrameRate = round(min(maxT / 10, 120)); % dynamically set the frame rate to have a video length of 10s maximum
    vidObj.Quality = 100;
    close(vidObj);

    channel_text_font_size = 12;
    %annotation_text_font_size = 14;
    offset_from_bot = 6;
    offset_from_left = 5;

    numframes = size(resized_data{1,1},3);

    % Reorder endpoint to get this information before the movie, then add 
    %mean_eGFP = 1204; % Manual epi GFP expression pull from ImageJ, (590 bgkd e488)

    %----------------------------------------------------------
    % KLS normalize each channel's LUT 0-1
    %----------------------------------------------------------
    normalize_data = cell(size(Individual_ROI_data));
    for i = 1:Channels
        normalize_data{i,1} = resized_data{i,1} - channel_LUTs{i,1}(1); % set minimum
        normalize_data{i,1} = normalize_data{i,1} ./ (channel_LUTs{i,1}(2) - channel_LUTs{i,1}(1)); % set max to 1, account for lower bound shift
    end

    %----------------------------------------------------------
    % Pull out unique timepoints in the impulse and response data
    % so as to correctly link tracks
    %----------------------------------------------------------
    % Extract the relevant data
    if ~isempty(Impulse_track_annotation_channels)
        impulse_tracking_data = normalize_data{Impulse_track_annotation_channels(1), 1};
            % Find the unique time points in the 3rd dimension using your logic
            [x, y, t] = size(impulse_tracking_data);
            reshapedData = reshape(impulse_tracking_data, x*y, t);
            [~, uniqueIndices] = unique(reshapedData', 'rows');
            unique_impulse_timepoints = sort(uniqueIndices);
    end
    
    %----------------------------------------------------------
    % Combine Individual image channel variables
    % into one variable
    %----------------------------------------------------------
    data = zeros([size(normalize_data{1,1},1), size(normalize_data{1,1},2)*Channels, size(normalize_data{1,1},3)]);
    for i = 1:Channels
        x = 1:size(normalize_data{1,1},1);
        y = ((i-1)*(size(normalize_data{1,1},2)))+1:i*(size(normalize_data{1,1},2));
        data(x, y, :) = normalize_data{Reorder_channels(i),1};
    end
    
    % WIP: KLS Find a way to systematically define activation frame
    
    %Timestamps = zero_based_time - (actframe); % actframe in time (s)
    switch deformation_flag
        case 1 % If there is a deformation annotation set that as activation time
            Tstamps_final = Timestamps - Timestamps(IRM_deformation_frame);
        case 0 % If there are no deformation annotations set time based on landing
            Tstamps_final = Timestamps - Timestamps(first_frame_landed);
    end

    %----------------------------------------------------------
    % Generate Video
    %----------------------------------------------------------
    close (vidObj);
    open(vidObj);

    clc
    close all
    figure
    for i = 1:numframes
        imshow(data(:,:,i),[0 1], 'Border', 'tight','InitialMagnification',200)
        %h = gcf;
        %set(h, 'visible', 'off');

        hold on;

        %----------------------------------------------------------
        % Scale Bar:
        %----------------------------------------------------------
        scaleLine([],10,0.157,[450 172],'hor',[1 1 1],4,channel_text_font_size,0); % 10 um, with text

        %----------------------------------------------------------
        % Timestamp:
        %----------------------------------------------------------
        time_str = KLS_format_seconds_to_time_string(Tstamps_final(i), Tstamps_final(end));
        switch deformation_flag
            case 1 % If there is a deformation annotation put '[]' around the time
                x = offset_from_left;
                y = 8;
                text(x,y,['[' time_str ']'],'color','k','fontsize',channel_text_font_size)
            case 0 % else no brackets around the timestamp
                x = offset_from_left;
                y = 8;
                text(x,y,time_str,'color','k','fontsize',channel_text_font_size)
        end

        %----------------------------------------------------------
        % Channel labels:
        %----------------------------------------------------------
        for ii = 1:Channels
            try
                if contains(lower(channel_labels{Reorder_channels(ii)}), '405') || contains(lower(channel_labels{Reorder_channels(ii)}), 'blue')
                    text_color = '#007bff'; % Blue
                elseif contains(lower(channel_labels{Reorder_channels(ii)}), '561') || contains(lower(channel_labels{Reorder_channels(ii)}), 'dnd')
                    text_color = '#ffad00'; % Orange
                elseif contains(lower(channel_labels{Reorder_channels(ii)}), '647') || contains(lower(channel_labels{Reorder_channels(ii)}), '640') || contains(lower(channel_labels{Reorder_channels(ii)}), 'deep red')
                    text_color = '#ff0000'; % Red
                else
                    text_color = channel_colors{Reorder_channels(ii)};
                    %text_color = '#ffffff'; % white
                end    
                x = offset_from_left+(size(data,2)/Channels)*(ii-1);
                y = size(data,1)-offset_from_bot;
                text(x,y,channel_labels{Reorder_channels(ii)},'fontsize',channel_text_font_size,'color',text_color,'rot',0,'FontWeight','bold')
            catch
                try
                    x = offset_from_left+(size(data,2)/Channels)*(ii-1);
                    y = size(data,1)-offset_from_bot;
                    text(x,y,channel_labels{Reorder_channels(ii)},'fontsize',channel_text_font_size,'color',text_color,'rot',0,'FontWeight','bold')
                catch
                    x = offset_from_left+(size(data,2)/Channels)*(ii-1);
                    y = size(data,1)-offset_from_bot;
                    text(x,y,['Ch' sprintf('%02d', ii)],'fontsize',channel_text_font_size,'color','#ffffff','rot',0,'FontWeight','bold')
                end
            end
        end

        %----------------------------------------------------------
        % Vertical Line between Timepoints:
        %----------------------------------------------------------
        for ii = 1:Channels-1
            line([(size(data,2)/Channels)*ii+1 (size(data,2)/Channels)*ii+1],[0 size(data,1)],'Color',[1 1 1],'LineWidth',2);
        end
        
        %----------------------------------------------------------
        % Add all cell boundary to all channels
        %----------------------------------------------------------
        % Works with multi ROIs in one Image
        [B,~] = bwboundaries(cell_mask(:,:,i) == ROI_n,'noholes');
        hold on
        for k = 1:length(B)
           boundary = B{k};
           x_outline = boundary(:,2);
           y_outline = boundary(:,1);
           for ii = 1:Channels
               x = x_outline+(size(data,2)/Channels*(ii-1));
               y = y_outline;
               plot(x,y,'Color',"#00FFFF",'LineStyle',':','LineWidth',2) % Third Channel
           end
        end   
        hold off

        clear x y
        %----------------------------------------------------------
        % Add Impulse localizations and Tracks to indicated channels
        %----------------------------------------------------------
        if no_Impulse_tracks_flag == 0
            % Binding Localizations
            hold on
            if i <= size(STLN_Tracks,2)
                if ~all(isnan(STLN_Tracks(:,i,1)))                   
                    [B,~] = bwboundaries(cell_mask(:,:,i) == ROI_n,'noholes');
                    hold on
                    for k = 1:length(B)
                       boundary = B{k};
                       x_outline = boundary(:,2);
                       y_outline = boundary(:,1);

                        x = STLN_Tracks(:,i,1);
                        y = STLN_Tracks(:,i,2);

                        % TrackMate indexs at 0,0. Matlab at 1,1;
                        x = x+1; 
                        y = y+1;

                        % Which localizations are inside the cell contact
                            % area?
                        [in,~] = inpolygon(x,y, x_outline, y_outline);
                        for ii = find(Reorder_channels == Impulse_track_annotation_channels)
                            x = x(in) + ((size(data,2)/Channels)*(ii-1));
                            y = y(in);
                            scatter(x,y,85,'o','MarkerEdgeColor',"#00FFFF",'LineWidth',1)
                        end
                    end 
                end
            end
            hold off

            
            % Define max_fade_in_t and identify frames to loop over using unique timepoints
            max_fade_in_t = 10; % # frame to trace back in time (in frames)
            curr_idx = find(unique_impulse_timepoints <= i,1,'last');
            first_idx = max([1 (curr_idx - max_fade_in_t + 1)]);
            frames_to_loop_over = unique_impulse_timepoints(first_idx:curr_idx);
            
            % Check if there are any valid frames to look for links
            if exist('in','var')
                % Use the same loop logic as before
                for nn_idx = length(frames_to_loop_over):-1:2
                    nn = frames_to_loop_over(nn_idx);
                    mm = frames_to_loop_over(nn_idx-1);
                    if nn >= size(STLN_Tracks,2)
                        continue
                    end
                    % Current frame localization
                    x_c = STLN_Tracks(:,nn,1) .* in;
                    y_c = STLN_Tracks(:,nn,2) .* in;
                    % Last frame localization
                    x_l = STLN_Tracks(:,mm,1) .* in;
                    y_l = STLN_Tracks(:,mm,2) .* in;
            
                    for track_n = 1:size(x_c,1)
                        x = [x_c(track_n,1) x_l(track_n,1)];
                        y = [y_c(track_n,1) y_l(track_n,1)];
                        if ~isnan(x)
                            %zz = find(flip(unique_impulse_timepoints) == nn);
                            hold on
                            for ii = find(Reorder_channels == Impulse_track_annotation_channels)
                                x_tracks = x + ((size(data,2)/Channels)*(ii-1)) + 1;
                                y_tracks = y + 1;
                                line_opacity = 1-((-1/(1/max_fade_in_t))*(nn_idx/max_fade_in_t)+abs(-1/(1/max_fade_in_t)))/max_fade_in_t;
                                plot(x_tracks,y_tracks, ...
                                    'Color',[0 1 1 line_opacity], ...
                                    'LineStyle','-','LineWidth',3) % fade old tracks
                            end
                            hold off
                        end
                    end
                end
            end     
        end
        
        clear x y x_c y_c x_l y_l in

        % Various figure details:
        h = gcf;     
        %screen_size = get(0,'ScreenSize');
        set(gcf,'Color',[0 0 0])
        %set(gcf, 'Position', [100 200 1150 screen_size(4)*0.28]);
        set(gca,'Visible','off')
        set(gcf,'NextPlot','add');
        set(gcf,'InvertHardCopy','off');
        set(gcf,'PaperPositionMode','auto')
        currFrame = getframe(h);
    
        writeVideo(vidObj, currFrame);
        %close(h)
        clf
    end
    close(h)
    close(vidObj);
    %}
end