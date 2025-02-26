function IRCE_gen_Nch_montage(num_montage_frames, Time_stamps_address, ROI_n, ...
    Individual_ROI_data, Dilated_label_forROIs, channel_freq, channel_LUTs, ...
    channel_labels, channel_colors, roi_IRM_interface, STLN_Tracks, ...
    ch_with_tracks, Montage_Name, Reorder_channels)

    if isempty(STLN_Tracks)
        no_tracks_flag = 1;
    else
        no_tracks_flag = 0;
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


    % Double check that the number of montage frames is =< than the
    % data itself
    if num_montage_frames > maxT
        num_montage_frames = size(resized_data{1,1},3);
    end
    
    last_frame = maxT;
    min_frame_start = last_frame - num_montage_frames;    
    
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

    channel_text_font_size = 12;
    %annotation_text_font_size = 14;
    offset_from_bot = 12;
    offset_from_left = 5;

    % Determine maximum label width in pixels
    label_widths = zeros(Channels,1);
    for ii = 1:Channels
        label = channel_labels{Reorder_channels(ii)};
        h_fig = figure('visible','off');
        h_ax = axes('Parent',h_fig,'Units','pixels');
        h_text = text(0,0,label,'FontSize', channel_text_font_size,'FontWeight','bold','Units','pixels');
        extent = get(h_text,'Extent');
        label_widths(ii) = extent(3);
        close(h_fig);
    end
    max_label_width = max(label_widths);

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
    % Combine Individual dip_image channel variables
    % into one dip_image variable
    %----------------------------------------------------------
    montage_height = size(normalize_data{1,1},1) * Channels;
    montage_width = size(normalize_data{1,1},2) * num_montage_frames;
    montage_data = zeros([montage_height, montage_width]);

    % Figure out what frames go into the montage
    if first_frame_landed > min_frame_start
        % If cell landing is less than the number of 
        % montage frames from the end, idx from the end to then end minus
        % the number of montage frames
        montage_idx = min_frame_start+1:last_frame;
    else
        % If landing is early enough in time, go 1 frame
        % before landing and proceed
        if first_frame_landed > 1
            montage_idx = round(linspace(first_frame_landed-1, last_frame, num_montage_frames));
        else
            montage_idx = round(linspace(first_frame_landed, last_frame, num_montage_frames));
        end
    end
    
    % Populate the montage
    for i = 1:num_montage_frames
        for ii = 1:Channels
            switch ii
                case 1
                    montage_data(1:size(normalize_data{1,1},1), ...
                        1+(size(normalize_data{1,1},2)*(i-1)):(size(normalize_data{1,1},2)*(i))) ...
                        = normalize_data{Reorder_channels(ii),1}(:,:,montage_idx(i));
                otherwise
                    montage_data(size(normalize_data{1,1},1)*(ii-1)+1:size(normalize_data{1,1},1)*(ii), ...
                        1+(size(normalize_data{1,1},2)*(i-1)):(size(normalize_data{1,1},2)*(i))) ...
                        = normalize_data{Reorder_channels(ii),1}(:,:,montage_idx(i));
            end
        end
   end

    % WIP: Find defining way to document the activation frame
    % At present the definition is the first IRM deformation annotation
    % Timestamps = zero_based_time - (actframe); % actframe in time (s)
    
    switch deformation_flag
        case 1 % If there is a deformation annotation set that as activation time
            Tstamps_final = Timestamps - Timestamps(IRM_deformation_frame);
        case 0 % If there are no deformation annotations set time based on landing
            Tstamps_final = Timestamps - Timestamps(first_frame_landed);
    end

    %----------------------------------------------------------
    % Generate Video
    %----------------------------------------------------------
    clc
    % Check if the longest channel labels fits in the first column of the
        % montage, if not, upscale the image so that the label fits
    if (max_label_width+offset_from_left) > size(normalize_data{1,1},2)
        image_scaling_factor = 100 * ((max_label_width+offset_from_left*2) / size(normalize_data{1,1},2));
        image_scaling_factor = ceil(image_scaling_factor);
    else
        image_scaling_factor = 100;
    end
    imshow(montage_data(:,:),[0 1], 'Border', 'tight','InitialMagnification',image_scaling_factor)

    hold on;
    %----------------------------------------------------------
    % Scale Bar:
    %----------------------------------------------------------
    scaleLine([],10,0.157,[450 172],'hor',[1 1 1],4,channel_text_font_size,0); % 10 um, with text
    
    %----------------------------------------------------------
    % Timestamp:
    %----------------------------------------------------------
    for i = 1:num_montage_frames
        time_str = KLS_format_seconds_to_time_string(Tstamps_final(montage_idx(i)), Tstamps_final(end));
        switch deformation_flag
            case 1 % If there is a deformation annotation put '[]' around the time
                x = offset_from_left+((i-1)*size(normalize_data{1,1},2));
                y = 8;
                text(x,y,['[' time_str ']'],'color','k','fontsize',channel_text_font_size)
            case 0 % else no brackets around the timestamp
                x = offset_from_left+((i-1)*size(normalize_data{1,1},2));
                y = 8;
                text(x,y,time_str,'color','k','fontsize',channel_text_font_size)
        end
    end

    %----------------------------------------------------------
    % Channel labels:
    %----------------------------------------------------------
    for ii = 1:Channels
        switch ii
            case 1
                try
                    text_color = channel_colors{Reorder_channels(ii)};

                    x = offset_from_left;
                    y = size(normalize_data{1,1},1)-offset_from_bot;
                    text(x,y,channel_labels{Reorder_channels(ii)},'fontsize',channel_text_font_size,'color',text_color,'rot',0,'FontWeight','bold')
                catch
                    x = offset_from_left;
                    y = size(normalize_data{1,1},1)-offset_from_bot;
                    text(x,y,['Ch' sprintf('%02d', Reorder_channels(ii))],'fontsize',channel_text_font_size,'color',[0 0 0],'rot',0,'FontWeight','bold')
                end
            otherwise
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
                    x = offset_from_left;
                    y = (size(normalize_data{1,1},1)*ii)-offset_from_bot;
                    text(x,y,channel_labels{Reorder_channels(ii)},'fontsize',channel_text_font_size,'color',text_color,'rot',0,'FontWeight','bold')
                catch
                    try
                        x = offset_from_left;
                        y = (size(normalize_data{1,1},1)*ii)-offset_from_bot;
                        text(x,y,channel_labels{Reorder_channels(ii)},'fontsize',channel_text_font_size,'color','#ffffff','rot',0,'FontWeight','bold')
                    catch
                        x = offset_from_left;
                        y = (size(normalize_data{1,1},1)*ii)-offset_from_bot;
                        text(x,y,['Ch' sprintf('%02d', ii)],'fontsize',channel_text_font_size,'color','#ffffff','rot',0,'FontWeight','bold')
                    end
                end
        end
    end

    %{
    % CAR GFP level label:
    text((size(data,1)/Channels)*3-offset_from_left,size(data,2)-offset_from_bot,['GFP = ' num2str(mean_eGFP)],'fontsize',channel_text_font_size/2,'color','green','rot',0,'FontWeight','bold','HorizontalAlignment','right')
    %}

    % Annotations:
    % Activation annoation
    %{
    if Tstamps(:,i) > landing_ann_start-actframe && Tstamps(:,i) < landing_ann_end-actframe 
        text((size(data,1)/Channels)/2,size(data,2)*3/4,'Landing','Color',[1 1 1],'FontSize',annotation_text_font_size,'FontWeight','bold','HorizontalAlignment','center')
        % Create arrow:
    %         annotation(gcf,'arrow',[0.07 0.1],[0.4 0.51],'linewidth',2,'Color',[1 1 1]);
    end
    %}

    % Polarization annoation
    %{
    if Tstamps(:,i) > polarization_ann_start-actframe && Tstamps(:,i) < polarization_ann_end-actframe 
        text((size(data,1)/Channels)/2+(size(data,1)/Channels)*2,size(data,2)*3/4,'Polarization','Color',[1 1 0],'FontSize',annotation_text_font_size,'FontWeight','bold','HorizontalAlignment','center')
    %         annotation(gcf,'arrow',[0.79 0.82],[0.4 0.51],'linewidth',2,'Color',[1 1 0]);
    end    
    %}

    %----------------------------------------------------------
    % Vertical Line between Timepoints:
    %----------------------------------------------------------
    for i = 1:num_montage_frames-1
        x = [(size(montage_data,2)/num_montage_frames)*i-1 (size(montage_data,2)/num_montage_frames)*i-1];
        y = [0 size(montage_data,1)];
        line(x,y,'Color',[1 1 1],'LineWidth',2);
    end
    %----------------------------------------------------------
    % Horizontal Line between Channels
    %----------------------------------------------------------
    for i = 1:Channels-1
        x = [1 size(montage_data,2)];
        y = [(size(montage_data,1)/Channels)*i+1 (size(montage_data,1)/Channels)*i+1];
        line(x,y,'Color',[1 1 1],'LineWidth',2);
    end
    

    %----------------------------------------------------------
    % Add all cell boundary to all channels
    %----------------------------------------------------------
    % Works with multi ROIs in one Image
    for i = 1:num_montage_frames
        [B,~] = bwboundaries(cell_mask(:,:,montage_idx(i)) == ROI_n,'noholes');
        hold on
        for k = 1:length(B)
           boundary = B{k};
           x_outline = boundary(:,2)+((i-1)*size(normalize_data{1,1},2));
           y_outline = boundary(:,1)+size(normalize_data{1,1},1);
           for ii = 1:Channels
               x = x_outline;
               y = y_outline+(size(normalize_data{1,1},1)*(ii-1));
               plot(x,y,'Color',"#00FFFF",'LineStyle',':','LineWidth',2) % Third Channel
           end
        end   
    end
    hold off

    %----------------------------------------------------------
    % Add localizations to preferred channels
    %----------------------------------------------------------
    if no_tracks_flag == 0
        % Binding Localizations
        hold on

        for i = 1:num_montage_frames
            if montage_idx(i) <= size(STLN_Tracks,2)
                if ~all(isnan(STLN_Tracks(:,montage_idx(i),1)))                    
                    [B,~] = bwboundaries(cell_mask(:,:,montage_idx(i)) == ROI_n,'noholes');
                    hold on
                    for k = 1:length(B)
                       boundary = B{k};
                       x_outline = boundary(:,2);
                       y_outline = boundary(:,1);

                        x = STLN_Tracks(:,montage_idx(i),1);
                        y = STLN_Tracks(:,montage_idx(i),2);

                        % TrackMate indexs at 0,0. Matlab at 1,1;
                        x = x+1; 
                        y = y+1;

                        % Which localizations are inside the cell contact
                            % area?
                        [in,~] = inpolygon(x,y, x_outline, y_outline);
                        for ii = find(Reorder_channels == ch_with_tracks)
                            x = x(in) + ((i-1)*size(normalize_data{1,1},2));
                            y = y(in) + (size(normalize_data{1,1},1)*(ii-1));
                            scatter(x,y,85,'o','MarkerEdgeColor',"#00FFFF",'LineWidth',1)
                        end
                    end 
                end
            end
        end
        hold off
    end

    % Various figure details:

    %screen_size = get(0,'ScreenSize');
    set(gcf,'Color',[0 0 0])
    %set(gcf, 'Position', [100 200 screen_size(3)*0.95 size(data,1)*2]);
    set(gca,'Visible','off')
    set(gcf,'NextPlot','add');
    set(gcf,'InvertHardCopy','off');
    set(gcf,'PaperPositionMode','auto')

    % Save the displayed figure as PNG
    saveas(gcf, [Montage_Name '.png'], 'png');

    % Save the displayed figure as JPEG
    %saveas(gcf, [Montage_Name '.jpeg'], 'jpeg');
    saveas(gcf, [Montage_Name '.fig'], 'fig');
    h = gcf;
    close(h)
end