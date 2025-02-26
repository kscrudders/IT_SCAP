function [no_output_gen_montage] = IRCE_gen_3ch_montage(Dilated_label_forROIs,Ch1_corr_IRM, IRM_LUT, Ch2_Impulse, Impulse_LUT, Ch3_Response, Response_LUT, Time_stamps_address, Montage_Name, ROI_n, channel_labels, STLN_Tracks, roi_IRM_interface)
    if nargin == 11
        STLN_Tracks = [];
        no_tracks_flag = 1;
    else
        no_tracks_flag = 0;
    end
        
    num_montage_frames = 10;

    %----------------------------------------------------------
    % Check if Ch1, Ch2 or Ch3 needs to be resized to the largest data in t
    %----------------------------------------------------------
    maxT = max([size(Ch1_corr_IRM,3), size(Ch3_Response,3), size(Ch2_Impulse,3)]);
    Data.IRM = KLS_resizeMatrix(Ch1_corr_IRM, maxT);
    Dilated_label_forROIs = KLS_resizeMatrix(Dilated_label_forROIs, maxT);
    Data.Lyso = KLS_resizeMatrix(Ch3_Response, maxT);
    Data.Binding = KLS_resizeMatrix(Ch2_Impulse, maxT);
    
    % Double check that the number of montage frames is not longer than the
    % data itself
    if num_montage_frames > size(Data.IRM,3)
        num_montage_frames = size(Data.IRM,3);
    end
    
    last_frame = size(Dilated_label_forROIs,3);
    min_frame_start = last_frame - num_montage_frames;    
    
    
    % Pull in the recorded timepoints
    fileID = fopen(Time_stamps_address, 'r');

    % Assuming time format in the file is like 12:34.567 (mm:ss.ms)
    % '%d' reads an integer, ':%d.%f' reads the seconds and milliseconds
    data = fscanf(fileID, '%d:%d.%d', [3, Inf]);

    switch floor(max(log10(data(3,:))))
        case -Inf
                % Convert data to total seconds
                % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
                Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 10;
        case 0
            % Convert data to total seconds
            % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
            Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 10;    
        case 1
            % Convert data to total seconds
            % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
            Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 100;                   
        case 2
            % Convert data to total seconds
            % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
            Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 1000;                   
        case 3
            % Convert data to total seconds
            % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
            Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 10000;                   
    end
    
    % Close the file
    fclose(fileID);


    ch1_data = Data.IRM;
    ch2_data = Data.Binding;
    ch3_data = Data.Lyso;
    cell_mask = Dilated_label_forROIs;

    % Identify when the cell lands (ie when the masking starts)
    size_vector = KLS_label2ROIsize(cell_mask == ROI_n, 0.157);
    first_frame_landed = find(size_vector > 0,1,'first'); % first frame the cell landed

    % First frame showing a IRM deformation annotation
    deformation_flag = 0;
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
    
    Channels = 3;
    ch1_LUT = IRM_LUT; % Contact Area - IRM
    ch2_LUT = Impulse_LUT; % Impulse - FOLR1-AF647 LL640
    ch3_LUT = Response_LUT; % Response - LysoBrite Blue

    channel_text_font_size = 12;
    %annotation_text_font_size = 14;
    offset_from_bot = 6;
    offset_from_left = 5;

    zero_based_time = Timing_seconds;


    % Reorder endpoint to get this information before the movie, then add 
    %mean_eGFP = 1204; % Manual epi GFP expression pull from ImageJ, (590 bgkd e488)

    %----------------------------------------------------------
    % KLS normalize each channel's LUT 0-1
    %----------------------------------------------------------
    ch1_data = ch1_data-ch1_LUT(1); % set minimum
    ch1_data  = ch1_data./(ch1_LUT(2)-ch1_LUT(1)); % set max to 1, account for lower bound shift

    ch2_data = ch2_data-ch2_LUT(1); % set minimum
    ch2_data = ch2_data./(ch2_LUT(2)-ch2_LUT(1)); % set max to 1, account for lower bound shift
    
    ch3_data = ch3_data-ch3_LUT(1); % set minimum
    ch3_data = ch3_data./(ch3_LUT(2)-ch3_LUT(1)); % set max to 1, account for lower bound shift

    %----------------------------------------------------------
    % Combine Individual dip_image channel variables
    % into one dip_image variable
    %----------------------------------------------------------
    montage_height = size(Data.IRM,1)*Channels;
    montage_width = size(Data.IRM,2)*num_montage_frames;
    data = zeros([montage_height, montage_width]);

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
        data(1:size(Data.IRM,1), 1+(size(Data.IRM,2)*(i-1)):(size(Data.IRM,2)*(i))) = ch1_data(:,:,montage_idx(i));
        data(size(Data.IRM,1)+1:size(Data.IRM,1)*2, 1+(size(Data.IRM,2)*(i-1)):(size(Data.IRM,2)*(i))) = ch2_data(:,:,montage_idx(i));
        data(size(Data.IRM,1)*2+1:size(Data.IRM,1)*3, 1+(size(Data.IRM,2)*(i-1)):(size(Data.IRM,2)*(i))) = ch3_data(:,:,montage_idx(i));
    end

    % Find defining way to document the activation frame
    %Timestamps = zero_based_time - (actframe); % actframe in time (s)
    Timestamps = zero_based_time; % actframe in time (s)
    Tstamps = Timestamps/60; % convert s to min
    
    switch deformation_flag
        case 1 % If there is a deformation annotation set that as activation time
            Tstamps_final = Tstamps - Tstamps(IRM_deformation_frame);
        case 0 % If there are no deformation annotations set time based on landing
            Tstamps_final = Tstamps - Tstamps(first_frame_landed);
    end
    %{
    % #001 --- Save Data Variables
    Data.IRM_mask = base_IRM_mask;

    cd(ResDir)
    % Check if the folder 'Data_save' exists in the current directory
    if ~exist('Data_save', 'dir')
        % If it doesn't exist, create the folder
        mkdir('Data_save');
    end
    % Move into 'Data_save' folder
    cd('Data_save');

    save([filestr(1:end-4) '_Data_save'],'Data','DataDir','ResDir','filestr','Timing_seconds' ...
        ,'Channels','actframe','landing_ann_start','landing_ann_end','polarization_ann_start','polarization_ann_end' ...
        ,'ch1_LUT','ch2_LUT','ch3_LUT','video_name','mean_eGFP')
    %}

    %----------------------------------------------------------
    % Generate Video
    %----------------------------------------------------------
    clc
    imshow(data(:,:),[0 1], 'Border', 'tight','InitialMagnification',200)
    %h = gcf;
    %set(h, 'visible', 'off');

    hold on;
    % Scale Bar:
    scaleLine([],10,0.157,[450 172],'hor',[1 1 1],4,channel_text_font_size,0); % 10 um, with text
    %diptruesize(gcf,180)

    %if size(Data.IRM,2) < 120
    %   channel_text_font_size = 8; 
    %end
    for i = 1:num_montage_frames    
        % Timestamp:
        switch deformation_flag
            case 1 % If there is a deformation annotation put '[]' around the time
                text(offset_from_left+((i-1)*size(Data.IRM,2)),8,['[' num2str(Tstamps_final(:,montage_idx(i)),'%.2f') '] min'],'color','k','fontsize',channel_text_font_size)
            case 0 % else no brackets around the timestamp
                text(offset_from_left+((i-1)*size(Data.IRM,2)),8,[num2str(Tstamps_final(:,montage_idx(i)),'%.2f') ' min'],'color','k','fontsize',channel_text_font_size)
        end
    end

    % Channel labels:
    % Ch1
    try
        text(offset_from_left,size(Data.IRM,1)-offset_from_bot,channel_labels.Ch1,'fontsize',channel_text_font_size,'color',[0 0 0],'rot',0,'FontWeight','bold')
    catch
        text(offset_from_left,size(Data.IRM,1)-offset_from_bot,'Ch1','fontsize',channel_text_font_size,'color',[0 0 0],'rot',0,'FontWeight','bold')
    end
    % Ch2
    try
        if contains(lower(channel_labels.Ch2), '561')
            text_color = '#ffad00'; % Orange
        elseif contains(lower(channel_labels.Ch2), '647') || contains(lower(channel_labels.Ch2), '640')
            text_color = '#ff0000'; % Red
        else
            text_color = '#ffffff'; % white
        end    
        text(offset_from_left,size(Data.IRM,1)*2-offset_from_bot,channel_labels.Ch2,'fontsize',channel_text_font_size,'color',text_color,'rot',0,'FontWeight','bold')
    catch
        text(offset_from_left,size(Data.IRM,1)*2-offset_from_bot,'Ch2','fontsize',channel_text_font_size,'color','#ff2100','rot',0,'FontWeight','bold')
    end
    % Ch3
    try
        if contains(lower(channel_labels.Ch3), 'blue')
            text_color = '#007bff'; % Blue
        elseif contains(lower(channel_labels.Ch3), 'dnd')
            text_color = '#ffad00'; % Orange
        elseif contains(lower(channel_labels.Ch3), 'dr') || contains(lower(channel_labels.Ch3), 'deep red')
            text_color = '#ff0000'; % Red
        else
            text_color = '#ffffff'; % white
        end     
        text(offset_from_left,size(Data.IRM,1)*3-offset_from_bot,channel_labels.Ch3,'fontsize',channel_text_font_size,'color',text_color,'rot',0,'FontWeight','bold')
    catch
        text(offset_from_left,size(Data.IRM,1)*3-offset_from_bot,'Ch3','fontsize',channel_text_font_size,'color','#c6ff00','rot',0,'FontWeight','bold')
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

    % Vertical Line between Timepoints:
    for i = 1:num_montage_frames
        line([(size(data,2)/10)*i-1 (size(data,2)/10)*i-1],[0 size(data,1)],'Color',[1 1 1],'LineWidth',2);
    end
    % Horizontal Line between Channels
    for i = 1:Channels-1
        line([1 size(data,2)],[(size(data,1)/Channels)*i+1 (size(data,1)/Channels)*i+1],'Color',[1 1 1],'LineWidth',2);
    end
    
    % Works with multi ROIs in one Image
    % Add all cell boundary to ch2 and ch3 (impulse and response) 
    for i = 1:num_montage_frames
        [B,~] = bwboundaries(cell_mask(:,:,montage_idx(i)) == ROI_n,'noholes');
        hold on
        for k = 1:length(B)
           boundary = B{k};
           x_outline = boundary(:,2)+((i-1)*size(Data.IRM,2));
           y_outline = boundary(:,1)+size(Data.IRM,1);
           plot(x_outline,y_outline ,'Color',"#00FFFF",'LineStyle',':','LineWidth',2) % Second Channel
           plot(x_outline,y_outline+size(Data.IRM,1) ,'Color',"#00FFFF",'LineStyle',':','LineWidth',2) % Third Channel
        end   
    end
    hold off

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
                       x_outline = boundary(:,2)+((i-1)*size(Data.IRM,2));
                       y_outline = boundary(:,1)+size(Data.IRM,1);

                        x = STLN_Tracks(:,montage_idx(i),1);
                        y = STLN_Tracks(:,montage_idx(i),2);
                        %x = x./0.157; % convert micron to px
                        %y = y./0.157; % convert micron to px
                        %x(x == 0) = nan;
                        %y(y == 0) = nan;

                        x = x+((i-1)*size(Data.IRM,2))+1;
                        y = y+size(Data.IRM,1)+1;
                        [in,~] = inpolygon(x,y, x_outline, y_outline);
                        scatter(x(in),y(in),85,'o','MarkerEdgeColor',"#00FFFF",'LineWidth',1)
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