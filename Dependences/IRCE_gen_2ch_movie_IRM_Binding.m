function [no_output_gen_video] = IRCE_gen_2ch_movie_IRM_Binding(Dilated_label_forROIs, Ch1_corr_IRM, IRM_LUT, Ch2_Impulse, Binding_LUT, Time_stamps_address, Video_Name, ROI_n,channel_labels, roi_IRM_interface)
    Data.IRM = Ch1_corr_IRM;
    Data.Binding = Ch2_Impulse;

    %----------------------------------------------------------
    % Check if Ch1, Ch2 needs to be resized to the largest data in t
    %----------------------------------------------------------
    maxT = max([size(Data.IRM,3), size(Data.Binding,3)]);
    Data.IRM = KLS_resizeMatrix(Data.IRM, maxT);
    Dilated_label_forROIs = KLS_resizeMatrix(Dilated_label_forROIs, maxT);
    Data.Binding = KLS_resizeMatrix(Data.Binding, maxT);
    
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
    ch3_data = Data.Binding;
    cell_mask = Dilated_label_forROIs;

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
    
    Channels = 2;
    ch1_LUT = IRM_LUT; % Contact Area - IRM
    ch3_LUT = Binding_LUT; % Response - LysoBrigth Blue

    % Generate video object
    vidObj = VideoWriter([Video_Name '.mp4'], 'MPEG-4');
    vidObj.FrameRate = round(min(maxT / 10, 120)); % dynamically set the frame rate to have a video length of 10s maximum
    vidObj.Quality = 100;
    close(vidObj);

    channel_text_font_size = 12;
    %annotation_text_font_size = 14;
    offset_from_bot = 6;
    offset_from_left = 5;

    zero_based_time = Timing_seconds;

    numframes = size(Data.IRM,3);

    % Reorder endpoint to get this information before the movie, then add 
    %mean_eGFP = 1204; % Manual epi GFP expression pull from ImageJ, (590 bgkd e488)

    %----------------------------------------------------------
    % KLS normalize each channel's LUT 0-1
    %----------------------------------------------------------
    ch1_data = ch1_data-ch1_LUT(1); % set minimum
    ch1_data  = ch1_data./(ch1_LUT(2)-ch1_LUT(1)); % set max to 1, account for lower bound shift

    ch3_data = ch3_data-ch3_LUT(1); % set minimum
    ch3_data = ch3_data./(ch3_LUT(2)-ch3_LUT(1)); % set max to 1, account for lower bound shift

    %----------------------------------------------------------
    % Combine Individual dip_image channel variables
    % into one dip_image variable
    %----------------------------------------------------------
    data = zeros([size(Data.IRM,1), size(Data.IRM,2)*Channels, size(Data.IRM,3)]);

    data(1:size(Data.IRM,1), 1:size(Data.IRM,2), :) = ch1_data;
    data(1:size(Data.IRM,1), 1*(size(Data.IRM,2))+1:2*(size(Data.IRM,2)), :) = ch3_data;

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
    close (vidObj);
    open(vidObj);

    clc
    close all
    figure
    for i = 1:numframes
        imshow(data(:,:,i),[0 1], 'Border', 'tight','InitialMagnification',300)
        %h = gcf;
        %set(h, 'visible', 'off');

        hold on;

        % Scale Bar:
        scaleLine([],10,0.157,[450 172],'hor',[1 1 1],4,channel_text_font_size,0); % 10 um, with text
        %diptruesize(gcf,180)

        % Timestamp:
        switch deformation_flag
            case 1 % If there is a deformation annotation put '[]' around the time
                text(offset_from_left,8,['[' num2str(Tstamps_final(:,i),'%.2f') '] min'],'color','k','fontsize',channel_text_font_size)
            case 0 % else no brackets around the timestamp
                text(offset_from_left,8,[num2str(Tstamps_final(:,i),'%.2f') ' min'],'color','k','fontsize',channel_text_font_size)
        end
        % Channel labels:
        % Ch1
        try
            text(offset_from_left+(size(data,2)/Channels)*0,size(data,1)-offset_from_bot,channel_labels.Ch1,'fontsize',channel_text_font_size,'color',[0 0 0],'rot',0,'FontWeight','bold')
        catch
            text(offset_from_left+(size(data,2)/Channels)*0,size(data,1)-offset_from_bot,'Ch1','fontsize',channel_text_font_size,'color',[0 0 0],'rot',0,'FontWeight','bold')
        end
        %Ch3
        try
            if contains(lower(channel_labels.Ch2), '405')
                text_color = '#007bff'; % Blue
            elseif contains(lower(channel_labels.Ch2), '561')
                text_color = '#ffad00'; % Orange
            elseif contains(lower(channel_labels.Ch2), '647') || contains(lower(channel_labels.Ch3), '640')
                text_color = '#ff0000'; % Red
            else
                text_color = '#ffffff'; % white
            end    
            text(offset_from_left+(size(data,2)/Channels)*1,size(data,1)-offset_from_bot,channel_labels.Ch2,'fontsize',channel_text_font_size,'color',text_color,'rot',0,'FontWeight','bold')
        catch
            text(offset_from_left+(size(data,2)/Channels)*1,size(data,1)-offset_from_bot,'Ch2','fontsize',channel_text_font_size,'color','#c6ff00','rot',0,'FontWeight','bold')
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

        % Line between channels:
        line([size(data,2)/2+1 size(data,2)/2+1],[0 size(data,1)],'Color',[1 1 1],'LineWidth',2);

        % Works with multi ROIs in one Image
        % Add all cell boundary to ch2 and ch3 (impulse and response) 
        [B,~] = bwboundaries(cell_mask(:,:,i) == ROI_n,'noholes');
        hold on
        for k = 1:length(B)
           boundary = B{k};
           %plot(boundary(:,2)+(size(data,2)/Channels)*0, boundary(:,1),'Color',"#00FFFF",'LineStyle',':','LineWidth',2)
           plot(boundary(:,2)+(size(data,2)/Channels)*1, boundary(:,1),'Color',"#00FFFF",'LineStyle',':','LineWidth',2)
        end
        hold off

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
    end
    close(h)
    close(vidObj);
    %}
end