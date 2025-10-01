function [no_output_gen_video] = IRCE_gen_IRMResponse_movie(Dilated_label_forROIs, Ch1_corr_IRM, IRM_LUT, Ch3_Response, Response_LUT, Time_stamps, Video_Name, ROI_n, roi_IRM_interface,channel_labels)
    Data.IRM = Ch1_corr_IRM;
    Data.Lyso = Ch3_Response;
    
    %----------------------------------------------------------
    % Check if Ch1 or Ch2 needs to be remapped to get paired
    % data for each frame for each channel
    %----------------------------------------------------------
    maxT = max([size(Data.IRM,3), size(Data.Lyso,3)]);
    Data.IRM = KLS_resizeMatrix(Data.IRM, maxT);
    Data.Lyso = KLS_resizeMatrix(Data.Lyso, maxT);
    Dilated_label_forROIs = KLS_resizeMatrix(Dilated_label_forROIs, maxT);
    roi_IRM_interface = KLS_resizeAnnotationsCell(roi_IRM_interface, maxT);
    
    Timing_seconds = Time_stamps;

    ch1_data = Data.IRM;
    ch3_data = Data.Lyso;
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
    ch3_LUT = Response_LUT; % Response - LysoBrigth Blue

    % Generate video object
    vidObj = VideoWriter([Video_Name '.mp4'], 'MPEG-4');
    vidObj.FrameRate = round(min(ceil(size(ch1_label,3) / 10), 120)); % dynamically set the frame rate to have a video length of 10s maximum
    vidObj.Quality = 100; % Optional: Adjust video quality if needed

    %vidObj.FileFormat = 'mp4';
    close(vidObj);

    channel_text_font_size = 24;
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
    if exist('h','var')
        close(h)
    end
    close all
    figure
    for i = 1:numframes
        hold off
        imshow(data(:,:,i),[0 1], 'Border', 'tight','InitialMagnification',300)
        %h = gcf;
        %set(h, 'visible', 'off');

        hold on;

        % Scale Bar:
        scaleLine([],10,0.157,[450 172],'hor',[1 1 1],4,channel_text_font_size,0); % 10 um, with text
        %diptruesize(gcf,180)

        % Timestamp:
        %text(offset_from_left,8,[num2str(Tstamps(:,i),'%.2f') ' [' num2str(Tstamps_final(:,i),'%.2f') '] min'],'color','k','fontsize',channel_text_font_size)
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
            text(offset_from_left+(size(data,2)/Channels)*1,size(data,1)-offset_from_bot,channel_labels.Ch3,'fontsize',channel_text_font_size,'color',text_color,'rot',0,'FontWeight','bold')
        catch
            text(offset_from_left+(size(data,2)/Channels)*1,size(data,1)-offset_from_bot,'Ch3','fontsize',channel_text_font_size,'color','#c6ff00','rot',0,'FontWeight','bold')
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
        line([size(data,2)/2 size(data,2)/2],[0 size(data,1)],'Color',[1 1 1],'LineWidth',2);

        % Works with multi ROIs in one Image
        % Add all cell boundary to ch2 and ch3 (impulse and response) 
        [B,~] = bwboundaries(cell_mask(:,:,i) == ROI_n,'noholes');
        hold on
        for k = 1:length(B)
           boundary = B{k};
           plot(boundary(:,2)+(size(data,2)/Channels)*0, boundary(:,1),'Color',"#00FFFF",'LineStyle',':','LineWidth',2)
           plot(boundary(:,2)+(size(data,2)/Channels)*1, boundary(:,1),'Color',"#00FFFF",'LineStyle',':','LineWidth',2)
        end
        hold off
        
        % IRM Deformation annoation
        for k = 1:size(roi_IRM_interface,3)
            current_cicle = roi_IRM_interface{1,i,k}; % [x1 y1; x2 y2]
            if isempty(current_cicle)
                break
            end
            x = (current_cicle(1)+current_cicle(2))/2;
            y = (current_cicle(3)+current_cicle(4))/2;
            r = sqrt(((current_cicle(2)-current_cicle(1))^2+(current_cicle(4)-current_cicle(3))^2))/2;
            viscircles([x+1 y+1],r,'Color','#FAFA33','LineStyle','--','Linewidth',2,'EnhanceVisibility',0);
            viscircles([x+1+(size(data,2)/Channels)*1 y+1],r,'Color','#FAFA33','LineStyle','--','Linewidth',2,'EnhanceVisibility',0);
        end

        % Various figure details:
        h = gcf;     
        %screen_size = get(0,'ScreenSize');
        %img_width = size(data,2);
        %img_height = size(data,1);
        %set(gcf, 'OuterPosition', [100, 100, img_width*2, img_height*2+33.5*2]);
        %set(gcf, 'OuterPosition', [100, 100, 720, 720*(img_height+33.5)/img_width]);
        set(gcf,'Color',[0 0 0])        
        set(gca,'Visible','off')
        set(gcf,'NextPlot','add');
        set(gcf,'InvertHardCopy','off');
        set(gcf,'PaperPositionMode','auto')
        currFrame = getframe(h);

        writeVideo(vidObj, currFrame);
        %close(h)
    end

    close(vidObj);
    close all
    %}
end
