function [no_output_gen_video] = IRCE_FinalMaskVideo_IRM(Save_individual_acq_dir, Specific_ROIs, roi_mask_Manual_Crops, Dilated_label_ROIs, processed_data, channel_LUTs, Contact_Channel, roi_corners)
cd(Save_individual_acq_dir)

if isempty(Specific_ROIs)
    Specific_ROIs = 1:size(roi_mask_Manual_Crops,1);
end

tic
%---------------------------------------------------------%
% Generate Label Video to Review
%---------------------------------------------------------%
cd(Save_individual_acq_dir)
    folderName = 'Label_Videos_PostManual_Edits'; % Generate a folder name for this data set
    if ~exist(fullfile(cd, folderName), 'dir')
        % If folder doesn't exist, create it
        mkdir(folderName);
    end
    % Move to the folder
    cd(folderName);

    % #001 --- Generate label video
    n = 1;
    while  n <= width(Specific_ROIs)
        i = Specific_ROIs(n);
        disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])

        ch1_label = Dilated_label_ROIs{i,1};
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);
        ch1_IRM = processed_data{Contact_Channel,1}(min(y):max(y),min(x):max(x),:);

        min_LUT = channel_LUTs{Contact_Channel,1}(1);
        max_LUT = channel_LUTs{Contact_Channel,1}(2);

        video_name = ['Label_Videos_PostManual_Edits_ROI_' num2str(i,'%03.f')];
        vidObj = VideoWriter([video_name '.mp4'], 'MPEG-4');
        vidObj.FrameRate = round(min(size(ch1_label,3) / 10, 120)); % dynamically set the frame rate to have a video length of 10s maximum
        vidObj.Quality = 100; % Optional: Adjust video quality if needed 

        close (vidObj);
        open(vidObj);
        clc
        numframes = size(ch1_label,3);

        if exist('h','var')
            close(h)
        end
        for t = 1:numframes
            L = ch1_label(:,:,t);
            img = ch1_IRM(:,:,t);

            % show that color coded labeled image
            imshow(img,[min_LUT max_LUT],'Border', 'tight','InitialMagnification',100)

            %----------------------------------------------------------
            % Add all cell boundary to all channels
            %----------------------------------------------------------
            % Works with multi ROIs in one Image
            [B,~] = bwboundaries(L > 0,'noholes');
            hold on
            for k = 1:length(B)
               boundary = B{k};
               x_outline = boundary(:,2);
               y_outline = boundary(:,1);
               x = x_outline;
               y = y_outline;
               plot(x,y,'Color',"#00FFFF",'LineStyle',':','LineWidth',2) % Third Channel
            end   

            hold on 
            s = regionprops(L, 'centroid');
            for ii = 1:numel(s)
                c = s(ii).Centroid; % Get the centroid coordinates
                
                % 20240506 Attempt to speed up labeling when there is a
                % label with a large number, i.e. 999
                if isnan(c(1))
                    continue
                end
                text(c(1), c(2), sprintf('%d', ii), 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', 'Color', 'black','FontSize',14);
            end
            % Various figure details:
            %screen_size = get(0,'ScreenSize');
            set(gcf,'Color',[0 0 0])
            %set(gcf, 'Position', [100 200 1150 screen_size(4)*0.28]);
            set(gca,'Visible','off')
            set(gcf,'NextPlot','add');
            set(gcf,'InvertHardCopy','off');
            set(gcf,'PaperPositionMode','auto')
            h = gcf;
            currFrame = getframe(h);

            writeVideo(vidObj, currFrame);
            clf
        end

        close(vidObj);
        n = n+1;
    end

    close all
end
