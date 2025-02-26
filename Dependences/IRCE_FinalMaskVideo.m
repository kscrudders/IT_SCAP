function [no_output_gen_video] = IRCE_FinalMaskVideo(Save_individual_acq_dir, Specific_ROIs, roi_mask_Manual_Crops, Dilated_label_ROIs)
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

        min_LUT = 1;
        max_LUT = max(ch1_label,[],'all');

        if max_LUT == 0
           n = n+1;
           continue;
        end

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
            % convert label image to rgb
            RGB2 = label2rgb(L,prism(max_LUT),'k');
            % show that color coded labeled image
            imshow(RGB2, 'Border', 'tight','InitialMagnification',200)

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
            %close(h)
            clf
        end

        close(vidObj);

        toc
        n = n+1;
    end

    close all
end
