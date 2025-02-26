function IRCE_gen_DIA_Video_ROI_Corners(DIA_stack, Save_individual_acq_dir, DIA_LUT, Specific_ROIs)
    cd(Save_individual_acq_dir)
    if exist('roi_corners.mat','file')
        load('roi_corners','-mat')
    else
        disp('You need to generate ROIs corners before you visualize them.')
        return;
    end
    
    if bitor(max(Specific_ROIs) > length(roi_corners), min(Specific_ROIs) < 1)
        disp(['Specific_ROIs to be edited must be between 1 and ' num2str(length(roi_corners))])
        return;
    end
    
    folderName = 'ROI_Corners_Video_Revew'; % Generate a folder name for this data set
    if ~exist(fullfile(cd, folderName), 'dir')
        % If folder doesn't exist, create it
        mkdir(folderName);
    end
    % Move to the folder
    cd(folderName);
    
    for i = Specific_ROIs
        
        disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);
        img = DIA_stack(min(y):max(y),min(x):max(x),:);
                   

        video_name = ['Review_ROICorners_' num2str(i,'%03.f')];
        vidObj = VideoWriter([video_name '.mp4'], 'MPEG-4');
        vidObj.FrameRate = round(min(size(img,3) / 10, 120)); % dynamically set the frame rate to have a video length of 10s maximum
        vidObj.Quality = 100; % Optional: Adjust video quality if needed 

        close (vidObj);
        open(vidObj);
        %clc
        numframes = size(img,3);

        if exist('h','var')
            close(h)
        end
        for t = 1:numframes
            imshow(img(:,:,t),DIA_LUT, 'Border', 'tight','InitialMagnification',200)
            
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
        end

        close(vidObj);
    end

    sound(exp(sin(1:1500))) 
end