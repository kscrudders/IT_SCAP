function IRCE_gen_IRM_Video_ROI_Corners(Ch1_corr_IRM, Save_individual_acq_dir, IRM_thres, IRM_LUT, Specific_ROIs)
    cd(Save_individual_acq_dir)
    if exist('roi_corners.mat','file')
        load('roi_corners','-mat')
    else
        disp('You need to generate ROIs corners before you edit them.')
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

    % Basic Segment of that cell ROI in time
    Base_label = Ch1_corr_IRM<IRM_thres;
    
    %---------------------------------------------------------%
    % Clean up that basic Segmentation
    %---------------------------------------------------------%
    % 304 px^2 ~= 7.5 µm^2
    P = 150; 

    small_obj_removed = zeros(size(Ch1_corr_IRM));
    t = 1;
    while t <= size(Base_label,3)
        CC = bwconncomp(Base_label(:,:,t), 8); % 8-way connectivity
        S = regionprops(CC, 'Area');
        L = labelmatrix(CC);
        small_obj_removed(:,:,t) = ismember(L, find([S.Area] >= P));
        t = t+1;
    end
    small_obj_removed = small_obj_removed>0;

    SE = strel('disk', 6); % Flat Structuring Element for Image dilation, disk of size 9
    SE_2 = strel('disk', 7); % Flat Structuring Element for Image dilation, disk of size 2

    t = 1;
    while t <= size(Base_label,3)
        Base_label(:,:,t) = imdilate(small_obj_removed(:,:,t),SE);
        t = t+1;
    end

    t = 1;
    while t <= size(Base_label,3)
        Base_label(:,:,t) = imfill(Base_label(:,:,t),'holes');
        t = t+1;
    end


    t = 1;
    while t <= size(Base_label,3)
        Base_label(:,:,t) = imerode(Base_label(:,:,t),SE_2);
        t = t+1;
    end

    
    for i = Specific_ROIs
        
        disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);
        img = Ch1_corr_IRM(min(y):max(y),min(x):max(x),:);
                    
        ch1_label = Base_label(min(y):max(y),min(x):max(x),:);

        video_name = ['Review_ROICorners_' num2str(i,'%03.f')];
        vidObj = VideoWriter([video_name '.mp4'], 'MPEG-4');
        vidObj.FrameRate = round(min(ceil(size(ch1_label,3) / 10), 120)); % dynamically set the frame rate to have a video length of 10s maximum
        vidObj.Quality = 100; % Optional: Adjust video quality if needed 

        close (vidObj);
        open(vidObj);
        %clc
        numframes = size(ch1_label,3);

        if exist('h','var')
            close(h)
        end
        for t = 1:numframes
            imshow(img(:,:,t),IRM_LUT, 'Border', 'tight','InitialMagnification',200)
            curr_label = ch1_label(:,:,t);
            
            [rows, cols] = size(curr_label);

            L = curr_label;
    
            [B,~] = bwboundaries(L>0,'noholes');
            hold on
            for k = 1:length(B)
                boundary = B{k};

                % Check if the boundary touches the edge of the image
                if any(boundary(:,1) == 1) || any(boundary(:,1) == rows) || any(boundary(:,2) == 1) || any(boundary(:,2) == cols)
                    % Boundary touches the edge, color red
                    plot(boundary(:,2), boundary(:,1), 'Color', 'red', 'LineStyle', ':', 'LineWidth', 2)
                else
                    % Boundary does not touch the edge, color yellow
                    plot(boundary(:,2), boundary(:,1), 'Color', 'cyan', 'LineStyle', ':', 'LineWidth', 2)
                end
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
        end

        close(vidObj);
    end
end