function [roi_handdrawn_mask_Post, Base_label_ROIs, Dilated_label_ROIs] = IRCE_HandDrawMask_Post(Specific_ROIs, roi_mask_Manual_Crops, roi_corners, Ch1_corr_IRM, IRM_LUT, Base_label_ROIs, Dilated_label_ROIs, Save_individual_acq_dir)
    tic    
    cd(Save_individual_acq_dir)
    
    if exist('roi_handdrawn_mask_Post.mat','file')
        load('roi_handdrawn_mask_Post.mat','roi_handdrawn_mask_Post')
    else
        roi_handdrawn_mask_Post = cell(size(roi_mask_Manual_Crops));
    end
    
    
    n = 1;
    while  n <= width(Specific_ROIs)
        i = Specific_ROIs(n);

        %Grab current ROI pixels: x and y
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);
        img = Ch1_corr_IRM(min(y):max(y),min(x):max(x),:);

        % Basic Segment of that cell ROI in time
        Base_label = zeros(size(img));

        % Populate the Non-edited Label for each cell ROI
        Base_label_ROIs{i,1} = Base_label;

        MakeAMask = 'n'; % initiate make mask check variable
        approve_remaining_flag = 0;
        ii = 1;
        while ii <= size(Ch1_corr_IRM,3) % Loop over time
            %masked_img = img(:,:,ii) .* Base_label(:,:,ii);
            close all
            
            if ii > 1 && all(img(:,:,ii-1) == img(:,:,ii),'all')
                roi_handdrawn_mask_Post{i,ii,1} = roi_handdrawn_mask_Post{i,ii-1,1};
                Base_label_ROIs{i,1}(:,:,ii) = roi_handdrawn_mask_Post{i,ii,1};
                
                ii = ii+1;
                continue
            end
            
            if approve_remaining_flag == 0 % skip correcting mask if user wants to auto approve remaining frames
                figure()
                %imshow(masked_img,IRM_LUT)
                imshow(img(:,:,ii),IRM_LUT, 'Border', 'tight','InitialMagnification',400)

                %g = gcf;
                %g.WindowState = 'maximized';

               try % Try using getKey
                    disp('Generate a Mask this frame? (k = yes, l = no, ; = go back 1 frame, a = skip remaining frames): ');
                    MakeAMask = getkey(1);
                    MakeAMask = lower(char(MakeAMask));

                    % Loop until the input is 'k','l'
                    while ~any(strcmp(MakeAMask, {'k','l',';','a'}))
                        disp('Invalid input. Please enter ''k - yes'', ''l - no'', '';'' or ''a''.');
                        disp('Generate a Mask this frame? (k = yes, l = no, ; = go back 1 frame, a = skip remaining frames): ');
                        MakeAMask = getkey(1);
                        MakeAMask = lower(char(MakeAMask));
                    end

                catch % if getKey not present, then use default Matlab functions
                    MakeAMask = lower(input('Generate a Mask this frame? (k = yes, l = no, ; = go back 1 frame, a = skip remaining frames): ','s'));
                    % Loop until the input is 'k','l',';'
                    while ~any(strcmp(MakeAMask, {'k','l',';','a'}))
                        disp('Invalid input. Please enter ''k - yes'', ''l - no'', '';'' or ''a''.');
                        MakeAMask = lower(input('Generate a Mask this frame? (k = yes, l = no, ; = go back 1 frame, a = skip remaining frames): ','s'));
                    end
               end

                if MakeAMask == 'k'
                    roi_interface{i,ii,1} = drawfreehand('color','w','Closed',1);
                    roi_handdrawn_mask_Post{i,ii,1} = createMask(roi_interface{i,ii,1});
                    Base_label_ROIs{i,1}(:,:,ii) = roi_handdrawn_mask_Post{i,ii,1};
                else

                    roi_handdrawn_mask_Post{i,ii,1} = zeros(size(Base_label,[1 2]));
                    Base_label_ROIs{i,1}(:,:,ii) = roi_handdrawn_mask_Post{i,ii,1};
                end
                
                if MakeAMask == 'a'
                    approve_remaining_flag = 1;
                end
            else
                    roi_handdrawn_mask_Post{i,ii,1} = zeros(size(Base_label,[1 2]));
                    Base_label_ROIs{i,1}(:,:,ii) = roi_handdrawn_mask_Post{i,ii,1};
            end
            
            if MakeAMask == ';'
               ii = ii - 2;
                switch ii
                    case -1
                        ii = 0;
                    case 0
                        ii = 0;
                    otherwise
                        while all(Base_label_ROIs{i,1}(:,:,ii+1) == Base_label_ROIs{i,1}(:,:,ii),'all')
                            ii = ii-1;
                            if ii == 0 
                                break
                            end
                        end

                        if all(Base_label_ROIs{i,1}(:,:,ii+1) ~= Base_label_ROIs{i,1}(:,:,ii),'all')
                            while all(Base_label_ROIs{i,1}(:,:,ii+1) == Base_label_ROIs{i,1}(:,:,ii),'all')
                                ii = ii-1;
                                if ii == 0 
                                    break
                                end
                            end
                        end
                end
            end
                    
            ii = ii+1;
        end
        
        save('roi_handdrawn_mask_Post.mat','roi_handdrawn_mask_Post','-v7.3')

        n = n+1;
    end

    %---------------------------------------------------------%
    % Dilate the base label to yeild the dilated mask
    %---------------------------------------------------------%
    SE_3 = strel('disk', 6); % Flat Structuring Element for Image dilation, disk of size 9

    tic
    n = 1;
    while  n <= width(Specific_ROIs)
        i = Specific_ROIs(n);
        ii = 1;
        while ii <= size(Ch1_corr_IRM,3)
            Dilated_label_ROIs{i,1}(:,:,ii) = imdilate(Base_label_ROIs{i,1}(:,:,ii),SE_3);
            ii = ii+1;
        end
        n = n+1;
    end
    toc_01 = toc;
    disp(['Labels dilated in ' num2str(toc_01/60) ' min'])
end