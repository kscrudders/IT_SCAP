function [roi_IRM_interface] = IRCE_IRMDeformationAnnotation_Automated(Save_individual_acq_dir, roi_corners, Ch1_corr_IRM, Base_label_ROIs, IRM_LUT, Dilated_label_ROIs,IRM_thres)
    %---------------------------------------------------------%
    % Maunally Edit Masks -- Ninja Slice ROIs apart
    %---------------------------------------------------------%
    cd(Save_individual_acq_dir)
    if exist('roi_mask_IRM_Deformation.mat','file')
        load('roi_mask_IRM_Deformation','-mat')
    else
        num_possible_deformations_per_frame = 10;
        % size of roi_IRM_interface = [# ROIs, # frames, # of deformations]
        roi_IRM_interface = cell([size(roi_corners,1) size(Ch1_corr_IRM,3) num_possible_deformations_per_frame]);

        %---------------------------------------------------------%
        % loop over each cell ROI
        %---------------------------------------------------------%
        for i = 1:size(roi_corners,1) % Parallel loop over manually selected ROIs
            clc
            disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])
            pause(0.5)
            %Grab current ROI pixels: x and y
            x = roi_corners{i,1}(:,1);
            y = roi_corners{i,1}(:,2);
            img = Ch1_corr_IRM(min(y):max(y),min(x):max(x),:);

            %---------------------------------------------------------%
            % Loop over every frame
            %---------------------------------------------------------%
            close all    
            figure()

            ii = 1;
            while ii <= size(Ch1_corr_IRM,3) % Loop over time
                % skip the current frame if the current cell ROI is absent
                if isempty(find(Base_label_ROIs{i,1}(:,:,ii) == i,1))
                    ii = ii+1;
                    continue
                end
                hold off
                imshow(img(:,:,ii),IRM_LUT)
                g = gcf;
                g.WindowState = 'maximized';            

                % Works with multi ROIs in one Image
                % Add all cell boundary to ch2 and ch3 (impulse and response) 
                [B,~] = bwboundaries(Dilated_label_ROIs{i,1}(:,:,ii) == i,'noholes');
                hold on
                for k = 1:length(B)
                   boundary = B{k};
                   plot(boundary(:,2), boundary(:,1),'Color',"#00FFFF",'LineStyle',':','LineWidth',2)
                end
                hold off

                try
                    disp('Generate a deformation annotation this frame? (k = yes, l = no, ; = go back 1 frame): ')
                    MakeADeformation = getkey(1);
                    MakeADeformation = lower(char(MakeADeformation));
                    
                    % Loop until the input is 'k','l',';'
                    while ~any(strcmp(MakeADeformation,{'k','l',';'}))
                        disp('Invalid input. Please enter ''k - yes'', ''l - no'' or '';''.');
                        
                        disp('Generate a deformation annotation this frame? (k = yes, l = no, ; = go back 1 frame): ')
                        MakeADeformation = getkey(1);
                        MakeADeformation = lower(char(MakeADeformation));
                    end

                    n = 1; % first deformation
                    while MakeADeformation == 'k'
                        temp_line_handle = drawline('Color','y');
                        roi_IRM_interface{i,ii,n} = temp_line_handle.Position;

                        n = n+1;

                        disp('Generate a deformation annotation this frame? (k = yes, l = no, ; = go back 1 frame): ')
                        MakeADeformation = getkey(1);
                        MakeADeformation = lower(char(MakeADeformation));
                        % Loop until the input is 'k','l',';'
                        while ~any(strcmp(MakeADeformation,{'k','l',';'}))
                            disp('Invalid input. Please enter ''k - yes'', ''l - no'' or '';''.');
                            
                            disp('Generate a deformation annotation this frame? (k = yes, l = no, ; = go back 1 frame): ')
                            MakeADeformation = getkey(1);
                            MakeADeformation = lower(char(MakeADeformation));
                        end
                    end   
                catch
                    MakeADeformation = lower(input('Generate a deformation annotation this frame? (k = yes, l = no, ; = go back 1 frame): ','s'));
                    % Loop until the input is 'k','l',';'
                    while ~any(strcmp(MakeADeformation,{'k','l',';'}))
                        disp('Invalid input. Please enter ''k - yes'', ''l - no'' or '';''.');
                        MakeADeformation = lower(input('Generate a deformation annotation this frame? (k = yes, l = no, ; = go back 1 frame): ','s'));
                    end

                    n = 1; % first deformation
                    while MakeADeformation == 'k'
                        temp_line_handle = drawline('Color','y');
                        roi_IRM_interface{i,ii,n} = temp_line_handle.Position;

                        n = n+1;

                        MakeADeformation = lower(input('Generate an additional deformation annotation this frame? (k = yes, l = no, ; = go back 1 frame): ','s'));
                        % Loop until the input is 'k','l',';'
                        while ~any(strcmp(MakeADeformation,{'k','l',';'}))
                            disp('Invalid input. Please enter ''k - yes'', ''l - no'' or '';''.');
                            MakeADeformation = lower(input('Generate an additional deformation annotation this frame? (k = yes, l = no, ; = go back 1 frame): ','s'));
                        end
                    end    
                end
                
                if MakeADeformation == ';'
                    ii = ii-2;
                end
                ii = ii+1;
            end
            save('roi_mask_IRM_Deformation.mat','roi_IRM_interface','-v7.3')
        end
    end

    
    % if more rois have been added since last time annotations were made
    if size(roi_IRM_interface,1) < size(roi_corners,1)
        num_possible_deformations_per_frame = 10;
        % size of roi_IRM_interface = [# ROIs, # frames, # of deformations]
        roi_IRM_interface = cell([size(roi_corners,1) size(Ch1_corr_IRM,3) num_possible_deformations_per_frame]);

        %---------------------------------------------------------%
        % loop over each cell ROI
        %---------------------------------------------------------%
        for i = size(roi_mask_IRM_Deformation,1)+1 : size(roi_corners,1) % Parallel loop over manually selected ROIs
            clc
            disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])
            pause(0.5)
            %Grab current ROI pixels: x and y
            x = roi_corners{i,1}(:,1);
            y = roi_corners{i,1}(:,2);
            img = Ch1_corr_IRM(min(y):max(y),min(x):max(x),:);

            %---------------------------------------------------------%
            % Loop over every frame
            %---------------------------------------------------------%
            close all    
            figure()

            ii = 1;
            while ii <= size(Ch1_corr_IRM,3) % Loop over time
                % skip the current frame if the current cell ROI is absent
                if isempty(find(Base_label_ROIs{i,1}(:,:,ii) == i,1))
                    ii = ii+1;
                    continue
                end

                imshow(img(:,:,ii),IRM_LUT)
                g = gcf;
                g.WindowState = 'maximized';            

                % Works with multi ROIs in one Image
                % Add all cell boundary to ch2 and ch3 (impulse and response) 
                [B,~] = bwboundaries(Dilated_label_ROIs{i,1}(:,:,ii) == i,'noholes');
                hold on
                for k = 1:length(B)
                   boundary = B{k};
                   plot(boundary(:,2), boundary(:,1),'Color',"#00FFFF",'LineStyle',':','LineWidth',2)
                end
                hold off

                MakeADeformation = lower(input('Generate a deformation annotation this frame? (k = yes, l = no, ; = go back 1 frame): ','s'));
                % Loop until the input is 'k','l',';'
                while ~any(strcmp(MakeADeformation,{'k','l',';'}))
                    disp('Invalid input. Please enter ''k - yes'', ''l - no'' or '';''.');
                    MakeADeformation = lower(input('Generate a deformation annotation this frame? (k = yes, l = no, ; = go back 1 frame): ','s'));
                end

                n = 1; % first deformation
                while MakeADeformation == 'k'
                    temp_line_handle = drawline('Color','y');
                    roi_IRM_interface{i,ii,n} = temp_line_handle.Position;

                    n = n+1;

                    MakeADeformation = lower(input('Generate an additional deformation annotation this frame? (k = yes, l = no, ; = go back 1 frame): ','s'));
                    % Loop until the input is 'k','l',';'
                    while ~any(strcmp(MakeADeformation,{'k','l',';'}))
                        disp('Invalid input. Please enter ''k - yes'', ''l - no'' or '';''.');
                        MakeADeformation = lower(input('Generate an additional deformation annotation this frame? (k = yes, l = no, ; = go back 1 frame): ','s'));
                    end
                end

                if MakeADeformation == ';'
                    ii = ii-2;
                end
                ii = ii+1;
            end
            save('roi_mask_IRM_Deformation.mat','roi_IRM_interface','-v7.3')
        end 
    end
    
    clc
    close all
    sound(exp(sin(1:1500)))
end