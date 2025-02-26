function [no_output_gen_image] = IRCE_MaskVideo_Post(Save_individual_acq_dir, Specific_ROIs, Dilated_label_ROIs)
    cd(Save_individual_acq_dir)

    if isempty(Specific_ROIs)
        Specific_ROIs = 1:size(Dilated_label_ROIs,1);
    end

    %---------------------------------------------------------%
    % Generate Label Video to Review
    %---------------------------------------------------------%
    cd(Save_individual_acq_dir)
    folderName = 'Label_Videos_Post_PreManual_Edits'; % Generate a folder name for this data set
    if ~exist(fullfile(cd, folderName), 'dir')
        % If folder doesn't exist, create it
        mkdir(folderName);
    end
    % Move to the folder
    cd(folderName);

    % #001 --- Generate label video
    n = 1;
    while  n <= width(Specific_ROIs)
        t = 1;
        
        i = Specific_ROIs(n);
        disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])

        ch1_label = Dilated_label_ROIs{i,1};

        max_LUT = max(ch1_label,[],'all');

        if max_LUT == 0
           n = n+1;
           continue;
        end

        % Generate video object
        img_name = ['Label_Img_Post_PreManual_Edits_ROI_' num2str(i,'%03.f')];

        L = ch1_label(:,:,t);
        % convert label image to rgb
        RGB2 = label2rgb(L, prism(max_LUT),'k');
        % show that color coded labeled image
        imshow(RGB2, 'Border', 'tight','InitialMagnification',200)

        hold on 
        s = regionprops(L, 'centroid');
        for ii = 1:numel(s)
            c = s(ii).Centroid; % Get the centroid coordinates
            
            % 20240506 Attempt to speed up labeling when there is a
            % label with a large number, i.e. 999
                % 20240519 KLS - this works. There is propably a better
                % method ...
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

        % Save the displayed figure as PNG
        saveas(gcf, [img_name '.png'], 'png');
        
        % Save the displayed figure as JPEG
        %saveas(gcf, [Montage_Name '.jpeg'], 'jpeg');
        %saveas(gcf, [img_name '.fig'], 'fig');

        h = gcf;
        close(h)

        n = n+1;
    end

    close all
end