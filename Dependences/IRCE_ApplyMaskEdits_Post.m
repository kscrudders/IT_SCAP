function [Dilated_label_ROIs,  Base_label_ROIs] = IRCE_ApplyMaskEdits_Post(roi_mask_Manual_Crops, Raw_label_ROIs, Ch1_corr_IRM, Specific_ROIs, Base_label_ROIs)
    if isempty(Specific_ROIs)
        Specific_ROIs = 1:size(roi_mask_Manual_Crops,1);
    end
    %---------------------------------------------------------%
    % Update Basic Labels with manual edits from S02c
    %---------------------------------------------------------%
    SE = strel('disk', 2); % Flat Structuring Element for Image dilation, disk of size 2

    for i = Specific_ROIs
        for ii = 1:size(Ch1_corr_IRM,3)
            Base_label_ROIs{i,1}(:,:,ii) = bitand(Raw_label_ROIs{i,1}(:,:,ii)>0,~imdilate(roi_mask_Manual_Crops{i,ii},SE));
        end                
    end

    for i = Specific_ROIs
        Base_label_ROIs{i,1} = bwlabel(Base_label_ROIs{i,1},4);
    end

    %---------------------------------------------------------%
    % Remove Edge touching ROIs
    %---------------------------------------------------------%

    for i = Specific_ROIs
        Base_label_ROIs{i,1} = KLS_remove_ROI_touching_edge(Base_label_ROIs{i,1});
    end

    %---------------------------------------------------------%
    % Remove small cell traces
    %---------------------------------------------------------%
    pixel_size = 0.157;
    min_size = 0; % minimium contact area a cell must have for the entire trace to be included

    for i = Specific_ROIs
        Base_label_ROIs{i,1} = KLS_RemoveSmallCellLabels(Base_label_ROIs{i,1}, pixel_size, min_size);
    end

    for i = Specific_ROIs
        Base_label_ROIs{i,1} = double(Base_label_ROIs{i,1});
    end
    
    %---------------------------------------------------------%
    % Dilate the base label to yeild the dilated mask
    %---------------------------------------------------------%
    SE_3 = strel('disk', 6); % Flat Structuring Element for Image dilation, disk of size 9
    Dilated_label_ROIs = Base_label_ROIs;

    for i = Specific_ROIs
        ii = 1;
        while ii <= size(Ch1_corr_IRM,3)
            Dilated_label_ROIs{i,1}(:,:,ii) = imdilate(Base_label_ROIs{i,1}(:,:,ii),SE_3);
            ii = ii+1;
        end
    end

    clear i  ii SE SE_2 SE_3 toc_01 min_size pixel_size
    sound(exp(sin(1:1500)))
end