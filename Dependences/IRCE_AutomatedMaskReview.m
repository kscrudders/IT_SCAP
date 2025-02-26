function [multi_mask_list, mask_gaps_list] = IRCE_AutomatedMaskReview(Specific_ROIs, Dilated_label_ROIs)
% Last update: 20241203 KLS
%     - Corrected error if a ROI does not have a label

    if isempty(Specific_ROIs)
        Specific_ROIs = 1:size(Dilated_label_ROIs,1);
    end
    
    no_mask_list = nan([1 size(Dilated_label_ROIs,1)]);
    idx_nomask = 1;
    one_mask_list = nan([1 size(Dilated_label_ROIs,1)]);
    idx_onemask = 1;
    multi_mask_list = nan([1 size(Dilated_label_ROIs,1)]);
    idx_multimask = 1;
    
    % Initialize a structure to store ROIs and labels with gaps
    labels_with_gaps_list = struct('ROI', {}, 'Labels', {});
    
    %---------------------------------------------------------%
    % Loop through each cell ROIs to find out how many masks there are in
    % that piece of data
    %---------------------------------------------------------%

    n = 1;
    while  n <= width(Specific_ROIs)
        i = Specific_ROIs(n);
        %disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])

        ch1_label = Dilated_label_ROIs{i,1};
        curr_labels = unique(ch1_label(:));
        curr_valid_labels = [0 i 999];

        %---------------------------------------------------------%
        % Check for gaps in label occurrence over time
        %---------------------------------------------------------%

        % Exclude labels 0 and 999
        labels_to_check = setdiff(curr_labels, [0, 999]);

        % Initialize a list to store labels with gaps
        labels_with_gaps = [];

        % Loop over each label
        for l = labels_to_check'
            if isempty(l)
                continue;
            end
            % Get the time points where label l occurs
            label_mask = (ch1_label == l);
            time_points_l = find(squeeze(any(any(label_mask,1),2)));
            
            % Check for gaps between first and last occurrence
            if ~isempty(time_points_l)
                t_min = min(time_points_l);
                t_max = max(time_points_l);
                t_full = t_min:t_max;
                
                if length(t_full) ~= length(time_points_l) || any(~ismember(t_full, time_points_l))
                    % There are gaps
                    labels_with_gaps = [labels_with_gaps; l];
                end
            end
        end

        % If any labels have gaps, collect the ROI number and labels
        if ~isempty(labels_with_gaps)
            labels_with_gaps_list(end+1).ROI = i;
            labels_with_gaps_list(end).Labels = labels_with_gaps;
        end

        max_LUT = max(curr_labels,[],'all');

        switch max_LUT
            case 0
                no_mask_list(idx_nomask) = i;
                idx_nomask = idx_nomask+1;        
            case 1
                one_mask_list(idx_onemask) = i;
                idx_onemask = idx_onemask+1;   
            otherwise 
                m = length(curr_labels);
                switch m
                    case 2 % Check if the current ROI has only labels 0 and i
                        if sum(ismember(curr_labels, curr_valid_labels(1:2))) == 2
                            one_mask_list(idx_onemask) = i;
                            idx_onemask = idx_onemask+1;   
                        else
                            multi_mask_list(idx_multimask) = i;
                            idx_multimask = idx_multimask+1;   
                        end
                    case 3 % Check if the current ROI has only labels 0, i and 999
                        if sum(ismember(curr_labels, curr_valid_labels)) == 3
                            one_mask_list(idx_onemask) = i;
                            idx_onemask = idx_onemask+1;   
                        else
                            multi_mask_list(idx_multimask) = i;
                            idx_multimask = idx_multimask+1;   
                        end
                    otherwise
                        multi_mask_list(idx_multimask) = i;
                        idx_multimask = idx_multimask+1; 
                end
        end
        n = n+1;
    end
    
    clc
    disp(['Number of ROIs = ' num2str(size(Dilated_label_ROIs,1))])

    % Check if any ROIs lack a mask
    if any(~isnan(no_mask_list))
        disp(' ')
        disp('ROIs with no mask:')
        no_mask_list = no_mask_list(~isnan(no_mask_list));
        disp(no_mask_list)   
    end

    % Check if any ROIs have only one mask (this is good)
    if isnan(one_mask_list)
        disp(' ')
        disp('ROIs with only one masks = 0')
    else
        disp(' ')
        disp('ROIs with only one mask:')
        one_mask_list = one_mask_list(~isnan(one_mask_list));
        disp(one_mask_list)
    end

    disp(' ')

    % Check if any ROIs have 2 or more masks (this is bad)
    % You need to manual identify which mask is the relevent cell
    if any(~isnan(multi_mask_list))
        disp(' ')
        disp('ROIs with 2 or more masks:')
        multi_mask_list = multi_mask_list(~isnan(multi_mask_list));
        disp(multi_mask_list)
    end

    %---------------------------------------------------------%
    % Display ROIs with labels that have gaps in time
    %---------------------------------------------------------%

    if ~isempty(labels_with_gaps_list)
        mask_gaps_list = zeros([length(labels_with_gaps_list) 1]);
        disp(' ')
        disp('ROIs with labels that have gaps in time:')
        for idx = 1:length(labels_with_gaps_list)
            roi_num = labels_with_gaps_list(idx).ROI;
            mask_gaps_list(idx) = roi_num;
            labels_in_gap = labels_with_gaps_list(idx).Labels;
            disp(['ROI ' num2str(roi_num) ', labels with gaps: ' num2str(labels_in_gap')]);
        end
    else
        mask_gaps_list = [];
        disp(' ')
        disp('No ROIs have labels with gaps in time.')
    end

    clear no_mask_list one_mask_list
    close all
end

%{
% Old version of the code 20241013
function multi_mask_list = IRCE_AutomatedMaskReview(Specific_ROIs, Dilated_label_ROIs)
    
    if isempty(Specific_ROIs)
        Specific_ROIs = 1:size(Dilated_label_ROIs,1);
    end
    
    no_mask_list = nan([1 size(Dilated_label_ROIs,1)]);
    idx_nomask = 1;
    one_mask_list = nan([1 size(Dilated_label_ROIs,1)]);
    idx_onemask = 1;
    multi_mask_list = nan([1 size(Dilated_label_ROIs,1)]);
    idx_multimask = 1;
    
    %---------------------------------------------------------%
    % Loop through each cell ROIs to find out how many masks there are in
    % that piece of data
    %---------------------------------------------------------%

    n = 1;
    while  n <= width(Specific_ROIs)
        i = Specific_ROIs(n);
        %disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])

        ch1_label = Dilated_label_ROIs{i,1};
        curr_labels = unique(ch1_label(:));
        curr_valid_labels = [0 i 999];

        max_LUT = max(curr_labels,[],'all');

        switch max_LUT
            case 0
                no_mask_list(idx_nomask) = i;
                idx_nomask = idx_nomask+1;        
            case 1
                one_mask_list(idx_onemask) = i;
                idx_onemask = idx_onemask+1;   
            otherwise 
                m = length(curr_labels);
                switch m
                    case 2 % Check if the current ROI has only labels 0 and i
                        if sum(ismember(curr_labels, curr_valid_labels(1:2))) == 2
                            one_mask_list(idx_onemask) = i;
                            idx_onemask = idx_onemask+1;   
                        else
                            multi_mask_list(idx_multimask) = i;
                            idx_multimask = idx_multimask+1;   
                        end
                    case 3 % Check if the current ROI has only labels 0, i and 999
                        if sum(ismember(curr_labels, curr_valid_labels)) == 3
                            one_mask_list(idx_onemask) = i;
                            idx_onemask = idx_onemask+1;   
                        else
                            multi_mask_list(idx_multimask) = i;
                            idx_multimask = idx_multimask+1;   
                        end
                    otherwise
                        multi_mask_list(idx_multimask) = i;
                        idx_multimask = idx_multimask+1; 
                end
        end
        n = n+1;
    end
    
    clc
    disp(['Number of ROIs = ' num2str(size(Dilated_label_ROIs,1))])

    % Check if any ROIs lack a mask
    if any(~isnan(no_mask_list))
        disp(' ')
        disp('ROIs with no mask:')
        no_mask_list = no_mask_list(~isnan(no_mask_list));
        disp(no_mask_list)   
    end


    % Check if any ROIs have only one mask (this is good)
    if isnan(one_mask_list)
        disp(' ')
        disp('ROIs with only one masks = 0')
    else
        disp(' ')
        disp('ROIs with only one mask:')
        one_mask_list = one_mask_list(~isnan(one_mask_list));
        disp(one_mask_list)
    end

    disp(' ')

    % Check if any ROIs have 2 or more masks (this is bad)
    % You need to manual identify which mask is the relevent cell
    if any(~isnan(multi_mask_list))
        disp(' ')
        disp('ROIs with 2 or more masks:')
        multi_mask_list = multi_mask_list(~isnan(multi_mask_list));
        disp(multi_mask_list)
    end

    clear no_mask_list one_mask_list
    close all
end
%}