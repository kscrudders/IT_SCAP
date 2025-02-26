%% "file type - data analysis/practice/etc"
% "Title"
%
% "Data source" by "Person genreating data"
%
% Started "date"
% Last updated: "date" by "person updating"
%
%
% Purpose: "purpose"
%
% Colors hexcode references:
    % '#007bff'; % Blue
    % '#ffad00'; % Orange
    % '#ff0000'; % Red
    % '#ffffff'; % white

	% Cobalt Blue - "#1446A0"
	% Razzmatazz - "#DB3069"
	
	% Slate Gray - "#627C85"
	% Imperial Red - "#F71735"
	% Russian Violet - "#2A1E5C"
	
	% Sea Green Crayola - "#44FFD1"	
	% Cyan - "#00FFFF"
% Helpful symbols
	% ± x̄
	% α ß Γ π Σ σ µ τ Φ
%---------------------------------------------------------%
%
%---------------------------------------------------------%


%% Section_00: Import Dependences --
clr

%---------------------------------------------------------%
% Shade data generated from lipophilic dye loaded SLBs
%---------------------------------------------------------%
Dependency_folder = 'H:\...\Support Dependencies'; % <--- Change me as neeeded
load(fullfile(Dependency_folder,'Shade_20240701.mat'),'-mat') % <--- Change me as neeeded

%---------------------------------------------------------%
% Set a default figure position
%---------------------------------------------------------%
Pos = get(0,'defaultfigureposition');

% Get the screen size
screenSize = get(0, 'ScreenSize'); % screenSize = [left, bottom, width, height]

% Define your figure size
figureWidth = Pos(3); % Width of the figure
figureHeight = Pos(4); % Height of the figure

% Calculate the position to place the figure at the far left edge of the screen
Pos = [5, screenSize(4) - figureHeight - 85, figureWidth, figureHeight]; % [left, bottom, width, height]

% Create the figure with the specified position
%figure('Position', Pos);

clear screenSize figureWidth figureHeight

KLS_Check_tif_Imwrite()

clear Dependency_folder

%% Section_01: Import Data -- 
%---------------------------------------------------------%
% Stuff to Change
%---------------------------------------------------------%
% Where is your data? 
data_dir = 'G:\...\FolderWithTargetImageData'; % <--- Change me
% Where do you want to save data?
save_data_in_this_folder = 'G:\...\FolderWithTargetImageData'; % <--- Change me

% What is the name of your data?
name_of_data = 'Example_image'; % <--- Change me

base_save_dir = [save_data_in_this_folder '\' name_of_data];

% What are you time stamps? Needs to be the same size as the image data
    % Raw time stampes from Nikon's NIS elements metadata mm:ss.ms format
try
    file_name = dir(fullfile(base_save_dir, '*_Time_Field*.txt')); 

    Time_stamps_address = fullfile(base_save_dir, file_name.name);
catch
    % If not please put the full path to the timings .txt file here
    disp(['Timestamp import is automatic if you save ' ... 
        'the ND2 time information in a txt file with ' ...
        '''_Time_Field'' in it''s name'])
    Time_stamps_address = ''; % <--- Change Me as needed = 
end

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%
[Raw_data, Save_individual_acq_dir, Raw_dir, Processed_dir, ROIs_dir] = IRCE_ImportFileSetup(data_dir, save_data_in_this_folder, name_of_data);

% Pull out the metadata from the ND2 file
try
    % Try the ND2 file in the data dir
    meta_data = KLS_ParseND2Metadata(fullfile(data_dir, [name_of_data '.nd2']));
catch
    % Sometime the metadata is corrupt, manually generate it
    meta_data = KLS_GetMetadata(base_save_dir);
end

KLS_Check_tif_Imwrite()
clear loc_of_data_in_dir curr_dir file_name

disp(' ')
disp(' ')
disp('Channels:')
for i = 1:length(meta_data)
    disp(['Ch' num2str(i) ' --- ' meta_data(i).Name])
end

clear i

%% Section_02a: Img Process -- Raw Data and  --     
%---------------------------------------------------------%
% Stuff to Change
%---------------------------------------------------------%
close all

num_ch = 3;

% Which channel denotes where the cell contact is?
Contact_Channel = 1;

% Is this Impulse-Response Data? Which channel denotes the response?
    % set to [] if no response
Response_Channel = 2;
Impulse_Channel = 3;

% Make sure that all of the following are

% channel frequencies, list must be num_ch long 
channel_freq = [7 1 4]; 
% channel labels, cell array must be num_ch long 
channel_labels = {'IRM'; 'LysoBrite Blue'; 'FOLR1 Binding';}; % 't647 Short Expo. Binding'
% channel colors, cell array must be num_ch long, hex-code or RGB acceptable 
channel_colors = {'#ffffff'; '#007bff'; '#ff0000';};
    % '#007bff'; % Blue
    % '#00ff00'; % Green
    % '#ffad00'; % Orange
    % '#ff0000'; % Red
    % '#ffffff';   % White
% channel LUT, cell array must be num_ch long, individual LUTs are
    % 2-element row vectors
channel_LUTs = {[3800 7800]; [0 50]; [0 0.75e4];};

% channel numbers that you want later tracking data to be overlayed on
Impulse_track_annotation_channels = [3]; % leave empty if you do not want localization/track overlays
Response_track_annotation_channels = [2]; % leave empty if you do not want localization/track overlays

% binary flag indicating which channels get median image filtered, row vector must be num_ch long 
median_filter_ch_flag = [1 0 0];
% Of the channels being median image filtered, which can serve as their own
    % background data? (i.e. a tiled data set with mutliple FOVs/timepoint), 
    % row vector must be num_ch long 
self_med_filter_ch_flag = [1 0 0];

med_filter_lower_thresholds = {[3 5]; []; [];};

% binary flag indicating which channels get shade corrected, row vector must be num_ch long 
shade_correct_ch_flag = [0 1 1];
% Shade_NDFin is a 4 by 3 cell array
    % row 1-4, are t405, t488, t561, t647
    % column 1-2 are TIRF direction 0°, 45°, 135°
NDFin_or_NDFout = 1; % in = 1, out = 0

% binary flag indicating which channels get converted from arbitrary 
    % fluorescence units to eletron volts (photons), row vector must be num_ch long
AU_to_Photon_flag = [0 1 1];

% binary flag indicating which channels get photobleach corrected using
    % the mean image intensity decay as the data for a single exponential 
    % fit, row vector must be num_ch long
bleach_correct_flag = [0 0 1];

%---------------------------------------------------------%
% If you need an external background reference image complete the details
% for where that data is below
%---------------------------------------------------------%
if any(median_filter_ch_flag .* ~self_med_filter_ch_flag)
    % If you have multiple channels in your background, let the program
    % know which slice is approperate for each data chennal you are median
    % image filtering.
    external_ch_num = [2 0 0]; % Must be a length of num_ch
    %---------------------------------------------------------%
    % Manual Median Image Correct the IRM data
    %---------------------------------------------------------%
    
    % What folder is your IRM bkgd image in?
    cd 'G:\...\FolderWithTargetImageData' % <--- Change me
    data_dir = dir;
    loc_or_name_of_median_data = [25 26 27]; % <--- Change Me (What file do you want to import?) = 
    
    temp_cell = cell([length(loc_or_name_of_median_data) 1]);
    for i = 1:length(loc_or_name_of_median_data)
        temp_cell{i} = KLS_ND2ImportAll(data_dir(loc_or_name_of_median_data(i)).name);
    end
end

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%
processed_data = cell([num_ch 1]); % a cell array {n_channels by 1}
median_threshold = zeros([num_ch 1]);

if exist('Raw_data','var') && ~exist('base_data','var')
    base_data = cell([num_ch 1]); % a cell array {n_channels by 1}

    for i = 1:num_ch
        base_data{i,1} = Raw_data(:,:,i:num_ch*channel_freq(i):end);

        % If this was not taken every frame, check if if you need to remove
        % repeating or empty (all zeros) frames. 
        % Update channel_freq for  the unique image data.
        if channel_freq(i) ~= 1
            [x, y, t] = size(base_data{i,1});
            reshapedData = reshape(base_data{i,1}, x*y, t);
            [uniqueColumns, ~, ~] = unique(reshapedData', 'rows', 'stable');
            uniqueData = reshape(uniqueColumns', x, y, size(uniqueColumns, 1));
            nonZeroFrames = any(any(uniqueData, 1), 2);
            base_data{i,1} = uniqueData(:, :, nonZeroFrames);
            channel_freq(i) = ceil(size(Raw_data,3) / (t*num_ch));
        end
        clear x y t reshapedData uniqueColumns uniqueData nonZeroFrames
    end
    clear Raw_data
end

% Let's do median image correction, AU to photon convertion, shade
% correction, and blech correct as per requested above
for i = 1:num_ch
    processed_data{i,1} = base_data{i,1};

    if median_filter_ch_flag(i) == 1
        switch self_med_filter_ch_flag(i)
            case 1
                [processed_data{i,1}, ~, median_threshold(i), ~] = KLS_RICM_bkgd_correction(processed_data{i,1}, med_filter_lower_thresholds{i});
            case 0
                temp_bkgd = zeros([size(temp_cell{1}, 1) size(temp_cell{1}, 2) length(loc_or_name_of_median_data)]);

                for ii = 1:length(loc_or_name_of_median_data)
                    temp_bkgd(:,:,ii) = temp_cell{ii}(:,:,external_ch_num(i));
                end
                
                [~, median_img, median_threshold(i), ~] = KLS_RICM_bkgd_correction(temp_bkgd, med_filter_lower_thresholds{i});
                % Pull bkgd IRM image from the large data set.
                    % If no large image data available construct a background from different
                    % IRM image data sets taken that day.
                    % Check if processed_data is 512x512 in the first
                    % two dimensions, then resize median_img.
                rows = round(size(processed_data{i,1},1)/512);
                cols = round(size(processed_data{i,1},2)/512);

                if floor(rows) ~= rows || floor(cols) ~= cols
                    disp('Current script only handles image data that is divisible by 512.')
                    return;
                end
                if rows + cols > 2
                    median_img = repmat(median_img,[rows cols 1]);
                end
                
                processed_data{i,1} = (processed_data{i,1} - median_img) + round(mean(median_img,'all'));
        end
    end

    if AU_to_Photon_flag(i) == 1
            gain = meta_data(i).Multiplier;
            [conversion, offset] = KLS_gain_basic(gain);
    
            processed_data{i,1} = (processed_data{i,1}-offset) .* conversion;
    end
    
    if shade_correct_ch_flag(i) == 1
        rows = round(size(base_data{i,1}(:,:,1),1)/512);
        cols = round(size(base_data{i,1}(:,:,1),2)/512);

    if floor(rows) ~= rows || floor(cols) ~= cols
        disp('Current script only handles image data that is divisible by 512.')
        return;
    end
    if isscalar(meta_data(i).ExWavelength)
            switch meta_data(i).ExWavelength
                case 405
                    WL_num = 1;
                case 488
                    WL_num = 2;
                case 561
                    WL_num = 3;
                case 640
                    WL_num = 4;
                otherwise
                    % Check meta_data(i).Name for wavelengths
                    if contains(meta_data(i).Name, '405')
                        WL_num = 1;
                    elseif contains(meta_data(i).Name, '488')
                        WL_num = 2;
                    elseif contains(meta_data(i).Name, '561')
                        WL_num = 3;
                    elseif contains(meta_data(i).Name, '640') || contains(meta_data(i).Name, '647')
                        WL_num = 4;
                    else
                        WL_num = 2; % Default to 488
                        warning('ExWavelength not recognized. Setting WL_num to 2 (488).');
                    end
            end
    else
            % Check meta_data(i).Name for wavelengths
            if contains(meta_data(i).Name, '405')
                WL_num = 1;
            elseif contains(meta_data(i).Name, '488')
                WL_num = 2;
            elseif contains(meta_data(i).Name, '561')
                WL_num = 3;
            elseif contains(meta_data(i).Name, '640') || contains(meta_data(i).Name, '647')
                WL_num = 4;
            else
                WL_num = 2; % Default to 488
                warning('ExWavelength not recognized. Setting WL_num to 2 (488).');
            end
    end
    
    if isempty(meta_data(i).TIRF_Direction)
        Direction_num = 1;
    else
        switch meta_data(i).TIRF_Direction
            case 0
                Direction_num = 1;
            case 45
                Direction_num = 2;
            case 135
                Direction_num = 3;
            otherwise
                Direction_num = 1;
        end
    end
        if NDFin_or_NDFout == 1
            shade_img = Shade_NDFin{WL_num,Direction_num};
        else
            shade_img = Shade_NDFout{WL_num,Direction_num};
        end

        Resized_shade_img = repmat(shade_img,[rows cols]);
    
        ii = 1;
        while ii <= size(processed_data{i,1},3)
            processed_data{i,1}(:,:,ii) = processed_data{i,1}(:,:,ii)./Resized_shade_img;
            ii = ii+1;
        end
    end

    if bleach_correct_flag(i) == 1 && size(processed_data{i,1},3) > 3
        %---------------------------------------------------------%
        % Basic Bleach Correction using Mean Image Intensity Bleaching Rate
        %---------------------------------------------------------%
        figure()
        y = median(processed_data{i,1},[1 2]);
        y = squeeze(unique(y,'stable'));
        x = 0:size(y,1)-1;

        scatter(x,y)
        [fitFunction, gof, fit_str, lifetime_tau] = KLS_Exponentialfit_and_plot(processed_data{i,1}, 1);
        
        xlabel('Frame')
        if AU_to_Photon_flag(i) == 1
            ylabel('Mean Intensity (Photons)')
        else
            ylabel('Mean Intensity (AU)')
        end

        legend('Data',fit_str,'location','best')
        title({['Ch_' num2str(i) ' Bleach Correction'], ['Tau = ' num2str(lifetime_tau,3)]})

        x = 0:size(processed_data{i,1},3)-1;
        y = fitFunction(x);
        y = KLS_NormStack(y);

        for ii = 1:size(processed_data{i,1},3)
            processed_data{i,1}(:,:,ii) = processed_data{i,1}(:,:,ii) ./ y(ii);
        end
    end
end

filtered_IRM_data = processed_data{Contact_Channel,1};
IRM_thres = median_threshold(Contact_Channel);

for i = 1:num_ch
    figure('Position',Pos + (i-1)*[Pos(3) 0 0 0])
    histogram(processed_data{i,1},'Normalization','PDF')
    title(['Histogram: ' channel_labels{i}])

    box off
end

clear i ii AU_to_Photon_flag bleach_correct_flag cols Direction_num gain 
clear median_filter_ch_flag NDFin_or_NDFout offset rows self_med_filter_ch_flag 
clear shade_correct_ch_flag WL_num shade_img 

%% Section_03a: -- Select ROIs --   
close all
Add_more_ROIs_flag = 0; % Want to add more ROIs?

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%

roi_corners = IRCE_CellROIs(filtered_IRM_data, Save_individual_acq_dir, Add_more_ROIs_flag, IRM_thres, channel_LUTs{Contact_Channel});

close all
clear Add_more_ROIs_flag

%% Section_03a(I): [Optional] [Video Output] -- Video to Review ROI Corners --   
close all
Specific_ROIs = [];

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%
if isempty(Specific_ROIs)
    Specific_ROIs = 1:length(roi_corners);
    IRCE_gen_IRM_Video_ROI_Corners(filtered_IRM_data, Save_individual_acq_dir, IRM_thres, channel_LUTs{Contact_Channel}, Specific_ROIs)
else
	IRCE_gen_IRM_Video_ROI_Corners(filtered_IRM_data, Save_individual_acq_dir, IRM_thres, channel_LUTs{Contact_Channel}, Specific_ROIs)
end

close all
clear Specific_ROIs

%% Section_03a(II): [As Needed] -- Edit -- Revise Specific ROIs
close all
Specific_ROIs = [];

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%
filtered_IRM_data = processed_data{Contact_Channel,1};
IRM_thres = median_threshold(Contact_Channel);

if ~isempty(Specific_ROIs)
    roi_corners = IRCE_Edit_CellROIs(filtered_IRM_data, Save_individual_acq_dir, IRM_thres, roi_corners, Specific_ROIs);
end

close all
    IRCE_gen_IRM_Video_ROI_Corners(filtered_IRM_data, Save_individual_acq_dir, IRM_thres, channel_LUTs{Contact_Channel}, Specific_ROIs)

close all

%% Section_03a(III): [As Needed] -- Edit -- Remove Specific ROIs
close all
Specific_ROIs = [];

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%
filtered_IRM_data = processed_data{Contact_Channel,1};

if ~isempty(Specific_ROIs)
    roi_corners = IRCE_Remove_CellROIs(filtered_IRM_data, Save_individual_acq_dir, roi_corners, Specific_ROIs);
end

%% Section_03b: -- Automated Mask Cells & Manually Seperate Connecting Cells --    

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%
Base_label_ROIs = IRCE_MaskSeperation(roi_corners, filtered_IRM_data, median_threshold(Contact_Channel), channel_LUTs{Contact_Channel}, Save_individual_acq_dir);
Base_label_ROIs = IRCE_RemoveSmallCellTraces(Base_label_ROIs, Save_individual_acq_dir);

%% Section_03c(I): [Video Output] -- Label Video to Review Masking --  
Specific_ROIs = [];
[multi_mask_list, mask_gaps_list] = IRCE_AutomatedMaskReview(Specific_ROIs, Base_label_ROIs);

if all(isnan(multi_mask_list))
    disp('No ROIs with multiple labels, check for no mask ROIs')
else
    Specific_ROIs = multi_mask_list;

    %IRCE_MaskVideo(Save_individual_acq_dir, Specific_ROIs, roi_mask_Manual_Crops, Base_label_ROIs);
    IRCE_MaskVideo_IRMoverlay(Save_individual_acq_dir, Specific_ROIs, Base_label_ROIs, filtered_IRM_data, channel_LUTs{Contact_Channel}, roi_corners)
end

if isempty(find(~ismember(mask_gaps_list, multi_mask_list),1))
    disp('No other ROIs have gaps.')
else
    Specific_ROIs = mask_gaps_list(~ismember(mask_gaps_list, multi_mask_list));

    %IRCE_MaskVideo(Save_individual_acq_dir, Specific_ROIs, roi_mask_Manual_Crops, Base_label_ROIs);
    IRCE_MaskVideo_IRMoverlay(Save_individual_acq_dir, Specific_ROIs, Base_label_ROIs, filtered_IRM_data, channel_LUTs{Contact_Channel}, roi_corners)
end

IRCE_AutomatedMaskReview([], Base_label_ROIs);

clear mask_gaps_list multi_mask_list ans

%% Section_03c(II): [As Needed] -- Edit -- Regenerate Automated Masking for Specific ROI Masks -- 
cd(Save_individual_acq_dir)
%---------------------------------------------------------%
% Run Section S02e again when done
%---------------------------------------------------------%
Specific_ROIs = [];

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%
if ~isempty(Specific_ROIs)
    % Redo the ninja slice for a particular ROI(s)
    Base_label_ROIs = IRCE_RedoMaskSeperation(Specific_ROIs, roi_corners, filtered_IRM_data, median_threshold(Contact_Channel), channel_LUTs{Contact_Channel}, Save_individual_acq_dir);    
    
    IRCE_MaskVideo_IRMoverlay(Save_individual_acq_dir, Specific_ROIs, Base_label_ROIs, filtered_IRM_data, channel_LUTs{Contact_Channel}, roi_corners)
end

IRCE_AutomatedMaskReview(Specific_ROIs, Base_label_ROIs);



%% Section_03c(II): [As Needed] -- Edit -- Mask Specific ROIs Using CellPose
% 2024102 WIP KLS

cd(Save_individual_acq_dir)
%---------------------------------------------------------%
% Run Section S02e again when done
%---------------------------------------------------------%
Specific_ROIs = [];

cp = cellpose(Model='D:\01_Matlab\01_KLS_ImpulseResponse_CellROI_Extractionf_IRCE\Support Dependencies\IRM_cellpose_Model');

if ~isempty(Specific_ROIs)
    %Redo the ninja slice for a particular ROI(s)
    %[roi_mask_Manual_Crops, Base_label_ROIs] = IRCE_RedoMaskActiveContour(Specific_ROIs, roi_corners, filtered_IRM_data, roi_mask_Manual_Crops, channel_LUTs{Contact_Channel}, Base_label_ROIs);
    Base_label_ROIs = IRCE_RedoMaskCellpose(Specific_ROIs, roi_corners, filtered_IRM_data, cp, Base_label_ROIs, channel_LUTs{Contact_Channel}, Save_individual_acq_dir);


    IRCE_MaskVideo_IRMoverlay(Save_individual_acq_dir, Specific_ROIs, Base_label_ROIs, filtered_IRM_data, channel_LUTs{Contact_Channel}, roi_corners)
end

IRCE_AutomatedMaskReview([], Base_label_ROIs);



%% Section_03c(II): [As Needed] -- Edit -- Manually Draw Specific ROI Masks -- 
cd(Save_individual_acq_dir)
Specific_ROIs = [];

if ~isempty(Specific_ROIs)
    Base_label_ROIs = IRCE_HandDrawMask(Specific_ROIs, roi_corners, filtered_IRM_data, channel_LUTs{Contact_Channel}, Base_label_ROIs, Save_individual_acq_dir);

    IRCE_MaskVideo_IRMoverlay(Save_individual_acq_dir, Specific_ROIs, Base_label_ROIs, filtered_IRM_data, channel_LUTs{Contact_Channel}, roi_corners)
end

IRCE_AutomatedMaskReview(Specific_ROIs, Base_label_ROIs);



%% Section_03c(III): -- Manually Connect Any Remaining Disconnected Masks --
%---------------------------------------------------------%
% Manually connect cell rois
%---------------------------------------------------------%
% KLS_generateBrackets(n) <-- This command generates a blank set of brackets

% List out each Label # that represents the target cell
% i.e. [1 2] if the target cell is both label mask 1 and 2
% i.e. [4] if the target cell is label mask 4
% [] is equivalent to [1]
labelsForTheTargetCell = {
    []; ...
    []; ...
    []; ...
    []; ...
    []; ...
        []; ...
        []; ...
        []; ...
        []; ...
        []; ... %010 
    []; ...
    []; ...
};


%---------------------------------------------------------%
% Combine ROIS that are one cell
%---------------------------------------------------------%
i = 1;
while i <= length(labelsForTheTargetCell)
    % Check if this has section has already already been run before
    %final_ROI = length(labelsForTheTargetCell);
    %uniqueVals = unique(Base_label_ROIs{final_ROI,1}(:));
    %expectedVals = [0; final_ROI];

    uniqueVals = unique(Base_label_ROIs{i,1}(:));
    expectedVals = [0; i];
    if isequal(sort(uniqueVals), sort(expectedVals))
        disp(['ROI ' num2str(i) ' label assignment appears complete.']);

        %---------------------------------------------------------%
        % Dilate the final label to yeild the dilated mask
        %---------------------------------------------------------%
        SE_3 = strel('disk', 6); % Flat Structuring Element for Image dilation, disk of size 9
        
        Dilated_label_ROIs = Base_label_ROIs;
        
        ii = 1;
        while ii <= size(Base_label_ROIs{1},3)
            Dilated_label_ROIs{i,1}(:,:,ii) = imdilate(Base_label_ROIs{i,1}(:,:,ii),SE_3);
            ii = ii+1;
        end

        clear j labelsOnTheTargetCell targetLabel currentSet idx_target_dilatedLabel idx_target_baseLabel SE_3 n ii expectedVals uniqueVals final_ROI

        i = i + 1;
        continue;
    end

    
    % Extract the current set of labels
    currentSet = labelsForTheTargetCell{i};
    if isempty(currentSet)
        currentSet(1) = 1;
    end

    % If the only values in the mask are
    current_valid_labels = [0 i 999];
    current_labels = unique(Base_label_ROIs{i,1}(:));
    n = length(current_labels);
    switch n
        case 2 % Check if the current ROI has only labels 0 and i
            if sum(ismember(current_labels, current_valid_labels(1:2))) == 2
                i = i+1;
                continue
            end
        case 3 % Check if the current ROI has only labels 0, i and 999
            if sum(ismember(current_labels, current_valid_labels)) == 3
                i = i+1;
                continue
            end
    end

    % Iterate through each label in the current set, starting from the second label
    idx_target_baseLabel = (Base_label_ROIs{i,1} == currentSet(1));
    
    for j = 2:length(currentSet)
        % combined idx of the target label
        idx_target_baseLabel = bitor(Base_label_ROIs{i,1} == currentSet(j),idx_target_baseLabel);
    end
    % for all the pixels not associated with the target cells, set them to
        % 999
    Base_label_ROIs{i,1}(bitand(~idx_target_baseLabel, Base_label_ROIs{i,1}>0)) = 999;
    
    % for all the pixels that are associated with the target cells, set them to
        % ROI#, i
    Base_label_ROIs{i,1}(idx_target_baseLabel) = i;
    
    i = i + 1;
end

%---------------------------------------------------------%
% Dilate the final label to yeild the dilated mask
%---------------------------------------------------------%
SE_3 = strel('disk', 6); % Flat Structuring Element for Image dilation, disk of size 9

Dilated_label_ROIs = Base_label_ROIs;

for  i = 1:length(Base_label_ROIs)
    ii = 1;
    while ii <= size(Base_label_ROIs{1},3)
        Dilated_label_ROIs{i,1}(:,:,ii) = imdilate(Base_label_ROIs{i,1}(:,:,ii),SE_3);
        ii = ii+1;
    end
end

clear i j labelsOnTheTargetCell targetLabel currentSet idx_target_dilatedLabel idx_target_baseLabel SE_3 n ii expectedVals uniqueVals labelsForTheTargetCell


%% Section_03c(IV): [Video Output] -- Review Mask Assignments --    
Specific_ROIs = [];
IRCE_FinalMaskVideo_IRMoverlay(Save_individual_acq_dir, Specific_ROIs, Base_label_ROIs, filtered_IRM_data, channel_LUTs{Contact_Channel}, roi_corners);

clear multi_mask_list

%% Section_O3d: ROI channel data -- Gen channel data for each ROI --  
% Pull out data for each channel corrisponding to the Cell ROIs
base_data_ROIs = cell([num_ch 1]); % a cell array {n_channels by 1}
processed_data_ROIs = cell([num_ch 1]); % a cell array {n_channels by 1}

for i = 1:num_ch
    base_data_ROIs{i,1} = IRCE_ROIchannelcrops(Dilated_label_ROIs, base_data{i,1}, roi_corners);
end

for i = 1:num_ch
    processed_data_ROIs{i,1} = IRCE_ROIchannelcrops(Dilated_label_ROIs, processed_data{i,1}, roi_corners);
end

[Base_label_ROIs, ~] = IRCE_MaskAndROICorrnersCrop(Dilated_label_ROIs, Base_label_ROIs, roi_corners);
[Dilated_label_ROIs, roi_corners] = IRCE_MaskAndROICorrnersCrop(Dilated_label_ROIs, Dilated_label_ROIs, roi_corners);

cd(Save_individual_acq_dir)
save('Base_label_ROIs.mat', 'Base_label_ROIs', '-v7.3')
save('roi_corners.mat', 'roi_corners', '-v7.3')

clear i

%% Section_04a: [Optional] IRM Deformation Annotation --    

[roi_IRM_interface] = IRCE_IRMDeformationAnnotation(Save_individual_acq_dir, roi_corners, filtered_IRM_data, Base_label_ROIs, channel_LUTs{Contact_Channel}, Dilated_label_ROIs);

%% Section_04b: [As Needed] -- IRM Deformation Mask Automated --    
% WIP as of 20240716 - KLS
[roi_IRM_interface] = IRCE_IRMDeformationAnnotation_Automated(Save_individual_acq_dir, roi_corners, filtered_IRM_data, Base_label_ROIs, channel_LUTs{Contact_Channel}, Dilated_label_ROIs,median_threshold(Contact_Channel));

%% Section_04b: [As Needed] -- Redo Specific ROI IRM Deformation Annotation -- 
IRM_label_ROIs = cell([size(roi_corners,1) 1]);
ROIs_to_edit_flag = [];

%---------------------------------------------------------%
% Maunally Generate Masks -- Line through deformations
%---------------------------------------------------------%
cd(Save_individual_acq_dir)
if exist('roi_mask_IRM_Deformation.mat','file')
    load('roi_mask_IRM_Deformation','-mat')
    close all
    
    %---------------------------------------------------------%
    % edit ROIs
    %---------------------------------------------------------%
    num_possible_deformations_per_frame = 10;
    for n = 1:length(ROIs_to_edit_flag) % Loop over manually selected ROIs
        i = (ROIs_to_edit_flag(n));
        
        % size of roi_IRM_interface = [# ROIs, # frames, # of deformations]
        % Clear out old data for roi_IRM_interface
        for x = i
            for y = 1:size(filtered_IRM_data,3)
                for z = 1:num_possible_deformations_per_frame
                    roi_IRM_interface{x,y,z} = [];
                end
            end
        end
    
        clc
        disp(['***** Current ROI = ' num2str(i,'%03.f') ' *****'])
        %---------------------------------------------------------%n

        % loop over each cell ROI
        %---------------------------------------------------------%
        %Grab current ROI pixels: x and y
        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);
        img = filtered_IRM_data(min(y):max(y),min(x):max(x),:);

        %---------------------------------------------------------%
        % Loop over every frame
        %---------------------------------------------------------%
        ii = 1;
        while ii <= size(filtered_IRM_data,3) % Loop over time
            % slip the current frame if the current cell ROI is absence
            if isempty(find(Base_label_ROIs{i,1}(:,:,ii) == i,1))
                ii = ii+1;
                continue
            end

            %close all

            %figure()
            imshow(img(:,:,ii),channel_LUTs{Contact_Channel}, 'Border', 'tight','InitialMagnification',400)

            % Works with multi ROIs in one Image
            % Add all cell boundary to ch2 and ch3 (impulse and response) 
            [B,~] = bwboundaries(Dilated_label_ROIs{i,1}(:,:,ii) == i,'noholes');
            hold on
            for k = 1:length(B)
               boundary = B{k};
               plot(boundary(:,2), boundary(:,1),'Color',"#00FFFF",'LineStyle',':','LineWidth',2)
            end
            hold off

            g = gcf;
            g.WindowState = 'maximized';

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
close all

clc


%% Section_05a: [Optional] [TIF Output] -- Entire FOV (Raw and Processed) Data --   
clc

%---------------------------------------------------------%
% Save the seperate raw data from the channels
%---------------------------------------------------------%
for i = 1:num_ch
    file_name = ['raw_' channel_labels{i,1}]; % <--- Change Me as needed =
    KLS_save_double2tif(base_data{i,1}, file_name, Raw_dir);
end

%---------------------------------------------------------%
% Save the seperate processed data from the channels
%---------------------------------------------------------%
for i = 1:num_ch
    file_name = ['processed_' channel_labels{i,1}]; % <--- Change Me as needed = 
    KLS_save_double2tif(processed_data{i,1}, file_name, Processed_dir);
end

clear i file_name file_path


%% Section_05b: [Optional] [TIF Output] -- Save Cells ROI (Raw, Processed) Data  
%---------------------------------------------------------%
% Save individual Cell ROIs
%---------------------------------------------------------%
clc 

Specific_ROIs = [];

if isempty(Specific_ROIs)
    Specific_ROIs = 1:size(Dilated_label_ROIs,1);
end

for n = 1:length(Specific_ROIs) % Loop over manually selected ROIs
    ROI_n = (Specific_ROIs(n));
    
    folderName = ['Cell_' num2str(ROI_n,'%03.f')];
    folder_path = fullfile(ROIs_dir, folderName);
    if ~exist(folder_path, 'dir')
        % If folder doesn't exist, create it
        mkdir(folder_path);
    end

    %---------------------------------------------------------%
    % Save Labeled Data as tif files in the cell ROIs
    %---------------------------------------------------------%    
    file_name = ['Cell_' num2str(ROI_n,'%03.f') '_processed_label_dilated']; % <--- Change Me as needed =
    KLS_save_double2tif(Dilated_label_ROIs{ROI_n,1}, file_name, folder_path);

    file_name = ['Cell_' num2str(ROI_n,'%03.f') '_processed_label_base']; % <--- Change Me as needed =
    KLS_save_double2tif(Base_label_ROIs{ROI_n,1}, file_name, folder_path);

    %---------------------------------------------------------%
    % Save Processed Data as tif files in the cell ROIs
    %---------------------------------------------------------%    
    for i = 1:num_ch
        file_name = ['Cell_' num2str(ROI_n,'%03.f') '_processed_' channel_labels{i,1}]; % <--- Change Me as needed =
        KLS_save_double2tif(processed_data_ROIs{i,1}{ROI_n,1}, file_name, folder_path);
    end

    %---------------------------------------------------------%
    % If channels have different frequencies resize the less frequent data
    % for purpose of matching the frequency of timestamps
    %---------------------------------------------------------%
    largest_ch_idx = find(channel_freq == min(channel_freq));
    maxT = size(processed_data_ROIs{largest_ch_idx(1),1}{1,1},3);   
        
    for i = 1:num_ch
        % Check if the current channel is the same duration as the most
        % frequent channel
        if maxT ~= size(processed_data_ROIs{i,1}{ROI_n,1},3)
            Resized_mask = KLS_resizeMatrix(Dilated_label_ROIs{ROI_n,1} == ROI_n, maxT);

            Resized_ch = KLS_resizeMatrix(processed_data_ROIs{i,1}{ROI_n,1}, maxT);

            file_name = ['Cell_' num2str(ROI_n,'%03.f') '_processed_forTracking_' channel_labels{i,1}]; % <--- Change Me as needed =
            KLS_save_double2tif(Resized_ch.*Resized_mask, file_name, folder_path);                   
        else
            % Double check if the cell mask is the same duration as the
            % current channel
            if size(Dilated_label_ROIs{ROI_n,1},3) == maxT
                file_name = ['Cell_' num2str(ROI_n,'%03.f') '_processed_forTracking_' channel_labels{i,1}]; % <--- Change Me as needed =
                KLS_save_double2tif(processed_data_ROIs{i,1}{ROI_n,1} .* (Dilated_label_ROIs{ROI_n,1} == ROI_n), file_name, folder_path);
            else
                Resized_mask = KLS_resizeMatrix(Dilated_label_ROIs{ROI_n,1} == ROI_n, maxT);

                file_name = ['Cell_' num2str(ROI_n,'%03.f') '_processed_forTracking_' channel_labels{i,1}]; % <--- Change Me as needed =
                KLS_save_double2tif(processed_data_ROIs{i,1}{ROI_n,1} .* Resized_mask, file_name, folder_path);  
            end
        end
    end
end
 
close all


pause(1)

clear i file_name file_path Specific_ROIs ROI_n folderName folder_path largest_ch_idx maxT Resized_mask Resized_ch n

%% Section_06a(i): [Optional] -- Aggregate files for batch TrackMate Tracking 
clc 
%file_name_pattern = '*_processed_LysoBrite Blue.tif';
%file_name_pattern = '*_processed_forTracking_LysoTracker DND-99.tif';

file_name_pattern = '*_processed_forTracking_Binding tirf 647.tif';
file_name_pattern = '*_processed_forTracking_LysoTracker DND-99.tif';


file_name_pattern = '*_processed_forTracking_FOLR1 Binding.tif';
file_name_pattern = '*_processed_forTracking_LysoBrite Blue.tif';

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%

% Define the root directory to search in
rootDir = ROIs_dir; % Change this to your selected directory

% Define the target directory where files are temporarily saved
%---------------------------------------------------------%
% Move cd to data parent data save folder
%---------------------------------------------------------%
folderName = 'Aggregated';
aggregatedDir = fullfile(Save_individual_acq_dir, folderName);
if ~exist(aggregatedDir, 'dir')
    % If folder doesn't exist, create it
    mkdir(aggregatedDir);
end

% Search for all folders in the root directory
folders = dir(rootDir);
folders = folders([folders.isdir]); % Filter out anything that's not a directory

% Initialize a container to store original file paths
originalPaths = containers.Map;

% Loop through each folder to find and copy files
for i = 1:length(folders)
    % Skip the '.' and '..' folders
    if strcmp(folders(i).name, '.') || strcmp(folders(i).name, '..')
        continue;
    end
    
    % Construct the folder path
    folderPath = fullfile(rootDir, folders(i).name);
    
    % Search for files with a specific pattern in the current folder
    files = dir(fullfile(folderPath, file_name_pattern));
    
    % Copy each found file to the target directory and store its original path
    for j = 1:length(files)
        srcFile = fullfile(folderPath, files(j).name);
        dstFile = fullfile(aggregatedDir, files(j).name);
        copyfile(srcFile, dstFile);
        
        % Store the original path with the file name as the key
        originalPaths(files(j).name) = folderPath;
    end
end

disp('Files copied successfully.');


clear folders i folderPath j srcFile dstFile copyfile file_name_pattern




%% Section_06a(ii): [Optional] -- Add empty trackfiles | Rename | Move TrackMate Tracking files
% Need to have run section_06a(i)
clc

% Define file_name_pattern and file_name_pattern_new as cell arrays
file_name_pattern = {
    '*_processed_forTracking_Binding tirf 647-all-spots.csv',  ...
    '*_processed_forTracking_LysoBrite Blue-all-spots.csv', ...
    '*_processed_forTracking_t647 Short Expo. Binding-all-spots.csv', ...
    '*_processed_forTracking_LysoTracker DND-99-all-spots.csv', ...
    '*_processed_forTracking_FOLR1 Binding-all-spots.csv', ...
};

file_name_pattern_new = {
    '_All_Spots_Binding.csv', ...
    '_All_Spots_Response.csv', ...
    '_All_Spots_Binding.csv', ...
    '_All_Spots_Response.csv', ...
    '_All_Spots_Binding.csv', ...
};

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%

% Define the template CSV file path
templateCSVPath = 'G:\010_SC_Data\Cell_001_emptyTrackMate.csv';

%---------------------------------------------------------%
% Rename Files
%---------------------------------------------------------%

% Check if the template CSV file exists
if ~exist(templateCSVPath, 'file')
    error('Template CSV file does not exist: %s', templateCSVPath);
end

% Get a list of all .tif files in the folder
tifFiles = dir(fullfile(aggregatedDir, '*.tif'));

% Loop over each .tif file
for k = 1:length(tifFiles)
    % Get the .tif file name
    tifFileName = tifFiles(k).name;

    % Exclude files that contain '-all-spots' in the name, i.e., the CSV files
    if contains(tifFileName, '-all-spots')
        continue; % Skip this file
    end

    % Get the base name without extension
    [~, baseName, ~] = fileparts(tifFileName);

    % Construct the corresponding .csv file name
    csvFileName = [baseName, '-all-spots.csv'];
    csvFilePath = fullfile(aggregatedDir, csvFileName);

    % Check if the .csv file exists
    if ~exist(csvFilePath, 'file')
        % Copy the template CSV file and rename it
        copyfile(templateCSVPath, csvFilePath);
        fprintf('Copied template CSV to: %s\n', csvFilePath);
    end
end

% Loop over each pattern in the cell arrays
for p = 1:length(file_name_pattern)
    % Get the current patterns
    current_pattern = file_name_pattern{p};
    current_new = file_name_pattern_new{p};

    % Get a list of all CSV files in the folder matching the current pattern
    fileList = dir(fullfile(aggregatedDir, current_pattern));

    % Loop through each file
    for i = 1:length(fileList)
        % Get the old filename
        oldFileName = fileList(i).name;

        % Extract the cell number from the old filename (e.g., '001', '002')
        cellNumber = regexp(oldFileName, 'Cell_(\d+)_', 'tokens', 'once');

        % Skip if cell number is not found
        if isempty(cellNumber)
            warning('No cell number found in file: %s', oldFileName);
            continue;
        end

        % Construct the new filename
        newFileName = ['Cell_', cellNumber{1}, current_new];

        % Rename the file
        movefile(fullfile(aggregatedDir, oldFileName), fullfile(aggregatedDir, newFileName));
    end
end

clear fileList cellNumber newFileName i p current_pattern current_new
clc

%---------------------------------------------------------%
% Move the renamed files
%---------------------------------------------------------%

% Get the list of renamed files in the Aggregated folder
renamedFiles = dir(fullfile(aggregatedDir, '*.csv')); % Assuming the renamed files have .csv extension

% Loop through each renamed file and move it back to the correct folder
for i = 1:length(renamedFiles)
    % Get the renamed file name
    renamedFileName = renamedFiles(i).name;

    % Extract the number from the file name (e.g., '001' from 'Cell_001_All_Spots_Response.csv')
    cellNumber = regexp(renamedFileName, 'Cell_(\d+)_', 'tokens', 'once');
    if isempty(cellNumber)
        warning('No cell number found in file: %s', renamedFileName);
        continue;
    end

    % Search for the folder with the matching cell number in the rootDir
    folderName = ['Cell_', cellNumber{1}];  % Assuming the folder name matches the cell number exactly

    matchingFolder = fullfile(rootDir, folderName);

    % Check if the folder exists
    if isfolder(matchingFolder)
        % Move the renamed file to the original folder
        srcFile = fullfile(aggregatedDir, renamedFileName);
        dstFile = fullfile(matchingFolder, renamedFileName);
        movefile(srcFile, dstFile);

        % Display confirmation
        fprintf('Moved: %s -> %s\n', renamedFileName, dstFile);
    else
        warning('No matching folder found for cell number: %s', cellNumber{1});
    end
end

disp('Renamed files moved back to their respective folders successfully.');

%---------------------------------------------------------%
% Delete all files in the Aggregated folder
%---------------------------------------------------------%

% Confirm deletion with the user
deleteConfirmation = input('Do you want to delete all files in the Aggregated folder? (y/n): ', 's');
if strcmpi(deleteConfirmation, 'y')
    % Get a list of all files and directories in the Aggregated folder
    allItems = dir(fullfile(aggregatedDir, '*'));
    % Exclude '.' and '..'
    allItems = allItems(~ismember({allItems.name}, {'.', '..'}));

    % Loop through and delete each file or directory
    for i = 1:length(allItems)
        % Construct full file path
        itemPath = fullfile(aggregatedDir, allItems(i).name);

        % Check if it's a directory
        if allItems(i).isdir
            % Remove directory and its contents
            rmdir(itemPath, 's');
            fprintf('Deleted directory: %s\n', itemPath);
        else
            % Delete the file
            delete(itemPath);
            fprintf('Deleted file: %s\n', itemPath);
        end
    end

    disp('All files in the Aggregated folder have been deleted.');
else
    disp('Deletion canceled. No files were deleted.');
end

clear renamedFiles renamedFileName folderName matchingFolder srcFile dstFile i cellNumber
clear aggregatedDir file_name_pattern file_name_pattern_new file_name k files 
clear originalPaths templateCSVPath tifFiles


%% Section_06b: [Optional] Tracking Import

largest_ch_idx = find(channel_freq == min(channel_freq));
maxT = size(processed_data_ROIs{largest_ch_idx(1),1}{1},3);

% Get a list of all subfolders in the ROIs directory
subFolders = dir(ROIs_dir);
subFolders = subFolders([subFolders.isdir] & ~ismember({subFolders.name}, {'.', '..'}));

% Initialize a cell array to store the results
Impulse_Tracking_ROIs = cell([str2double(regexp(subFolders(end).name, '\d+','match')) 1]);
Response_Tracking_ROIs = cell([str2double(regexp(subFolders(end).name, '\d+','match')) 1]);

% Loop through each subfolder
for k = 1:length(subFolders)
    subFolderName = subFolders(k).name;
    subFolderPath = fullfile(ROIs_dir, subFolderName);
    
    folderNumber = str2double(regexp(subFolderName, '\d+', 'match')); % What ROI number is this?
    
    % Look for files containing 'All_Spots_Binding' in their name
    bindingFile = dir(fullfile(subFolderPath, '*All_Spots_Binding*.csv'));
    
    if ~isempty(bindingFile)
        % Define the full path of the found CSV file
        csvFilePath = fullfile(subFolderPath, bindingFile.name);
        
        % Read the matrix from the CSV file
        All_Spots = readmatrix(csvFilePath);
        All_Spots = All_Spots(2:end,:);

        % Perform the tracking operation
        STLN_Tracks = KLS_TMCSV_2_Tracks(All_Spots, maxT);
        
        % Store the results in the cell array
        Impulse_Tracking_ROIs{folderNumber, 1}.All_Spots = All_Spots;
        Impulse_Tracking_ROIs{folderNumber, 1}.STLN_Tracks = STLN_Tracks;
    end
    
    % Look for files containing 'All_Spots_Response' in their name
    responseFile = dir(fullfile(subFolderPath, '*All_Spots_Response*.csv'));
    
    if ~isempty(responseFile)
        % Define the full path of the found CSV file
        csvFilePath = fullfile(subFolderPath, responseFile.name);
        
        % Read the matrix from the CSV file
        All_Spots = readmatrix(csvFilePath);
        All_Spots = All_Spots(2:end,:);

        % Perform the tracking operation
        STLN_Tracks = KLS_TMCSV_2_Tracks(All_Spots, maxT);
        
        % Store the results in the cell array
        Response_Tracking_ROIs{folderNumber, 1}.All_Spots = All_Spots;
        Response_Tracking_ROIs{folderNumber, 1}.STLN_Tracks = STLN_Tracks;
    end
end

clear All_Spots STLN_Tracks subFolders k subFolderName subFolderPath csvFilePath folderNumber
clear bindingFile responseFile maxT largest_ch_idx

%% Section_06c: [Optional] [.mat Output] -- Save Cells ROI (Raw, Processed and Tracks under cell) Data  
%---------------------------------------------------------%
% Save individual Cell ROIs
%---------------------------------------------------------%
clc 

Specific_ROIs = [];

if isempty(Specific_ROIs)
    Specific_ROIs = 1:size(Dilated_label_ROIs,1);
end

for n = 1:length(Specific_ROIs) % Loop over manually selected ROIs
    ROI_n = (Specific_ROIs(n));
    
    folderName = ['Cell_' num2str(ROI_n,'%03.f')];
    folder_path = fullfile(ROIs_dir, folderName);
    if ~exist(folder_path, 'dir')
        % If folder doesn't exist, create it
        mkdir(folder_path);
    end

    %---------------------------------------------------------%
    % Save Impulse Tracks/Localizations that are under the cell
    %---------------------------------------------------------%
    if exist('Impulse_Tracking_ROIs','var')
        if isfield(Impulse_Tracking_ROIs{ROI_n,1},'STLN_Tracks')

            if length(unique(channel_freq)) > 1
                largest_ch_idx = find(channel_freq == min(channel_freq));
                maxT = size(processed_data_ROIs{largest_ch_idx(1),1}{1,1},3);
            else
                maxT = size(processed_data_ROIs{1,1}{1,1},3);
            end
    
            binary_img = KLS_resizeMatrix(processed_data_ROIs{Contact_Channel,1}{ROI_n,1}, maxT);
              
            tracks = Impulse_Tracking_ROIs{ROI_n}.STLN_Tracks;

            tracks = KLS_Tracks_in_Mask(binary_img, tracks);
                
            Tracks_in_mask = tracks;

            var_name = ['Cell_' num2str(ROI_n,'%03.f') '_Binding_STLN_Tracks' '.mat'];
            file_path = fullfile(folder_path, var_name);
            save(file_path,'tracks', '-v7.3')
        end
    end

    %---------------------------------------------------------%
    % Save Response Tracks/Localizations that are under the cell
    %---------------------------------------------------------%
    if exist('Response_Tracking_ROIs','var')
        if isfield(Response_Tracking_ROIs{ROI_n,1},'STLN_Tracks')

            if length(unique(channel_freq)) > 1
                largest_ch_idx = find(channel_freq == min(channel_freq));
                maxT = size(processed_data_ROIs{largest_ch_idx(1),1}{1,1},3);
            else
                maxT = size(processed_data_ROIs{1,1}{1,1},3);
            end
    
            binary_img = KLS_resizeMatrix(processed_data_ROIs{Contact_Channel,1}{ROI_n,1}, maxT);
              
            tracks = Response_Tracking_ROIs{ROI_n}.STLN_Tracks;

            tracks = KLS_Tracks_in_Mask(binary_img, tracks);
                
            Tracks_in_mask = tracks;
            var_name = ['Cell_' num2str(ROI_n,'%03.f') '_Response_STLN_Tracks' '.mat'];
            file_path = fullfile(folder_path, var_name);
            save(file_path,'tracks', '-v7.3')
        end
    end

end
 
close all

pause(1)

clear n ROI_n Specific_ROIs folderName folder_path binary_img tracks Tracks_in_mask
clear var_name largest_ch_idx maxT


%% Section_07a(i): [Optional] [Montage Output] Save Cell ROI Montages
%---------------------------------------------------------%
% Save individual Cell ROIs
%---------------------------------------------------------%
clc 

Specific_ROIs = [];
Reorder_channels = [1 3 2];
num_montage_frames = 11;

if isempty(Specific_ROIs)
    Specific_ROIs = 1:size(Dilated_label_ROIs,1);
end


for n = 1:length(Specific_ROIs) % Loop over manually selected ROIs
    ROI_n = (Specific_ROIs(n));
    
    folderName = ['Cell_' num2str(ROI_n,'%03.f')];
    folder_path = fullfile(ROIs_dir, folderName);
    if ~exist(folder_path, 'dir')
        % If folder doesn't exist, create it
        mkdir(folder_path);
    end

    %---------------------------------------------------------%
    % Save Montage Data for each cell ROI
    %---------------------------------------------------------%
    % n-channel data:
    curr_ROI_data = cell(size(processed_data_ROIs));
    for ii = 1:num_ch
        curr_ROI_data{ii,1} = processed_data_ROIs{ii,1}{ROI_n};
    end

    if ~isempty(Impulse_Tracking_ROIs{ROI_n,1})
        if exist('roi_IRM_interface','var')
            Montage_Name = fullfile(folder_path, ['Cell_' num2str(ROI_n,'%03.f') '_IRM_Impulse_Response_Montage_withTracks']); % <--- Change Me as needed = 
            
            IRCE_gen_Nch_montage(num_montage_frames, ...
                Time_stamps_address, ROI_n, curr_ROI_data, ...
                Dilated_label_ROIs{ROI_n,1}, channel_freq, channel_LUTs, ...
                channel_labels, channel_colors, roi_IRM_interface(ROI_n,:,:), ...
                Impulse_Tracking_ROIs{ROI_n,1}.STLN_Tracks, ...
                Impulse_track_annotation_channels, Montage_Name, Reorder_channels);
            
        else
            Montage_Name = fullfile(folder_path, ['Cell_' num2str(ROI_n,'%03.f') '_IRM_Impulse_Response_Montage_withTracks']); % <--- Change Me as needed = 
            
            IRCE_gen_Nch_montage(num_montage_frames, ...
                Time_stamps_address, ROI_n, curr_ROI_data, ...
                Dilated_label_ROIs{ROI_n,1}, channel_freq, channel_LUTs, ...
                channel_labels, channel_colors, [], ...
                Impulse_Tracking_ROIs{ROI_n,1}.STLN_Tracks, Impulse_track_annotation_channels, ...
                Montage_Name, Reorder_channels);
            
        end
    else
        if exist('roi_IRM_interface','var')

            Montage_Name = fullfile(folder_path, ['Cell_' num2str(ROI_n,'%03.f') '_IRM_Impulse_Response_Montage_noTracks']); % <--- Change Me as needed = 

            IRCE_gen_Nch_montage(num_montage_frames, ...
                Time_stamps_address, ROI_n, curr_ROI_data, ...
                Dilated_label_ROIs{ROI_n,1}, channel_freq, channel_LUTs, ...
                channel_labels, channel_colors, roi_IRM_interface(ROI_n,:,:), ...
                [], Impulse_track_annotation_channels, Montage_Name, Reorder_channels);
            
        else
            Montage_Name = fullfile(folder_path, ['Cell_' num2str(ROI_n,'%03.f') '_IRM_Impulse_Response_Montage_noTracks']); % <--- Change Me as needed = 
            
            IRCE_gen_Nch_montage(num_montage_frames, ...
                Time_stamps_address, ROI_n, curr_ROI_data, ...
                Dilated_label_ROIs{ROI_n,1}, channel_freq, channel_LUTs, ...
                channel_labels, channel_colors, [], ...
                [], Impulse_track_annotation_channels, Montage_Name, Reorder_channels);
        end
    end
end
 
close all


pause(1)

clear Specific_ROIs Reorder_channels num_montage_frames n ROI_n folderName
clear folder_path ii curr_ROI_data Montage_Name

%% Section_07a(ii): [Optional]  [Video Output] Save Cell ROI Videos
%---------------------------------------------------------%
% Save individual Cell ROIs
%---------------------------------------------------------%
clc 

Specific_ROIs = [];
Reorder_channels = [1 3 2];

if isempty(Specific_ROIs)
    Specific_ROIs = 1:size(Dilated_label_ROIs,1);
end


for n = 1:length(Specific_ROIs) % Loop over manually selected ROIs
    ROI_n = (Specific_ROIs(n));
    
    folderName = ['Cell_' num2str(ROI_n,'%03.f')];
    folder_path = fullfile(ROIs_dir, folderName);
    if ~exist(folder_path, 'dir')
        % If folder doesn't exist, create it
        mkdir(folder_path);
    end

    %---------------------------------------------------------%
    % Save Video data as mp4 files for each cell ROI
    %---------------------------------------------------------%
    % n-channel data:
    curr_ROI_data = cell(size(processed_data_ROIs));
    for ii = 1:num_ch
        curr_ROI_data{ii,1} = processed_data_ROIs{ii,1}{ROI_n};
    end

    if ~isempty(Impulse_Tracking_ROIs{ROI_n,1})
        if exist('roi_IRM_interface','var')
            Video_Name = fullfile(folder_path, ['Cell_' num2str(ROI_n,'%03.f') '_IRM_Impulse_Response_Video_withTracks']); % <--- Change Me as needed = 

            IRCE_gen_Nch_movie(Time_stamps_address, ROI_n, ...
                curr_ROI_data, Dilated_label_ROIs{ROI_n,1}, ...
                channel_freq, channel_LUTs, channel_labels, channel_colors, ...
                roi_IRM_interface(ROI_n,:,:), Impulse_Tracking_ROIs{ROI_n,1}.STLN_Tracks, ...
                Impulse_track_annotation_channels, Video_Name, Reorder_channels);
        else
            Video_Name = fullfile(folder_path, ['Cell_' num2str(ROI_n,'%03.f') '_IRM_Impulse_Response_Video_withTracks']); % <--- Change Me as needed =          

            IRCE_gen_Nch_movie(Time_stamps_address, ROI_n, ...
                curr_ROI_data, Dilated_label_ROIs{ROI_n,1}, ...
                channel_freq, channel_LUTs, channel_labels, channel_colors, ...
                [], Impulse_Tracking_ROIs{ROI_n,1}.STLN_Tracks, ...
                Impulse_track_annotation_channels, Video_Name, Reorder_channels);
        end
    else
        if exist('roi_IRM_interface','var')
            Video_Name = fullfile(folder_path, ['Cell_' num2str(ROI_n,'%03.f') '_IRM_Impulse_Response_Video_noTracks']); % <--- Change Me as needed = 

            IRCE_gen_Nch_movie(Time_stamps_address, ROI_n, ...
                curr_ROI_data, Dilated_label_ROIs{ROI_n,1}, ...
                channel_freq, channel_LUTs, channel_labels, channel_colors, ...
                roi_IRM_interface(ROI_n,:,:), [], ...
                Impulse_track_annotation_channels, Video_Name, Reorder_channels);
        else
            Video_Name = fullfile(folder_path, ['Cell_' num2str(ROI_n,'%03.f') '_IRM_Impulse_Response_Video_noTracks']); % <--- Change Me as needed =   

            IRCE_gen_Nch_movie(Time_stamps_address, ROI_n, ...
                curr_ROI_data, Dilated_label_ROIs{ROI_n,1}, ...
                channel_freq, channel_LUTs, channel_labels, channel_colors, ...
                [], [], ...
                Impulse_track_annotation_channels, Video_Name, Reorder_channels);  
        end
    end
end
 
close all


pause(1)

clear Specific_ROIs Reorder_channels n ROI_n folderName
clear folder_path ii curr_ROI_data Video_Name

%% Section_07c: [Optional] Aggrogate Videos and Montages Into Two Seperate Folders --  

rootDir = ROIs_dir;

folderName = 'Montages';
folder_path = fullfile(Save_individual_acq_dir, folderName);
if ~exist(folder_path, 'dir')
    % If folder doesn't exist, create it
    mkdir(folder_path);
end
targetDirMontage = folder_path;


folderName = 'Videos';
folder_path = fullfile(Save_individual_acq_dir, folderName);
if ~exist(folder_path, 'dir')
    % If folder doesn't exist, create it
    mkdir(folder_path);
end
targetDirVideo = folder_path;


IRCE_aggrogate_montage_and_Video(rootDir,targetDirMontage,targetDirVideo)

clear rootDir folderName folder_path targetDirMontage targetDirVideo

%% Section_08a: ROI Stats 
%---------------------------------------------------------%
% Do not change remaining code in this section, 
    % unless you want to add new stats or modify old ones
%---------------------------------------------------------%

%---------------------------------------------------------%
% Save individual Cell ROIs
%---------------------------------------------------------%
max_ROI = size(Base_label_ROIs, 1);

folderName = 'Stats';
folder_path = fullfile(Save_individual_acq_dir, folderName);
if ~exist(folder_path,'dir')
    % If folder doesn't exist, create it
    mkdir(folder_path);
end

file_path = fullfile(folder_path, 'Stats_ROIs.mat');
if isfile(file_path)
    disp('Old stats file found: loading...');
    load(file_path)
    disp('Done')
    have_old_stats_flag = 1;
else
    disp('No existing stats file found. Let''s make one shall we.')
    % Initialize Stats_ROIs cell array
    Stats_ROIs = cell(max_ROI, 1);

    have_old_stats_flag = 0;
end

% Pull in the recorded timepoints
Timing_seconds = IRCE_ND2_TimeStamps(Time_stamps_address);

ROI_n = 1;
while ROI_n <= max_ROI
    % Get current ROI data
    currentROI = Base_label_ROIs{ROI_n};
    [x, y, t] = size(currentROI);
    
    % Initialize stats structure
    stats = struct();
    stats.Centroid = nan(t, 2); % in px
    stats.Area = nan(t, 1); % in micron^2
    stats.Speed = nan(t, 1); % in micron/s
    stats.LandingIdx = nan; % in frame #
    stats.TotalDistance_micron = 0; % in micron
    stats.FinalDistance_micron = nan; % in micron
    stats.Sphericality = nan(t, 1); % in normalized units 0 to 1
    stats.Timing_sec = Timing_seconds; % in seconds
    stats.MajorAxisLength = nan(t, 1); % in normalized units 0 to 1

    folderName = ['Cell_' num2str(ROI_n,'%03.f')];
    folder_path = fullfile(ROIs_dir, folderName);
    stats.folder_path = folder_path; % Save the current folder to the stats file

    px_to_micron = 0.157; % px size in micron

    prevCentroid = nan(1, 2);
    initialCentroid = nan(1, 2);
    totalDistance = 0;
    
    for t_idx = 1:t
        mask = currentROI(:, :, t_idx) == ROI_n;
        if any(mask(:)) % If there is a mask at this time point
            props = regionprops(mask, 'Centroid', 'Area', 'Perimeter', 'Circularity','MajorAxisLength');
            
              
            
            if numel(props) > 1
                centroids = cat(1, props.Centroid);
                centroid = mean(centroids,1);
                stats.Centroid(t_idx, :) = centroid;
                
                area = sum([props.Area]) * (px_to_micron^2);
                stats.Area(t_idx) = area;
                
                
                [~, idx] = max([props.Area]);
                stats.Sphericality(t_idx) = props(idx).Circularity;

                [~, idx] = max([props.MajorAxisLength]);
                stats.MajorAxisLength(t_idx) = props(idx).MajorAxisLength .* px_to_micron;
            else
                centroid = props.Centroid;
                area = props(1).Area * (px_to_micron^2); % Convert area to micron squared  
                
                stats.Centroid(t_idx, :) = centroid;
                stats.Area(t_idx) = area;
                stats.Sphericality(t_idx) = props.Circularity;
                stats.MajorAxisLength(t_idx) = props.MajorAxisLength .* px_to_micron ;
            end

            
            if isnan(stats.LandingIdx)
                stats.LandingIdx = t_idx; % Note the first t index with a mask
                initialCentroid = centroid; % Set the initial centroid
            end
            
            if ~isnan(prevCentroid(1))
                % Calculate speed (distance/time)
                distance = sqrt(sum((centroid - prevCentroid).^2)) * px_to_micron; % in micron/frame
                totalDistance = totalDistance + distance;
                speed = distance / (stats.Timing_sec(t_idx)-stats.Timing_sec(t_idx-1)); % micron/s
                stats.Speed(t_idx) = speed;
            end
            
            prevCentroid = centroid;
        end
    end
    
    % Calculate final distance from initial position
    if ~isnan(initialCentroid(1))
        finalDistance = sqrt(sum((prevCentroid - initialCentroid).^2)) * px_to_micron;
        stats.FinalDistance_micron = finalDistance;
    end
    
    stats.TotalDistance_micron = totalDistance;
    
    % Save stats in the cell array
    if have_old_stats_flag == 1
        Stats_ROIs{ROI_n} = KLS_updateStruct(Stats_ROIs{ROI_n}, stats);
    else
        Stats_ROIs{ROI_n} = stats;
        
    end
    
    ROI_n = ROI_n + 1;
end
have_old_stats_flag = 1;

% Update Needed: If size(IRM,3) doesn't match the most frequent channel
% then the pulled time is incorrect
if exist('roi_IRM_interface','var') && ~isempty(roi_IRM_interface)
    t = size(roi_IRM_interface,2);

    ROI_n = 1;
    while ROI_n <= max_ROI
    
        %---------------------------------------------------------%
        %
        %---------------------------------------------------------%

        % Initialize variables for each 'n'
        found_first_nonempty = false;
        first_nonempty_value = nan;

        for j = 1:t
            nonempty_count = 0;
            
            % Check each cell in the third dimension
            for k = 1:size(roi_IRM_interface,3)
                if ~isempty(roi_IRM_interface{ROI_n, j, k})
                    nonempty_count = nonempty_count + 1;
                    
                    if ~found_first_nonempty
                        first_nonempty_value = j;

                        found_first_nonempty = true;
                    end
                end
            end
        end
        
        % Calculate the timing difference if first_nonempty_value is found
        % Update Here
        if ~isnan(first_nonempty_value)
            Stats_ROIs{ROI_n}.FirstDeformationIdx = first_nonempty_value;
            Stats_ROIs{ROI_n}.Time2Deformation_sec = Stats_ROIs{ROI_n}.Timing_sec(first_nonempty_value) - Stats_ROIs{ROI_n}.Timing_sec(Stats_ROIs{ROI_n}.LandingIdx);
        else
            Stats_ROIs{ROI_n}.FirstDeformation_sec = [];
            Stats_ROIs{ROI_n}.Time2Deformation_sec = [];
        end
        
        %---------------------------------------------------------%
        % Run some code to do rudimentary connection of IRM annotations
        %---------------------------------------------------------%

        % Initialize lists to store connecting and non-connecting data points and their time points
        connecting_chains = {};
        non_connecting_points = {};
        non_connecting_times = [];
    
        % Loop through the selected 'n' dimension
        for j = 1:t
            for k = 1:size(roi_IRM_interface,3)
                line1 = roi_IRM_interface{ROI_n, j, k};
                if isempty(line1)
                    continue;
                end
                is_connected = false;
    
                % Check for existing chains to extend
                for idx = 1:length(connecting_chains)
                    chain = connecting_chains{idx};
    
                    if j == chain{end, 2} + 1 && IRCE_check_overlap(line1, chain{end, 1})
                        connecting_chains{idx} = [connecting_chains{idx}; {line1, j, k}];
                        is_connected = true;
                        break;
                    end
                end
    
                % Start a new chain if no connection found
                if ~is_connected
                    connecting_chains{end+1} = {line1, j, k};
                end
            end
        end
    
        % Separate chains into connecting and non-connecting, filtering out single-point chains
        final_connecting_chains = {};
        for idx = 1:length(connecting_chains)
            chain = connecting_chains{idx};
            if size(chain, 1) > 1
                final_connecting_chains{end+1} = chain;
            else
                non_connecting_points{end+1} = chain;
            end
        end

        Stats_ROIs{ROI_n}.IRMAnnotations_Connected = final_connecting_chains;
        Stats_ROIs{ROI_n}.IRMAnnotations_Isolated = non_connecting_points;

        ROI_n = ROI_n + 1;
    end
end

sigma_threshold = 1.5;
if ~isempty(Response_Channel) && Response_Channel <= num_ch 
    Stats_ROIs = IRCE_Process_data(Base_label_ROIs, processed_data_ROIs{Response_Channel,1}, sigma_threshold, Stats_ROIs);

    ROI_n = 1;
    while ROI_n <= max_ROI
        % Initialize stats structure
        stats = struct();
    
        %---------------------------------------------------------%
        % Save Response Tracks/Localizations that are under the cell
        %---------------------------------------------------------%
        if exist('Response_Tracking_ROIs','var')
            if isfield(Response_Tracking_ROIs{ROI_n,1},'STLN_Tracks')
    
                if length(unique(channel_freq)) > 1
                    largest_ch_idx = find(channel_freq == min(channel_freq));
                    maxT = size(processed_data_ROIs{largest_ch_idx(1),1}{1,1},3);
                else
                    maxT = size(processed_data_ROIs{1,1}{1,1},3);
                end
        
                binary_img = KLS_resizeMatrix(processed_data_ROIs{Contact_Channel,1}{ROI_n,1}, maxT);
                  
                tracks = Response_Tracking_ROIs{ROI_n}.STLN_Tracks;
    
                tracks = KLS_Tracks_in_Mask(binary_img, tracks);
                    
                stats.Response_Tracking.All_Spots = Response_Tracking_ROIs{ROI_n}.All_Spots;
                stats.Response_Tracking.STLN_Tracks = tracks;
                stats.Response_Tracking.threshold_above_bkgd = sigma_threshold;
            end
        end
    
        % Save stats in the cell array
        if have_old_stats_flag == 1
            Stats_ROIs{ROI_n} = KLS_updateStruct(Stats_ROIs{ROI_n}, stats);
        else
            Stats_ROIs{ROI_n} = stats;
        end

        ROI_n = ROI_n + 1;
    end

    if exist('Response_Tracking_ROIs','var') && isfield(Stats_ROIs{1},'Response_Tracking')
        Stats_ROIs = IRCE_Process_Response_data_withTracks(Base_label_ROIs, processed_data_ROIs{Response_Channel,1}, sigma_threshold, Stats_ROIs, Response_Tracking_ROIs);
    end
end

sigma_threshold = 3;
if ~isempty(Impulse_Channel) && Impulse_Channel <= num_ch 
    Stats_ROIs = IRCE_Process_data_Impulse(Base_label_ROIs, processed_data_ROIs{Impulse_Channel,1}, sigma_threshold, Stats_ROIs);

    ROI_n = 1;
    while ROI_n <= max_ROI
        
        % Initialize stats structure
        stats = struct();
    
        %---------------------------------------------------------%
        % Save Impulse Tracks/Localizations that are under the cell
        %---------------------------------------------------------%
        if exist('Impulse_Tracking_ROIs','var')
            if isfield(Impulse_Tracking_ROIs{ROI_n,1},'STLN_Tracks')
    
                if length(unique(channel_freq)) > 1
                    largest_ch_idx = find(channel_freq == min(channel_freq));
                    maxT = size(processed_data_ROIs{largest_ch_idx(1),1}{1,1},3);
                else
                    maxT = size(processed_data_ROIs{1,1}{1,1},3);
                end
        
                binary_img = KLS_resizeMatrix(processed_data_ROIs{Contact_Channel,1}{ROI_n,1}, maxT);
                  
                tracks = Impulse_Tracking_ROIs{ROI_n}.STLN_Tracks;
    
                tracks = KLS_Tracks_in_Mask(binary_img, tracks);
    
                stats.Impulse_Tracking.All_Spots = Impulse_Tracking_ROIs{ROI_n}.All_Spots;
                stats.Impulse_Tracking.STLN_Tracks = tracks;
                stats.Impulse_Tracking.threshold_above_bkgd = sigma_threshold;
            end
        end
    
        % Save stats in the cell array
        if have_old_stats_flag == 1
            Stats_ROIs{ROI_n} = KLS_updateStruct(Stats_ROIs{ROI_n}, stats);
        else
            Stats_ROIs{ROI_n} = stats;
        end
        
        ROI_n = ROI_n + 1;
    end

    if exist('Impulse_Tracking_ROIs','var') && isfield(Stats_ROIs{1},'STLN_Tracks')
        Stats_ROIs = IRCE_Process_Impulse_data_withTracks(Base_label_ROIs, processed_data_ROIs{Impulse_Channel,1}, sigma_threshold, Stats_ROIs, Impulse_Tracking_ROIs);
    end
end

save(file_path,'Stats_ROIs','-v7.3')

close all


clear max_ROI folderName folder_path Timing_seconds ROI_n sigma_threshold currentROI
clear x y t stats t_idx mask px_to_micron prevCentroid initialCentroid totalDistance
clear props centroids area speed finalDistance found_first_nonempty first_nonempty_value
clear j connecting_chains non_connecting_points non_connecting_times k line1
clear is_connected idx chain final_connecting_chains distance centroid

%% Section_08b: [Optional] [.mat Output] -- Save to Cells ROI (individual stats files, and timing variable) Data  
%---------------------------------------------------------%
% Save individual Cell ROIs
%---------------------------------------------------------%
clc 

Specific_ROIs = [];

if isempty(Specific_ROIs)
    Specific_ROIs = 1:size(Dilated_label_ROIs,1);
end

%---------------------------------------------------------%
% Save tracks to the main Stats file
%---------------------------------------------------------%
max_ROI = size(Base_label_ROIs, 1);

folderName = 'Stats';
folder_path = fullfile(Save_individual_acq_dir, folderName);
if ~exist(folder_path,'dir')
    % If folder doesn't exist, create it
    mkdir(folder_path);
end

file_path = fullfile(folder_path, 'Stats_ROIs.mat');
if isfile(file_path)
    disp('Old stats file found: loading...');
    load(file_path)
    disp('Done')
    have_old_stats_flag = 1;
else
    disp('No existing stats file found. Let''s make one shall we.')
    % Initialize Stats_ROIs cell array
    Stats_ROIs = cell(max_ROI, 1);

    have_old_stats_flag = 0;
end

close all

for n = 1:length(Specific_ROIs) % Loop over manually selected ROIs
    ROI_n = (Specific_ROIs(n));
    
    folderName = ['Cell_' num2str(ROI_n,'%03.f')];
    folder_path = fullfile(ROIs_dir, folderName);
    if ~exist(folder_path, 'dir')
        % If folder doesn't exist, create it
        mkdir(folder_path);
    end


    ROI_stats = Stats_ROIs{ROI_n};

    ROI_stats.folder_path = folder_path; % Save the current folder to the stats file

    var_name = ['Cell_' num2str(ROI_n,'%03.f') '_Stats' '.mat'];
    file_path = fullfile(folder_path, var_name);
    save(file_path,'ROI_stats', '-v7.3')
end
 
close all

clear max_ROI folderName folder_path ROI_n currentROI x y t stats largest_ch_idx
clear maxT binary_img tracks var_name file_path n Specific_ROIs

%% Section_09a: [Optional] General Figures for ROI stats
close all
rolling_avg = []; % movmean for the Response data, if you have it. Run once with an empty [].
max_ROI = size(Base_label_ROIs, 1);

folderName = 'Stats';
folder_path = fullfile(Save_individual_acq_dir, folderName);
if ~exist(folder_path,'dir')
    % If folder doesn't exist, create it
    mkdir(folder_path);
end

% check if IRM is not the most frequent channel;
if channel_freq(Contact_Channel) ~= 1
    % if IRM is the not most freqent data, then count 1:IRM_freq:max_frames
    Ch1_Time_idx = 1:channel_freq(Contact_Channel):size(processed_data_ROIs{find(channel_freq == min(channel_freq),1),1}{1,1},3);
else
    % if IRM is the most freqent data, then count from 1 to # IRM frames
    Ch1_Time_idx = 1:channel_freq(Contact_Channel):size(processed_data_ROIs{Contact_Channel,1}{1,1},3);
end

%---------------------------------------------------------%
% Individual ROI area as it changes in time
%---------------------------------------------------------%
figure;
fig = gcf;
fig.Position = Pos + [0 0 0 0];
hold on;

for i = 1:max_ROI
    stats = Stats_ROIs{i};
    landingIdx = stats.LandingIdx;
    
    if ~isnan(landingIdx)
        % Align the time points
        IRM_Timing = stats.Timing_sec(Ch1_Time_idx);
        
        alignedTime = (IRM_Timing - IRM_Timing(landingIdx)) / 60; % Convert to minutes
        areaData = stats.Area;
        
        % Plot each ROI's area

        plot(alignedTime(linspace(1,size(alignedTime,2), size(areaData,1))), areaData, 'DisplayName', ['Cell ' num2str(i, '%03.f')]);
    end
end

xlabel('Time (minutes)');
ylabel('Area (\mum^2)');
title('Area of Each ROI Over Time');
legend off;
hold off;

file_path = fullfile(folder_path, 'ROI area over time.fig');
savefig(file_path)
file_path = fullfile(folder_path, 'ROI area over time.eps');
saveas(gcf, file_path)
file_path = fullfile(folder_path, 'ROI area over time.png');
saveas(gcf,file_path)

%---------------------------------------------------------%
% Mean ROI area in time
%---------------------------------------------------------%
if max_ROI > 1
    % Determine the maximum number of time points across all ROIs
    maxTimePoints = max(cellfun(@(s) length(s.Area), Stats_ROIs));

    % Initialize arrays for mean and SEM calculation
    allAreas = nan(max_ROI, maxTimePoints);
    alignedTimes = nan(max_ROI, maxTimePoints);

    for i = 1:max_ROI
        stats = Stats_ROIs{i};
        landingIdx = stats.LandingIdx;

        if ~isnan(landingIdx)
            % Align the time points
            areaData = stats.Area(landingIdx:end);

            % Fill the allAreas and alignedTimes arrays
            allAreas(i, 1:length(areaData)) = areaData;
        end
    end

    % Calculate mean and SEM
    meanArea = mean(allAreas, 1,'omitnan');
    %semArea = std(allAreas,1,'omitnan') ./ sqrt(sum(~isnan(allAreas), 1));
    semArea = SEM_calc(allAreas,0.05); % Calculate 95% CI using SEM, assumes guassian population
    semArea(semArea > max(allAreas,[],'all')) = nan;

    % Find common time points for plotting
    commonTimePoints = IRM_Timing ./ 60;

    % Plot mean area with SEM
    figure;
    fig = gcf;
    fig.Position = Pos + [Pos(3) 0 0 0];
    hold on;
    commonTimePoints = commonTimePoints - commonTimePoints(1);
    t = commonTimePoints(linspace(1,size(commonTimePoints,2),size(meanArea,2)));
    plot(t, meanArea, 'k', 'LineWidth', 2, 'DisplayName', 'Mean Area');

    % Plot SEM as shaded area'
    non_nan_idx = ~isnan(meanArea) & ~isnan(semArea);

    fill([t(non_nan_idx), fliplr(t(non_nan_idx))], ...
         [meanArea(non_nan_idx) + semArea(non_nan_idx), fliplr(meanArea(non_nan_idx) - semArea(non_nan_idx))], ...
         'cyan','FaceAlpha',0.8);

    xlabel('Time (minutes)');
    ylabel('Mean Area (\mum^2)');
    title('Mean Area Over Time with SEM');
    legend off;
    hold off;

    file_path = fullfile(folder_path, 'Mean ROI area over time.fig');
    savefig(file_path)
    file_path = fullfile(folder_path,'Mean ROI area over time.eps');
    saveas(gcf,file_path)
    file_path = fullfile(folder_path,'Mean ROI area over time.png');
    saveas(gcf,file_path)    
end

%---------------------------------------------------------%
% Displacement in time, ploted as 0,0 start point
%---------------------------------------------------------%
figure;
fig = gcf;
fig.Position = Pos + [Pos(3)*2 0 0 0];
hold on;

px_to_micron = 0.157; % px to micron conversion factor

for i = 1:max_ROI
    stats = Stats_ROIs{i};
    
    if ~isnan(stats.LandingIdx)
        % Extract centroid positions
        xPos = stats.Centroid(:, 1) * px_to_micron; % Convert to microns
        yPos = stats.Centroid(:, 2) * px_to_micron; % Convert to microns
        
        % Remove NaN values
        validIndices = ~isnan(xPos) & ~isnan(yPos);
        xPos = xPos(validIndices);
        yPos = yPos(validIndices);
        
        % Normalize positions by subtracting the initial position
        if ~isempty(xPos) && ~isempty(yPos)
            xPos = xPos - xPos(1);
            yPos = yPos - yPos(1);
        
            % Plot the normalized centroid positions
            plot(xPos, yPos, 'DisplayName', ['Cell ' num2str(i, '%03.f')]);
        end
    end
end

xlabel('X Position (\mum)');
ylabel('Y Position (\mum)');
title('Normalized Centroid Positions Over Time');
legend off;
hold off;

file_path = fullfile(folder_path, 'ROI Displacement over time.fig');
savefig(file_path)
file_path = fullfile(folder_path, 'ROI Displacement over time.eps');
saveas(gcf, file_path)
file_path = fullfile(folder_path, 'ROI Displacement over time.png');
saveas(gcf, file_path)

%---------------------------------------------------------%
% Distance traveled
%---------------------------------------------------------%
% Extract data from Stats_ROIs
totalDistances = nan(max_ROI, 1);
finalDistances = nan(max_ROI, 1);
maxCellWidth = nan(max_ROI, 1);
speeds = [];

for i = 1:max_ROI
    stats = Stats_ROIs{i};
    
    if ~isnan(stats.TotalDistance_micron)
        totalDistances(i) = stats.TotalDistance_micron;
    end
    
    if ~isnan(stats.FinalDistance_micron)
        finalDistances(i) = stats.FinalDistance_micron;
    end
    
    speeds = [speeds; stats.Speed(~isnan(stats.Speed))];

    maxCellWidth(i) = max(stats.MajorAxisLength);
end

% Create the combined plot
figure;
fig = gcf;
fig.Position = Pos + [Pos(3)*3 0 0 0];
hold on;

% Scatter plots
scatter(ones(size(totalDistances)), totalDistances, 'filled', 'DisplayName', 'Total Distance Traveled');
scatter(2 * ones(size(finalDistances)), finalDistances, 'filled', 'DisplayName', 'Distance from Initial Point');

% Box-and-whisker plots
boxplot(totalDistances, 'positions', 1, 'Colors', 'k', 'Widths', 0.3);
boxplot(finalDistances, 'positions', 2, 'Colors', 'k', 'Widths', 0.3);

% Set plot labels and title
xticks([1 2]);
xticklabels({'Cumulative Distance', 'Final Displacement'});
ylabel('Distance (µm)');
legend off;
box off

hold off;
file_path = fullfile(folder_path, 'ROI Distance Traveled.fig');
savefig(file_path)
file_path = fullfile(folder_path, 'ROI Distance Traveled.eps');
saveas(gcf, file_path)
file_path = fullfile(folder_path, 'ROI Distance Traveled.png');
saveas(gcf, file_path)

%---------------------------------------------------------%
% Distance traveled (Final Disp./Cellwidth)
%---------------------------------------------------------%
% Just Final Displacement Relative to Max Cell Width
% Create the combined plot
figure;
fig = gcf;
fig.Position = Pos + [Pos(3)*3 -Pos(4)-90 0 0];
hold on;

% Scatter plots
plotSpread(finalDistances ./ maxCellWidth,[],[],[], 4)
%scatter(1 * ones(size(finalDistances)), finalDistances ./ maxCellWidth, 'filled', 'DisplayName', 'Distance from Initial Point');

% Box-and-whisker plots
boxplot(finalDistances ./ maxCellWidth, 'positions', 1, 'Colors', 'k', 'Widths', 0.3);

% Set plot labels and title
xticks([1]);
xticklabels({'Final Displacement/Max Cell Width'});
ylabel('Distance as a ratio of cell width (µm/µm)');
legend off;
box off

hold off;
file_path = fullfile(folder_path, 'ROI Distance Traveled as ratio of cell width.fig');
savefig(file_path)
file_path = fullfile(folder_path, 'ROI Distance Traveled as ratio of cell width.eps');
saveas(gcf, file_path)
file_path = fullfile(folder_path, 'ROI Distance Traveled as ratio of cell width.png');
saveas(gcf, file_path)

if isfield(Stats_ROIs{1,1},'integrated_response_above_background_guass')
    %---------------------------------------------------------%
    % Mean Integrated Lyso in time
    %---------------------------------------------------------%
    figure()
    fig = gcf;
    fig.Position = Pos + [Pos(3)*0 -Pos(4)-90 0 0];


    if length(Stats_ROIs{1,1}.Area) < length(Stats_ROIs{1,1}.integrated_response_above_background_guass)
        idx_map_to_time_idx = 1:channel_freq(Contact_Channel):length(Stats_ROIs{1,1}.integrated_response_above_background_guass);

        IRCE_Plotting_mean_ResponseSignal_sorted(Stats_ROIs, "#DB3069", [], idx_map_to_time_idx, rolling_avg)
    else
        IRCE_Plotting_mean_ResponseSignal_sorted(Stats_ROIs, "#DB3069", [], [], rolling_avg)
    end

    xlabel('Time (minutes)');
    ylabel('Integrated Lyso (Photons)');
    title('Mean Lyso Response Over Time');

    if ~isempty(rolling_avg)
        file_path = fullfile(folder_path, 'Mean Lyso Response (movmean).fig');
        savefig(file_path)
        file_path = fullfile(folder_path, 'Mean Lyso Response (movmean).eps');
        saveas(gcf, file_path)
        file_path = fullfile(folder_path, 'Mean Lyso Response (movmean).png');
        saveas(gcf, file_path)   
    else
        file_path = fullfile(folder_path, 'Mean Lyso Response.fig');
        savefig(file_path)
        file_path = fullfile(folder_path, 'Mean Lyso Response.eps');
        saveas(gcf, file_path)
        file_path = fullfile(folder_path, 'Mean Lyso Response.png');
        saveas(gcf, file_path)   
    end
 
end

if isfield(Stats_ROIs{1,1},'integrated_response_above_background_guass')
    %---------------------------------------------------------%
    % Single Cell Traces of Integrated Lyso in time
    %---------------------------------------------------------%
    figure()
    fig = gcf;
    fig.Position = Pos + [Pos(3)*1 -Pos(4)-90 0 0];

    if length(Stats_ROIs{1,1}.Area) < length(Stats_ROIs{1,1}.integrated_response_above_background_guass)
        idx_map_to_time_idx = 1:channel_freq(Contact_Channel):length(Stats_ROIs{1,1}.integrated_response_above_background_guass);

        IRCE_Plotting_SC_ResponseSignal_sorted(Stats_ROIs, [0 0 1], [], idx_map_to_time_idx)
    else
        IRCE_Plotting_SC_ResponseSignal_sorted(Stats_ROIs, [0 0 1], [], [])
    end

    xlabel('Time (minutes)');
    ylabel('Integrated Lyso (Photons)');
    title('Cell Traces: Lyso Response Over Time');

    
    file_path = fullfile(folder_path, 'Cell Traces - Lyso Response.fig');
    savefig(file_path)
    file_path = fullfile(folder_path, 'Cell Traces - Lyso Response.eps');
    saveas(gcf, file_path)
    file_path = fullfile(folder_path, 'Cell Traces - Lyso Response.png');
    saveas(gcf, file_path)    
 end

if isfield(Stats_ROIs{1,1},'integrated_response_above_background_guass_withlocalization')
    %---------------------------------------------------------%
    % Mean Integrated Lyso in time
    %---------------------------------------------------------%
    figure()
    fig = gcf;
    fig.Position = Pos + [Pos(3)*0 -Pos(4)-90-150 0 0];

    if length(Stats_ROIs{1,1}.Area) < length(Stats_ROIs{1,1}.integrated_response_above_background_guass_withlocalization)
        idx_map_to_time_idx = 1:channel_freq(Response_Channel):length(Stats_ROIs{1,1}.Timing_sec);

        IRCE_Plotting_mean_ResponseSignal_sorted_withLocalizationFilter(Stats_ROIs, "#DB3069", [], idx_map_to_time_idx)
    else
        IRCE_Plotting_mean_ResponseSignal_sorted_withLocalizationFilter(Stats_ROIs, "#DB3069", [], [])
    end

    xlabel('Time (minutes)');
    ylabel('Integrated Lyso (Photons)');
    title({'Mean Lyso Response Over Time', 'localization filtered'});
 
    file_path = fullfile(folder_path, 'Mean Lyso Response - localization filtered.fig');
    savefig(file_path)
    file_path = fullfile(folder_path, 'Mean Lyso Response - localization filtered.eps');
    saveas(gcf, file_path)
    file_path = fullfile(folder_path, 'Mean Lyso Response - localization filtered.png');
    saveas(gcf, file_path)
end

if isfield(Stats_ROIs{1,1},'integrated_response_above_background_guass_withlocalization')
    %---------------------------------------------------------%
    % Single Cell Traces of Integrated Lyso in time
    %---------------------------------------------------------%
    figure()
    fig = gcf;
    fig.Position = Pos + [Pos(3)*1 -Pos(4)-90-150 0 0];

    if length(Stats_ROIs{1,1}.Area) < length(Stats_ROIs{1,1}.integrated_response_above_background_guass_withlocalization)
        idx_map_to_time_idx = 1:channel_freq(Response_Channel):length(Stats_ROIs{1,1}.Timing_sec);

        IRCE_Plotting_SC_ResponseSignal_sorted_withLocalizationFilter(Stats_ROIs, [0 0 1], [], idx_map_to_time_idx)
    else
        IRCE_Plotting_SC_ResponseSignal_sorted_withLocalizationFilter(Stats_ROIs, [0 0 1], [], [])
    end

    xlabel('Time (minutes)');
    ylabel('Integrated Lyso (Photons)');
    title({'Cell Traces: Lyso Response Over Time', 'localization filtered'});

    file_path = fullfile(folder_path, 'Cell Traces - Lyso Response - localization filtered.fig');
    savefig(file_path)
    file_path = fullfile(folder_path, 'Cell Traces - Lyso Response - localization filtered.eps');
    saveas(gcf, file_path)
    file_path = fullfile(folder_path, 'Cell Traces - Lyso Response - localization filtered.png');
    saveas(gcf, file_path)
 end

if isfield(Stats_ROIs{1,1},'integrated_impulse_above_background')
    %---------------------------------------------------------%
    % Mean Integrated Lyso in time
    %---------------------------------------------------------%
    figure()
    fig = gcf;
    fig.Position = Pos + [Pos(3)*2 -Pos(4)-90 0 0];

    if length(Stats_ROIs{1,1}.Area) < length(Stats_ROIs{1,1}.integrated_impulse_above_background)
        idx_map_to_time_idx = 1:channel_freq(Impulse_Channel):length(Stats_ROIs{1,1}.Timing_sec);

        IRCE_Plotting_mean_ImpulseSignal_sorted(Stats_ROIs, "#DB3069", [], idx_map_to_time_idx)
    else
        IRCE_Plotting_mean_ImpulseSignal_sorted(Stats_ROIs, "#DB3069", [], [])
    end


    xlabel('Time (minutes)');
    ylabel('Integrated Binding Signal (Photons)');
    title('Mean Binding Signal Over Time');


    file_path = fullfile(folder_path, 'Mean Binding Signal.fig');
    savefig(file_path)
    file_path = fullfile(folder_path, 'Mean Binding Signal.eps');
    saveas(gcf, file_path)
    file_path = fullfile(folder_path, 'Mean Binding Signal.png');
    saveas(gcf, file_path)
end

if isfield(Stats_ROIs{1,1},'integrated_impulse_above_background')
    %---------------------------------------------------------%
    % Single Cell Traces of Integrated Lyso in time
    %---------------------------------------------------------%
    figure()
    fig = gcf;
    fig.Position = Pos + [Pos(3)*2 -Pos(4)-90-150 0 0];

    % if the impulse ch sampling freq > 1
    if length(Stats_ROIs{1,1}.Timing_sec) > length(Stats_ROIs{1,1}.integrated_impulse_above_background)
        idx_map_to_time_idx = 1:channel_freq(Impulse_Channel):length(Stats_ROIs{1,1}.Timing_sec);

        IRCE_Plotting_SC_ImpulseSignal_sorted(Stats_ROIs, [0 0 1], [], idx_map_to_time_idx)
    else
        IRCE_Plotting_SC_ImpulseSignal_sorted(Stats_ROIs, [0 0 1], [], [])
    end

    xlabel('Time (minutes)');
    ylabel('Integrated Binding Signal (Photons)');
    title('Cell Traces: Binding Signal Over Time');

    file_path = fullfile(folder_path, 'Cell Traces - Binding Impulse.fig');
    savefig(file_path)
    file_path = fullfile(folder_path, 'Cell Traces - Binding Impulse.eps');
    saveas(gcf, file_path)
    file_path = fullfile(folder_path, 'Cell Traces - Binding Impulse.png');
    saveas(gcf, file_path)
 end



clear rolling_avg folderName folder_path Ch1_Time_idx fig i stats landingIdx 
clear IRM_Timing alignedTime areaData maxTimePoints allAreas meanArea semArea 
clear commonTimePoints t non_nan_idx file_path px_to_micron xPos yPos validIndices 
clear totalDistances finalDistances maxCellWidth speeds idx_map_to_time_idx max_ROI
clear alignedTimes 

%% Section_09b: [Optional] 3D Spagetti Plot for 1 specific ROI --
% WIP needs to be updated for Nch formating
    % 20240919 working
clc
ROI_n = 2;

%---------------------------------------------------------%
% Some basic settings
%---------------------------------------------------------%
outline_width = 3;
spot_width = 1.5;
localization_marker_size = 75;
view_3D = [45 15];

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%

% n-channel Video:
curr_ROI_data = cell(size(processed_data_ROIs));
for ii = 1:num_ch
    curr_ROI_data{ii,1} = processed_data_ROIs{ii,1}{ROI_n};
end

% 3-D visualization of the above sections (spagetti plot)
close all

%----------------------------------------------------------
% Check all channels to make sure to resize all data to the max frequency
%----------------------------------------------------------
largest_ch_idx = find(channel_freq == min(channel_freq));
maxT = size(curr_ROI_data{largest_ch_idx(1),1},3);

resized_data = cell(size(curr_ROI_data));
for i = 1:num_ch
    resized_data{i,1} = KLS_resizeMatrix(curr_ROI_data{i,1}, maxT);
end

binary_img = KLS_resizeMatrix(Dilated_label_ROIs{ROI_n,1}, maxT);
photo_LL640 = KLS_resizeMatrix(resized_data{Impulse_track_annotation_channels, 1}, maxT);

disp(['Max Frame = ' num2str(size(photo_LL640,3))])

    frames_with_contacts = squeeze(sum(binary_img,[1 2])) > 0;
    landing_idx = find(frames_with_contacts, 1, 'first');
    last_contact_idx = find(frames_with_contacts, 1, 'last');
z_limits = [landing_idx last_contact_idx]; % <----- Make sure to set me

% Time_stamps_address from above ^^^
Timing_seconds = IRCE_ND2_TimeStamps(Time_stamps_address);

z = Timing_seconds/60; % make z axis in minutes
z = sort(z,'ascend');

Planes_to_add = round(linspace(z_limits(1), z_limits(2),4));
z = z-z(z_limits(1)); % for timing make first plane zero

   
if exist('Impulse_Tracking_ROIs','var')
    tracks = Impulse_Tracking_ROIs{ROI_n}.STLN_Tracks;
    tracks(tracks == 0) = nan; % remove zero, zero coords
    tracks = KLS_Fill_TrackZeros(tracks); % fill track gaps    
end

%---------------------------------------------------------%
% Get bounding boxes for each labeled region
%---------------------------------------------------------%
temp_label = Dilated_label_ROIs{ROI_n,1};
temp_label = double(temp_label == ROI_n);
temp_label(temp_label == 0) = nan;
min_label = min(temp_label,[],3);
min_label(isnan(min_label)) = 0;
stats = regionprops(min_label, 'BoundingBox','Centroid');

%---------------------------------------------------------%
% Filter Tracks for entirely within in the cell ROI
%---------------------------------------------------------%
%tracks = KLS_Tracks_in_Mask(binary_img, tracks); % WIP logic is off - KLS
    % 20241105, since tracks are already localizized with the cell mask applied,
    % then this is redundent

% For checking image size
img = squeeze(photo_LL640(:,:,1));

x1 = 1;
x2 = size(img,2); % x is columns
y1 = 1;
y2 = size(img,1); % y is rows

if exist('tracks','var')
    x = tracks(:,:,1);
    y = tracks(:,:,2);
    x(x > 0) = x(x > 0);
    y(y > 0) = y(y > 0);
end

figure()

%---------------------------------------------------------%
% Add each tracks in 3D with more than k links
%---------------------------------------------------------%
k = 1;

hold on
if exist('tracks','var')
    for i = 1:size(tracks,1)
        x_plot = x(i,:,1);
        y_plot = y(i,:,1);
        x_plot(x_plot == 0) = nan;
        y_plot(y_plot == 0) = nan;
        if size(find(x_plot>0),2) > k
            plot3(x_plot,y_plot,z(1:Impulse_Tracking_ROIs{ROI_n}.All_Spots(end,8)+1),'LineWidth',3)
        end
        C = colororder;
    end
end
hold off

%{
%---------------------------------------------------------%
% Project each track with more than k links on to the last frame
%---------------------------------------------------------%
% Last frame is defined by z(z_limits(2))
hold on
for i = 1:size(tracks,1)
    
    x_plot = x(i,:,1);
    y_plot = y(i,:,1);
    x_plot(x_plot == 0) = nan;
    y_plot(y_plot == 0) = nan;
    
    if size(find(x_plot>0),2) > k
        plot3(x_plot,y_plot,z(z_limits(2)).*ones(size(x)),'LineWidth',2,'Color','y')
    end
end
hold off
%}
%---------------------------------------------------------%
% Add several planes of imagings and their localizations
%---------------------------------------------------------%

hold on
for i =  Planes_to_add(1:end-1)
    if exist('Impulse_Tracking_ROIs','var')
        if i <= Impulse_Tracking_ROIs{ROI_n}.All_Spots(end,8)+1 % Tracks can end before the IRM data
            x_plot = x(:,i,1);
            y_plot = y(:,i,1);
            x_plot(x_plot == 0) = nan;
            y_plot(y_plot == 0) = nan;
            z_plot = ones(size(x_plot));
            scatter3(x_plot,y_plot,z(i).*ones(size(x_plot)),localization_marker_size,'MarkerEdgeColor',"#44FFD1",'Clipping','off','LineWidth',spot_width)
        end
    end
    %   boundary on top of the plot
    [row,col] = find(binary_img(:,:,i),1,'first');
    if ~isempty(row)
        boundary = bwtraceboundary(binary_img(:,:,i),[row, col],'E');
        plot3(boundary(:,2),boundary(:,1),z(i).*ones(size(boundary(:,2))),'Color',"#00FFFF",'LineStyle',':','LineWidth',outline_width);
    end
    
    [X,Y] = meshgrid(x1:x2,y1:y2);
    Z = z(i).*ones(size(X));
    img = squeeze(photo_LL640(:,:,i));
    wp_02 = warp(X,Y,Z,img,channel_LUTs{Impulse_track_annotation_channels});
    wp_02.FaceAlpha = 0.7;    
end
hold off

%---------------------------------------------------------%
% Do something special for the last plane
%---------------------------------------------------------%
% What is that special thing?
% Max project the binary to get the full cell area?
% Max projection of LL640 to see all?
    % For the moment I don't want to do anything special
hold on
for i =  Planes_to_add(end)
    if exist('Impulse_Tracking_ROIs','var')
        if i <= Impulse_Tracking_ROIs{ROI_n}.All_Spots(end,8)+1 % Tracks can end before the IRM data
            x_plot = x(:,i,1);
            y_plot = y(:,i,1);
            x_plot(x_plot == 0) = nan;
            y_plot(y_plot == 0) = nan;
            z_plot = ones(size(x_plot));
            scatter3(x_plot,y_plot,z(i).*ones(size(x_plot)),localization_marker_size,'MarkerEdgeColor',"#44FFD1",'Clipping','off','LineWidth',spot_width)
        end        
    end

    %{
    max_proj_binary = max(binary_img,[],3);
    [row,col] = find(max_proj_binary,1,'first');
    if ~isempty(row)
        boundary = bwtraceboundary(max_proj_binary,[row, col],'E');
        plot3(boundary(:,2),boundary(:,1),z(i).*ones(size(boundary(:,2))),'Color','y','LineStyle',':','LineWidth',outline_width);
    end
    %}
    %   boundary on top of the plot
    [row,col] = find(binary_img(:,:,i),1,'first');
    if ~isempty(row)
        boundary = bwtraceboundary(binary_img(:,:,i),[row, col],'E');
        plot3(boundary(:,2),boundary(:,1),z(i).*ones(size(boundary(:,2))),'Color',"#00FFFF",'LineStyle',':','LineWidth',outline_width);
    end
    
    [X,Y] = meshgrid(x1:x2,y1:y2);
    Z = z(i).*ones(size(X));
    hold on
%---------------------------------------------------------%
% Max projection
    %img = max(photo_LL640,[],3);
    %wp_02 = warp(X,Y,Z,img,[Impulse_LUT(1) Impulse_LUT(2)+1000]);
    img = squeeze(photo_LL640(:,:,i));
    wp_02 = warp(X,Y,Z,img,channel_LUTs{Impulse_track_annotation_channels});
    wp_02.FaceAlpha = 0.9;    
end
hold off


% ----------------------------------------------------------- %
% Formating
% Formating

box off
set(gca,'ZDir','reverse')
%set(gcf, 'color', 'none');    
%set(gca, 'color', 'none');

xlim([1 size(img,2)]);
xlabel('x (µm)')

scaling_factor = 0.157; % micron/px

% Calculate the length of the image in micrometers
img_length_um = size(img,2) * scaling_factor;

% Determine the xticks and xticklabels based on the length of the image
if img_length_um < 10
    xticks_vals = 5 / scaling_factor;
    xticklabels_vals = {'5'};
elseif img_length_um < 20
    xticks_vals = [5 / scaling_factor, 10 / scaling_factor, 15 / scaling_factor];
    xticklabels_vals = {'5', '10', '15'};
else
    % For lengths greater than or equal to 20 microns, use increments of 10 microns
    max_tick = floor(img_length_um / 10) * 10;
    xticks_vals = (10:10:max_tick) / scaling_factor;
    xticklabels_vals = arrayfun(@num2str, 10:10:max_tick, 'UniformOutput', false);
end

xticks(xticks_vals);
xticklabels(xticklabels_vals)


ylim([1 size(img,1)]);
ylabel('y (µm)')
% Calculate the length of the image in micrometers
img_width_um = size(img,1) * scaling_factor;

% Determine the xticks and xticklabels based on the length of the image
if img_width_um < 10
    yticks_vals = 5 / scaling_factor;
    yticklabels_vals = {'5'};
elseif img_width_um < 20
    yticks_vals = [5 / scaling_factor, 10 / scaling_factor, 15 / scaling_factor];
    yticklabels_vals = {'5', '10', '15'};
else
    % For lengths greater than or equal to 20 microns, use increments of 10 microns
    max_tick = floor(img_length_um / 10) * 10;
    yticks_vals = (10:10:max_tick) / scaling_factor;
    yticklabels_vals = arrayfun(@num2str, 10:10:max_tick, 'UniformOutput', false);
end

yticks(yticks_vals);
yticklabels(yticklabels_vals)

zlabel('Time (min)');
zlim([z(z_limits(1)) z(z_limits(2))])


%view(gca,[-45 38]);
view(gca,view_3D);

Pos = get(0,'defaultfigureposition');
fig = gcf;
fig.Position = Pos+[0 -750 Pos(3) Pos(4)*1.5];

%---------------------------------------------------------%
% Save figure in
%---------------------------------------------------------%
cd(ROIs_dir);

folderName = ['Cell_' num2str(ROI_n,'%03.f')];
if ~exist(fullfile(cd, folderName), 'dir')
    % If folder doesn't exist, create it
    mkdir(folder_path);
end
% Move to the folder
cd(folderName);
    
if exist('tracks','var')
    Tracks_in_mask = tracks;
    var_name = ['Cell_' num2str(ROI_n,'%03.f') '_STLN_Tracks' '.mat'];
    save(var_name,'tracks', '-v7.3')
end

set(gcf,'Color',[1 1 1])
%set(gcf, 'Position', [100 200 screen_size(3)*0.95 size(data,1)*2]);
set(gca,'Visible','on')
set(gcf,'InvertHardCopy','off');
set(gcf,'PaperPositionMode','auto')

% Save the displayed figure as PNG
saveas(gcf,['Cell_' num2str(ROI_n,'%03.f') '_SpagettiPlot' '.png'], 'png');
saveas(gcf,['Cell_' num2str(ROI_n,'%03.f') '_SpagettiPlot' '.fig'], 'fig');




%% Section_09c: [Optional] [Video Output] Save Cell ROI Videos, All Possible Annotations
%---------------------------------------------------------%
% Save individual Cell ROIs
%---------------------------------------------------------%
clc 

%channel_LUTs = {[2200 6200]; [0 250]; [0 500]; [0 2500];};

Specific_ROIs = [];
Reorder_channels = [1 3 2];

if isempty(Specific_ROIs)
    Specific_ROIs = 1:size(Dilated_label_ROIs,1);
end


for n = 1:length(Specific_ROIs) % Loop over manually selected ROIs
    ROI_n = (Specific_ROIs(n));
    
    folderName = ['Cell_' num2str(ROI_n,'%03.f')];
    folder_path = fullfile(ROIs_dir, folderName);
    if ~exist(folder_path, 'dir')
        % If folder doesn't exist, create it
        mkdir(folder_path);
    end

    % n-channel Video:
    curr_ROI_data = cell(size(processed_data_ROIs));
    for ii = 1:num_ch
        curr_ROI_data{ii,1} = processed_data_ROIs{ii,1}{ROI_n};
    end

    %---------------------------------------------------------%
    % Save Video data as mp4 files for each cell ROI
    % Three channel Video: IRM, CH2 - Impulse, and Ch3 - Response
    %  If you have both Ch2 and Ch3 data this runs 
    %---------------------------------------------------------%
    Video_Name = fullfile(folder_path, ['Cell_' num2str(ROI_n,'%03.f') '_IRM_Impulse_Response_Video_extraAnnotations']);

    Impulse_integrated_flag = 1;
    Response_integrated_flag = 1;

    if isempty(Impulse_Tracking_ROIs{ROI_n,1})
        current_Impulse_Tracks = [];
    else
        current_Impulse_Tracks = Impulse_Tracking_ROIs{ROI_n,1}.STLN_Tracks;
    end
    if isempty(Response_Tracking_ROIs{ROI_n,1})
        current_Response_Tracks = [];
    else
        current_Response_Tracks = Response_Tracking_ROIs{ROI_n, 1}.STLN_Tracks;
    end
    if exist('roi_IRM_interface(','var')
        current_IRM_interface = roi_IRM_interface(ROI_n,:,:);
    else
        current_IRM_interface = [];
    end
    %{
    IRCE_gen_Nch_movie_ExtraAnnotations(Time_stamps_address, ROI_n, ...
        curr_ROI_data, Dilated_label_ROIs{ROI_n,1}, ...
        channel_freq, channel_LUTs, channel_labels, channel_colors, ...
        current_IRM_interface, current_Impulse_Tracks, ...
        current_Response_Tracks, ...
        Impulse_track_annotation_channels, Response_track_annotation_channels, ...
        Video_Name, Reorder_channels, Stats_ROIs{ROI_n,1}, ...
        Impulse_integrated_flag, Response_integrated_flag)
    %}
    % Let not do any low-pH vesicle tracking
    IRCE_gen_Nch_movie_ExtraAnnotations(Time_stamps_address, ROI_n, ...
        curr_ROI_data, Dilated_label_ROIs{ROI_n,1}, ...
        channel_freq, channel_LUTs, channel_labels, channel_colors, ...
        current_IRM_interface, current_Impulse_Tracks, ...
        [], ...
        Impulse_track_annotation_channels, Response_track_annotation_channels, ...
        Video_Name, Reorder_channels, Stats_ROIs{ROI_n,1}, ...
        Impulse_integrated_flag, Response_integrated_flag)
end
 
close all


pause(1)

clear Specific_ROIs Reorder_channels n ROI_n folderName folder_path curr_ROI_data 
clear ii Video_Name Impulse_integrated_flag Reseponse_integrated_flag

%% Section_09d: [Optional] [Video Output] -- Cell ROI IRM Deformation Annotation v. Ch3_Response --    
% WIP as for 20240910 | KLS needs to be updated for Nch formating of
% channel data (i.e. IRM as the contact ch, and Lyso as the response ch,
% not specific variables)
%---------------------------------------------------------%
% Save individual Cell ROIs
%---------------------------------------------------------%
Specific_ROIs = [];

if isempty(Specific_ROIs)
    Specific_ROIs = 1:size(Dilated_label_ROIs,1);
end

for n = 1:length(Specific_ROIs) % Loop over manually selected ROIs
    
    ROI_n = (Specific_ROIs(n));

    
    folderName = ['Cell_' num2str(ROI_n,'%03.f')];
    folder_path = fullfile(ROIs_dir, folderName);
    if ~exist(folder_path, 'dir')
        % If folder doesn't exist, create it
        mkdir(folder_path);
    end

% If you need to reorder any channel use this
%{
    % Reorder LysoBlue so that it matach IRM
    original_matrix = Ch3_Response_ROIs{ROI_n,1};
    repeated_matrix = zeros(size(Ch1_corr_IRM_ROIs{ROI_n,1}));
    for ii = 1:10
        iii = 1;
        while iii <= 7
            repeated_matrix(:,:,(7*(ii-1))+iii) = original_matrix(:,:,ii);
            iii = iii +1;
        end
    end
%}
    temp_ch_labels.Ch1 = channel_labels{1};
    temp_ch_labels.Ch3 = channel_labels{3};
    Video_Name = fullfile(folder_path, ['Cell_' num2str(ROI_n,'%03.f') '_IRM_Deformation_Response_Video']); % <--- Change Me as needed = 
    IRCE_gen_IRMResponse_movie(Dilated_label_ROIs{ROI_n,1}, ...
        processed_data_ROIs{ROI_n,1}, channel_LUTs{1}, ...
        processed_data_ROIs{ROI_n,3}, channel_LUTs{3}, ...
        Time_stamps_address, Video_Name, ...
        ROI_n, roi_IRM_interface(ROI_n,:,:), temp_ch_labels); 
end

   

%%
%%
%%
%%
%%
%% [WIP] Date__S06(IIa): IRM Deformation Mapped to Binding
close all

%---------------------------------------------------------%
% Save figure in
%---------------------------------------------------------%
cd(ROIs_dir);

folderName = ['Cell_' num2str(ROI_n,'%03.f')];
if ~exist(fullfile(cd, folderName), 'dir')
    % If folder doesn't exist, create it
    mkdir(folderName);
end
% Move to the folder
cd(folderName);

% Pull in the recorded timepoints
fileID = fopen(Time_stamps_address, 'r');

% Assuming time format in the file is like 12:34.567 (mm:ss.ms)
% '%d' reads an integer, ':%d.%f' reads the seconds and milliseconds
data = fscanf(fileID, '%d:%d.%f', [3, Inf]);

% Close the file
fclose(fileID);

% Convert data to total seconds
switch floor(max(log10(data(3,:))))
    case -Inf
            % Convert data to total seconds
            % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
            Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 10;
    case 0
        % Convert data to total seconds
        % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
        Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 10;    
    case 1
        % Convert data to total seconds
        % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
        Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 100;                   
    case 2
        % Convert data to total seconds
        % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
        Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 1000;                   
    case 3
        % Convert data to total seconds
        % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
        Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 10000;                   
end
    
knocking_flag = 0;
patch_flag = 0;
remove_bkgd_IRM = 0;

ROI_n = 1;
max_z = 650;
min_z = -1450;
cell_mask = Dilated_label_ROIs{ROI_n,1};
ch1_data = Ch1_corr_IRM_ROIs{ROI_n,1};
ch2_data = Ch2_Impulse_ROIs{ROI_n,1};
ch3_data = Ch3_Response_ROIs{ROI_n,1};

IRM_bkgd_lvl = mean(ch1_data(~cell_mask),'all'); % figure out the local bkgd around the cell

ch1_data = ch1_data-IRM_bkgd_lvl; % scale things so bkgd == 0

%----------------------------------------------------------
% Check if Ch1 or Ch2 needs to be remapped to get paired
% data for each frame for each channel
%----------------------------------------------------------
maxT = max([size(ch1_data,3), size(ch2_data,3), size(ch3_data,3)]);
ch1_data = KLS_resizeMatrix(ch1_data, maxT);
ch2_data = KLS_resizeMatrix(ch2_data, maxT);
ch3_data = KLS_resizeMatrix(ch3_data, maxT);
cell_mask = KLS_resizeMatrix(cell_mask, maxT);
t_range = 1:800;
%t_range = 450:550;
y_range = 70:100;
x_range = 79:108;

t_range = 1:size(ch1_data,3);
y_range = 1:size(ch1_data,2);
x_range = 1:size(ch1_data,1);
ch1_data = ch1_data(x_range,y_range,t_range);
ch2_data = ch2_data(x_range,y_range,t_range);
ch3_data = ch3_data(x_range,y_range,t_range);
cell_mask = cell_mask(x_range,y_range,t_range);

% Generate video object
Video_Name = 'Test_IRM_warp';
vidObj = VideoWriter([Video_Name '.mp4'], 'MPEG-4');

vidObj.FrameRate = min(size(ch1_data,3) / 10, 120);
close(vidObj);
open(vidObj);
    
figure()
fig = gcf;
fig.Position = Pos+[0 -Pos(4)-90 Pos(3) Pos(4)];  
hold on

for i = 1:size(ch1_data,3)
    curr_mask = squeeze(cell_mask(:,:,i));
    curr_IRM = squeeze(ch1_data(:,:,i));
    IRM_bkgd = mean(curr_IRM(~curr_mask),'all');
    
    if remove_bkgd_IRM == 1
        curr_IRM(~curr_mask) = nan;
    end

    curr_wrap = squeeze(ch3_data(:,:,i));
    warp(curr_IRM,curr_wrap)
    %colormap(gcf,parula(500))
    %colormap(gcf,'jet')

   
    colormap(bipolar(201, 0.46,'linear'))
    
    xlim([1 size(ch1_data,2)])
    ylim([1 size(ch1_data,1)])
    zlim([min_z max_z])
    xticks('')
    yticks('')
    zticks('')
    %zlabel('IRM (AU)')
    box off
    caxis([450 7500])
    %caxis(Impulse_LUT)
    %view(gca,[-70 -50]);
    %title(['Frame = ' num2str(i)])
    text(1,1-25,max_z,[num2str(data(1, i)) 'min: ' num2str(data(2, i)) 's: ' num2str(data(3, i)) 'ms'],'horizontalalignment','center')
    set(gca, 'XDir','reverse')
    set(gca, 'YDir','reverse')
    
    hold on
    if patch_flag == 1
        
            %IRM floor
            patch([1 size(ch1_data,1) size(ch1_data,1) 1], [1 1 size(ch1_data,2) size(ch1_data,2)], [IRM_bkgd IRM_bkgd IRM_bkgd IRM_bkgd],'black') 

            %Walls above floor
            x = [1 size(ch1_data,1) size(ch1_data,1) 1];
            y = [1 1 size(ch1_data,2) size(ch1_data,2)];
            z = [IRM_bkgd IRM_bkgd IRM_bkgd IRM_bkgd];

            % Define the coordinates for the sides
            x_side1 = [1 1 1 1];
            y_side1 = [1 1 size(ch1_data,2) size(ch1_data,2)];
            z_side1 = [IRM_bkgd max_z max_z IRM_bkgd];

            x_side2 = [size(ch1_data,1) size(ch1_data,1) size(ch1_data,1) size(ch1_data,1)];
            y_side2 = [1 1 size(ch1_data,2) size(ch1_data,2)];
            z_side2 = [IRM_bkgd max_z max_z IRM_bkgd];

            x_side3 = [1 size(ch1_data,1) size(ch1_data,1) 1];
            y_side3 = [1 1 1 1];
            z_side3 = [IRM_bkgd IRM_bkgd max_z max_z];

            x_side4 = [1 size(ch1_data,1) size(ch1_data,1) 1];
            y_side4 = [size(ch1_data,2) size(ch1_data,2) size(ch1_data,2) size(ch1_data,2)];
            z_side4 = [IRM_bkgd IRM_bkgd max_z max_z];

            % Draw the side patches
            patch(x_side1, y_side1, z_side1, 'black');
            patch(x_side2, y_side2, z_side2, 'black');
            patch(x_side3, y_side3, z_side3, 'black');
            patch(x_side4, y_side4, z_side4, 'black');

            % Draw the sides of the box
            for ii = 1:4
                patch(x_sides(:,ii), y_sides(:,ii), z_sides1(:,ii),'black'); % Side patches
            end
    end
    hold off
    
    switch(knocking_flag)
        case 1
            for ii = 1:10
                view(gca,[0 60]);
                    currFrame = getframe(gcf);
                    writeVideo(vidObj, currFrame);
            end
            for ii = 60:90
                view(gca,[0 ii]);
                    currFrame = getframe(gcf);
                    writeVideo(vidObj, currFrame);
            end
            for ii = 1:10
                view(gca,[0 90]);
                    currFrame = getframe(gcf);
                    writeVideo(vidObj, currFrame);
            end
            for ii = 90:-1:60
                view(gca,[0 ii]);
                    currFrame = getframe(gcf);
                    writeVideo(vidObj, currFrame);
            end
            for ii = 1:10
                view(gca,[0 60]);
                    currFrame = getframe(gcf);
                    writeVideo(vidObj, currFrame);
            end            
        case 0
            for ii = 1:1
                view(gca,[-45 80]);
                    currFrame = getframe(gcf);
                    writeVideo(vidObj, currFrame);
            end 
            
            if i == 425
                annotation(gcf,'textarrow',[0.26 0.56],...
                    [0.62 0.50],'String',{'Rapid Drop in','Lyso Signal'});
            end
            if i == 480
                delete(findall(gcf,'type','annotation'))
            end
            if i == 490
                annotation(gcf,'textarrow',[0.36 0.55],...
                    [0.70 0.65],...
                    'String',{'Next IRM image:','Upward IRM Deformation'});
            end
            if i == 540
                delete(findall(gcf,'type','annotation'))
            end
    end
end
close(vidObj);
close all

%% [WIP] Date__S06(IIb): IRM Deformation Mapped to Binding -- Individual Images
% Generate video object
Video_Name = 'Test_IRM_only';
vidObj = VideoWriter([Video_Name '.mp4'], 'MPEG-4');
vidObj.FrameRate = 60;
close(vidObj);
open(vidObj);
    
figure()
fig = gcf;
fig.Position = Pos+[0 -Pos(4)-90 Pos(3) Pos(4)];  
hold on

for i = 1:size(ch1_data,3)
    curr_mask = squeeze(cell_mask(:,:,i));
    curr_IRM = squeeze(ch1_data(:,:,i));
    IRM_bkgd = mean(curr_IRM(~curr_mask),'all');


    %imshow(curr_IRM,channel_LUTs{Contact_Channel},'InitialMagnification',600)
    imshow(curr_IRM,[min_z max_z],'InitialMagnification',600)
            if any(i == 168:234)
                annotation(gcf,'textarrow',[0.313953488372093 0.5],...
                    [0.809027777777778 0.506944444444444],...
                    'String',{'Rapid Drop in','Lyso Signal'},'TextColor','white');
            end
            if any(i == 234:264)
                annotation(gcf,'textarrow',[0.387596899224806 0.445736434108527],...
                    [0.767361111111111 0.513888888888889],...
                    'String',{'Next IRM Image:','Upward IRM Deformation'},'TextColor','white');
            end
    currFrame = getframe(gcf);
    writeVideo(vidObj, currFrame);
end
hold off
close(vidObj);
close all
%
% Generate video object
Video_Name = 'Test_Lyso_only';
vidObj = VideoWriter([Video_Name '.mp4'], 'MPEG-4');
vidObj.FrameRate = 60;
close(vidObj);
open(vidObj);
    
figure()
fig = gcf;
fig.Position = Pos+[0 -Pos(4)-90 Pos(3) Pos(4)];  
hold on

for i = 1:size(ch1_data,3)
    curr_mask = squeeze(cell_mask(:,:,i));
    curr_Lyso = squeeze(ch3_data(:,:,i));
    IRM_bkgd = mean(curr_IRM(~curr_mask),'all');


    %imshow(curr_IRM,channel_LUTs{Contact_Channel},'InitialMagnification',600)
    imshow(curr_Lyso,[450 7500],'InitialMagnification',600)
    colormap(gcf,'jet')
            if any(i == 168:234)
                annotation(gcf,'textarrow',[0.313953488372093 0.5],...
                    [0.809027777777778 0.506944444444444],...
                    'String',{'Rapid Drop in','Lyso Signal'},'TextColor','white');
            end
            
            if any(i == 234:264)
                annotation(gcf,'textarrow',[0.387596899224806 0.445736434108527],...
                    [0.767361111111111 0.513888888888889],...
                    'String',{'Next IRM Image:','Upward IRM Deformation'},'TextColor','white');
            end
            
    currFrame = getframe(gcf);
    writeVideo(vidObj, currFrame);
end
hold off
close(vidObj);
close all

%% [WIP] Date__S06(IIc): IRM Deformation Mapped to Binding -- 3 subplots

% Define video object for combined video
Video_Name = 'Combined_Video';
vidObj = VideoWriter([Video_Name '.mp4'], 'MPEG-4');
vidObj.FrameRate = 120;
open(vidObj);

% Create figure for subplots
figure()
fig = gcf;
fig.Position = [100, 100, 1800, 600];  % Adjust the figure size as needed

y = 11:19;
x = ones(size(y)).*11;
for i = 1:size(ch1_data, 3)
    % Extract current frames
    curr_mask = squeeze(cell_mask(:, :, i));
    curr_IRM = squeeze(ch1_data(:, :, i));
    curr_Lyso = squeeze(ch3_data(:, :, i));

    % Calculate background
    IRM_bkgd = mean(curr_IRM(~curr_mask), 'all');

    % Create subplots
    subplot(1, 3, 1);
    imshow(curr_IRM, [min_z max_z], 'InitialMagnification', 600);
    colormap(gca,'grey')
    title('IRM');
    
    hold on;
    if any(i == 490)
        annotation('arrow',[0.14, 0.18], [0.54, 0.54], 'Color','yellow','HeadWidth',16,'LineWidth',8);
    end
    if i == 530
        delete(findall(gcf,'type','annotation'))
    end
    
    plot(x, y, 'b', 'LineWidth', 2); % Vertical red line
    hold off;

    subplot(1, 3, 2);
    imshow(curr_Lyso, [450 7500], 'InitialMagnification', 600);
    colormap(gca,'jet');
    title('Lyso');
    hold on;
    if any(i == 425)
        annotation('arrow',[0.43, 0.46], [0.52, 0.52], 'Color','yellow','HeadWidth',16,'LineWidth',8);
    end
    if i == 480
        delete(findall(gcf,'type','annotation'))
    end
    % Add a line to the second image for the intensity profile

   
    plot(x, y, 'r--', 'LineWidth', 2); % Vertical red line
    hold off;

    subplot(1, 3, 3);
    % Extract intensity values along the line in the second image
    lyso_profile = curr_Lyso(y, round(x(1)));
    IRM_proflie = curr_IRM(y, round(x(1))); 
    
    yyaxis left;
    plot((0:size(y,2)-1).*0.157,lyso_profile,'r--');
    title('Line Profile');
    xlabel('µm');
    ylabel('Lyso Intensity (AU)', 'Color', 'r');
    ax = gca;
    ax.YColor = 'r';
    ylim([450 12500])
    xlim([0.157 (size(y,2)-1).*0.157])
    
    yyaxis right;
    plot((0:size(y,2)-1).*0.157,IRM_proflie+2500,'b');
    ylabel('IRM (AU)', 'Color', 'b');
    ax = gca;
    ax.YColor = 'b';
    
    ylim([-1500+2500 800+2500])
    box off
    % Capture the frame for the video
    currFrame = getframe(gcf);
    writeVideo(vidObj, currFrame);
end

% Close the video object
close(vidObj);
close all;

%% [WIP - Working OK] Date__S06(IId): IRM Deformation (overlay blue/orange/magenta depending on Lyso dye) 
ROI_n = 1;

IRM_LUT_cutoffs = [-1500 1500]; % this is after making bkgd intenisty == 0;
Lyso_LUT_cutoffs = [0 150];

%---------------------------------------------------------%
% Save figure in
%---------------------------------------------------------%
cd(ROIs_dir);

folderName = ['Cell_' num2str(ROI_n,'%03.f')];
if ~exist(fullfile(cd, folderName), 'dir')
    % If folder doesn't exist, create it
    mkdir(folderName);
end
% Move to the folder
cd(folderName);

cell_mask = Dilated_label_ROIs{ROI_n,1};
ch1_data = Ch1_corr_IRM_ROIs{ROI_n,1};
ch2_data = Ch2_Impulse_ROIs{ROI_n,1};
ch3_data = Ch3_Response_ROIs{ROI_n,1};

% Pull in the recorded timepoints
fileID = fopen(Time_stamps_address, 'r');

% Assuming time format in the file is like 12:34.567 (mm:ss.ms)
% '%d' reads an integer, ':%d.%f' reads the seconds and milliseconds
data = fscanf(fileID, '%d:%d.%f', [3, Inf]);

switch floor(max(log10(data(3,:))))
    case -Inf
            % Convert data to total seconds
            % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
            Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 10;
    case 0
        % Convert data to total seconds
        % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
        Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 10;    
    case 1
        % Convert data to total seconds
        % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
        Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 100;                   
    case 2
        % Convert data to total seconds
        % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
        Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 1000;                   
    case 3
        % Convert data to total seconds
        % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
        Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 10000;                   
end

% Close the file
fclose(fileID);

% Convert data to total seconds
% data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
Tstamps = Timing_seconds/60; % convert s to min

size_vector = KLS_label2ROIsize(cell_mask == ROI_n, 0.157);
first_frame_landed = find(size_vector > 0,1,'first'); % first frame the cell landed
Tstamps_since_landed = Tstamps - Tstamps(first_frame_landed);  

IRM_bkgd_lvl = mean(ch1_data(~cell_mask),'all');

ch1_data = ch1_data-IRM_bkgd_lvl;
%----------------------------------------------------------
% Check if Ch1 or Ch2 needs to be remapped to get paired
% data for each frame for each channel
%----------------------------------------------------------
maxT = max([size(ch1_data,3), size(ch2_data,3), size(ch3_data,3)]);
ch1_data = KLS_resizeMatrix(ch1_data, maxT);
ch2_data = KLS_resizeMatrix(ch2_data, maxT);
ch3_data = KLS_resizeMatrix(ch3_data, maxT);
cell_mask = KLS_resizeMatrix(cell_mask, maxT);

t_range = 1:800;
%t_range = 450:550;
y_range = 70:100;
x_range = 79:108;


% The whole cell
t_range = 1:size(ch1_data,3);
y_range = 1:size(ch1_data,2);
x_range = 1:size(ch1_data,1);

ch1_data = ch1_data(:,:,t_range);
ch2_data = ch2_data(:,:,t_range);
ch3_data = ch3_data(:,:,t_range);
cell_mask = cell_mask(:,:,t_range);

% Generate video object
Video_Name = 'Overlay_IRM_with_Response';
vidObj = VideoWriter([Video_Name '.mp4'], 'MPEG-4');
vidObj.FrameRate = min(size(ch1_data,3) / 10, 120);
close(vidObj);
open(vidObj);
    
figure()
fig = gcf;
fig.Position = Pos+[0 -Pos(4)-90 Pos(3) Pos(4)];  
hold on


frameRepeatCount = 10;
i = 1;
while i < size(ch1_data,3)
    curr_mask = squeeze(cell_mask(:,:,i));
    curr_IRM = squeeze(ch1_data(:,:,i));
    curr_Lyso = squeeze(ch3_data(:,:,i));
    IRM_bkgd = mean(curr_IRM(~curr_mask),'all');

    curr_IRM = rescale(double(curr_IRM), 0, 1, 'InputMin', IRM_LUT_cutoffs(1), 'InputMax', IRM_LUT_cutoffs(2));
    curr_Lyso = rescale(double(curr_Lyso), 0, 1, 'InputMin', Lyso_LUT_cutoffs(1), 'InputMax', Lyso_LUT_cutoffs(2));
    
    curr_IRM_RGB = repmat(curr_IRM, [1, 1, 3]);
    
    
    if contains(lower(channel_labels.Ch3), 'blue')
        curr_Lyso_RGB = cat(3, zeros(size(curr_Lyso)), zeros(size(curr_Lyso)), curr_Lyso); % Combine only the blue channel
    elseif contains(lower(channel_labels.Ch3), 'dnd')
        curr_Lyso_RGB = cat(3, curr_Lyso, curr_Lyso * 0.5, zeros(size(curr_Lyso))); % Combine red and green channels for orange
    elseif contains(lower(channel_labels.Ch3), 'dr') || contains(lower(channel_labels.Ch3), 'deep red')
        curr_Lyso_RGB = cat(3, curr_Lyso, zeros(size(curr_Lyso)), curr_Lyso); % Combine red and blue channels for magenta
    else
        % Default magenta
        curr_Lyso_RGB = cat(3, curr_Lyso, zeros(size(curr_Lyso)), curr_Lyso); % Combine red and blue channels for magenta
    end   
    
    overlay_image = imadd(curr_IRM_RGB, curr_Lyso_RGB);
    imshow(overlay_image,[0 1],'InitialMagnification',1200)
    
    scaleLine([],10,0.157,[450 172],'hor',[1 1 1],4,14,0); % 10 um, with text
    
    text(5,5,[num2str(Tstamps(:,i),'%.2f') ' [' num2str(Tstamps_since_landed(:,i),'%.2f') '] min'],'color','white','fontsize',24)
    %{
    if any(i == 440:470)
        %annotation(gcf,'textarrow',[0.25 0.49],...
         %   [0.81 0.57],...
          %  'String',{'Rapid Drop in','Lyso Signal'},'TextColor','white');
        annotation(gcf,'arrow',[0.325 0.467],...
            [0.548 0.548],'Color','yellow','LineWidth',6,'HeadWidth',20);
    end
    if any(i == 475:535)
        %annotation(gcf,'textarrow',[0.28 0.49],...
         %   [0.90 0.57],...
        %    'String',{'Upward IRM Deformation'},'TextColor','white');
        annotation(gcf,'arrow',[0.325 0.467],...
            [0.548 0.548],'Color','yellow','LineWidth',6,'HeadWidth',20);
    end
            
       
    % Repeat frames 445 to 580 five times
    currFrame = getframe(gcf);
    if i >= 430 && i <= 540
        for repeat = 1:frameRepeatCount
            writeVideo(vidObj, currFrame);
        end
    else
        writeVideo(vidObj, currFrame);
    end
    %}
    
    currFrame = getframe(gcf);
    writeVideo(vidObj, currFrame);
    i = i+1;
end
hold off
close(vidObj);
close all
%



%% [WIP] Date__S06(II): IRM Deformation Mapped to Binding
close all
% Pull in the recorded timepoints
fileID = fopen(Time_stamps_address, 'r');

% Assuming time format in the file is like 12:34.567 (mm:ss.ms)
% '%d' reads an integer, ':%d.%f' reads the seconds and milliseconds
data = fscanf(fileID, '%d:%d.%f', [3, Inf]);

% Close the file
fclose(fileID);

% Convert data to total seconds
switch floor(max(log10(data(3,:))))
    case -Inf
            % Convert data to total seconds
            % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
            Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 10;
    case 0
        % Convert data to total seconds
        % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
        Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 10;    
    case 1
        % Convert data to total seconds
        % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
        Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 100;                   
    case 2
        % Convert data to total seconds
        % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
        Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 1000;                   
    case 3
        % Convert data to total seconds
        % data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
        Timing_seconds = data(1, :) * 60 + data(2, :) + data(3, :) / 10000;                   
end
    
cd 'D:\02_KLS_Tcell_Data\03_LyticGranule_Dynamics\20240529_SLB_LargeImage_sparseI-and-R_Fast-R-to-IRM_RegularI-R_Binding\03_SLB_60min_5nM-FOLR1-AF647-EC17_2nM-ICAM_CART_burstLyso-to-IRM\02_Landing_Burst-Lyso_I(647)_R(561)_SSL\Videos'
knocking_flag = 0;
patch_flag = 0;
remove_bkgd_IRM = 1;

ROI_n = 1;
max_z = 1500;
min_z = -1500;
cell_mask = Dilated_label_ROIs{ROI_n,1};
ch1_data = Ch1_corr_IRM_ROIs{ROI_n,1};
ch2_data = Ch2_Impulse_ROIs{ROI_n,1};
ch3_data = Ch3_Response_ROIs{ROI_n,1};

IRM_bkgd_lvl = mean(ch1_data(~cell_mask),'all');

ch1_data = ch1_data-IRM_bkgd_lvl;
%----------------------------------------------------------
% Check if Ch1 or Ch2 needs to be remapped to get paired
% data for each frame for each channel
%----------------------------------------------------------
maxT = max([size(ch1_data,3), size(ch2_data,3), size(ch3_data,3)]);
ch1_data = KLS_resizeMatrix(ch1_data, maxT);
ch2_data = KLS_resizeMatrix(ch2_data, maxT);
ch3_data = KLS_resizeMatrix(ch3_data, maxT);
cell_mask = KLS_resizeMatrix(cell_mask, maxT);
t_range = 430:514;
%t_range = 450:550;
y_range = 59:101;
x_range = 69:116;
ch1_data = ch1_data(x_range,y_range,t_range);
ch2_data = ch2_data(x_range,y_range,t_range);
ch3_data = ch3_data(x_range,y_range,t_range);
cell_mask = cell_mask(x_range,y_range,t_range);

% Generate video object
Video_Name = 'Test_IRM_warp';
vidObj = VideoWriter([Video_Name '.mp4'], 'MPEG-4');
vidObj.FrameRate = 60;
close(vidObj);
open(vidObj);
    
figure()
fig = gcf;
fig.Position = Pos+[0 -Pos(4)-90 Pos(3) Pos(4)];  
hold on

for i = 1:size(ch1_data,3)
    curr_mask = squeeze(cell_mask(:,:,i));
    curr_IRM = squeeze(ch1_data(:,:,i));
    IRM_bkgd = mean(curr_IRM(~curr_mask),'all');
    
    if remove_bkgd_IRM == 1
        curr_IRM(~curr_mask) = nan;
    end

    curr_wrap = squeeze(ch3_data(:,:,i));
    warp(curr_IRM,curr_wrap)
    %colormap(gcf,parula(500))
    colormap(gcf,'jet')
    xlim([1 size(ch1_data,2)])
    ylim([1 size(ch1_data,1)])
    zlim([min_z max_z])
    xticks('')
    yticks('')
    zticks('')
    %zlabel('IRM (AU)')
    box off
    caxis([450 7500])
    %caxis(Impulse_LUT)
    %view(gca,[-70 -50]);
    title(['Frame = ' num2str(i)])
    set(gca, 'XDir','reverse')
    set(gca, 'YDir','reverse')
    
    hold on
    if patch_flag == 1
        
            %IRM floor
            patch([1 size(ch1_data,1) size(ch1_data,1) 1], [1 1 size(ch1_data,2) size(ch1_data,2)], [IRM_bkgd IRM_bkgd IRM_bkgd IRM_bkgd],'black') 

            %Walls above floor
            x = [1 size(ch1_data,1) size(ch1_data,1) 1];
            y = [1 1 size(ch1_data,2) size(ch1_data,2)];
            z = [IRM_bkgd IRM_bkgd IRM_bkgd IRM_bkgd];

            % Define the coordinates for the sides
            x_side1 = [1 1 1 1];
            y_side1 = [1 1 size(ch1_data,2) size(ch1_data,2)];
            z_side1 = [IRM_bkgd max_z max_z IRM_bkgd];

            x_side2 = [size(ch1_data,1) size(ch1_data,1) size(ch1_data,1) size(ch1_data,1)];
            y_side2 = [1 1 size(ch1_data,2) size(ch1_data,2)];
            z_side2 = [IRM_bkgd max_z max_z IRM_bkgd];

            x_side3 = [1 size(ch1_data,1) size(ch1_data,1) 1];
            y_side3 = [1 1 1 1];
            z_side3 = [IRM_bkgd IRM_bkgd max_z max_z];

            x_side4 = [1 size(ch1_data,1) size(ch1_data,1) 1];
            y_side4 = [size(ch1_data,2) size(ch1_data,2) size(ch1_data,2) size(ch1_data,2)];
            z_side4 = [IRM_bkgd IRM_bkgd max_z max_z];

            % Draw the side patches
            patch(x_side1, y_side1, z_side1, 'black');
            patch(x_side2, y_side2, z_side2, 'black');
            patch(x_side3, y_side3, z_side3, 'black');
            patch(x_side4, y_side4, z_side4, 'black');

            % Draw the sides of the box
            for ii = 1:4
                patch(x_sides(:,ii), y_sides(:,ii), z_sides1(:,ii),'black'); % Side patches
            end
    end
    hold off
    
    switch(knocking_flag)
        case 1
            for ii = 1:10
                view(gca,[0 60]);
                    currFrame = getframe(gcf);
                    writeVideo(vidObj, currFrame);
            end
            for ii = 60:90
                view(gca,[0 ii]);
                    currFrame = getframe(gcf);
                    writeVideo(vidObj, currFrame);
            end
            for ii = 1:10
                view(gca,[0 90]);
                    currFrame = getframe(gcf);
                    writeVideo(vidObj, currFrame);
            end
            for ii = 90:-1:60
                view(gca,[0 ii]);
                    currFrame = getframe(gcf);
                    writeVideo(vidObj, currFrame);
            end
            for ii = 1:10
                view(gca,[0 60]);
                    currFrame = getframe(gcf);
                    writeVideo(vidObj, currFrame);
            end            
        case 0
            for ii = 1:1
                view(gca,[-20 34]);
                    currFrame = getframe(gcf);
                    writeVideo(vidObj, currFrame);
            end 
            
            if i == 10
                disp('annoate lysosignal drop')
                annotation(gcf,'textarrow',[0.259821428571429 0.517857142857143],...
                    [0.624 0.425],'String',{'Rapid Drop in','Lyso Signal'});
            end
            if i == 40
                disp('remove annoate lysosignal drop')
                delete(findall(gcf,'type','annotation'))
            end
            if i == 50
                disp('annotate IRM height increase')
                annotation(gcf,'textarrow',[0.4125 0.514285714285714],...
                    [0.745428571428571 0.489285714285714],...
                    'String',{'Next IRM image:','Upward IRM Deformation'});
            end
            if i == 90
                disp('remove IRM height increase')
                delete(findall(gcf,'type','annotation'))
            end
    end
end
close(vidObj);
close all

% Generate video object
Video_Name = 'Test_IRM_only';
vidObj = VideoWriter([Video_Name '.mp4'], 'MPEG-4');
vidObj.FrameRate = 60;
close(vidObj);
open(vidObj);
    
figure()
fig = gcf;
fig.Position = Pos+[0 -Pos(4)-90 Pos(3) Pos(4)];  
hold on

for i = 1:size(ch1_data,3)
    curr_mask = squeeze(cell_mask(:,:,i));
    curr_IRM = squeeze(ch1_data(:,:,i));
    IRM_bkgd = mean(curr_IRM(~curr_mask),'all');


    %imshow(curr_IRM,channel_LUTs{Contact_Channel},'InitialMagnification',600)
    imshow(curr_IRM,[min_z max_z],'InitialMagnification',600)
            if any(i == 10:40)
                annotation(gcf,'textarrow',[0.313953488372093 0.5],...
                    [0.809027777777778 0.506944444444444],...
                    'String',{'Rapid Drop in','Lyso Signal'},'TextColor','white');
            end
            if any(i == 500)
                annotation(gcf,'textarrow',[0.387596899224806 0.445736434108527],...
                    [0.767361111111111 0.513888888888889],...
                    'String',{'Next IRM Image:','Upward IRM Deformation'},'TextColor','white');
            end
    currFrame = getframe(gcf);
    writeVideo(vidObj, currFrame);
end
hold off
close(vidObj);
close all
%
% Generate video object
Video_Name = 'Test_Lyso_only';
vidObj = VideoWriter([Video_Name '.mp4'], 'MPEG-4');
vidObj.FrameRate = 60;
close(vidObj);
open(vidObj);
    
figure()
fig = gcf;
fig.Position = Pos+[0 -Pos(4)-90 Pos(3) Pos(4)];  
hold on

for i = 1:size(ch1_data,3)
    curr_mask = squeeze(cell_mask(:,:,i));
    curr_Lyso = squeeze(ch3_data(:,:,i));
    IRM_bkgd = mean(curr_IRM(~curr_mask),'all');


    %imshow(curr_IRM,channel_LUTs{Contact_Channel},'InitialMagnification',600)
    imshow(curr_Lyso,[450 7500],'InitialMagnification',600)
    colormap(gcf,'jet')
            if any(i == 500)
                annotation(gcf,'textarrow',[0.313953488372093 0.5],...
                    [0.809027777777778 0.506944444444444],...
                    'String',{'Rapid Drop in','Lyso Signal'},'TextColor','white');
            end
            
            if any(i == 500)
                annotation(gcf,'textarrow',[0.387596899224806 0.445736434108527],...
                    [0.767361111111111 0.513888888888889],...
                    'String',{'Next IRM Image:','Upward IRM Deformation'},'TextColor','white');
            end
            
    currFrame = getframe(gcf);
    writeVideo(vidObj, currFrame);
end
hold off
close(vidObj);
close all

%%
%%
%%
%%
%%
%% Section_Post_01: Import Post Data --
%---------------------------------------------------------%
% Stuff to Change
%---------------------------------------------------------%
% Where is your data? 
cd 'G:\...\FolderWithTargetImageData' % <--- Change me
% Where do you want to save data?
data_dir = dir;

loc_or_name_of_Post_data = '01_Endpoint';

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%
[Raw_Post_data, Save_individual_acq_dir, Raw_dir, Processed_dir, ROIs_dir] = IRCE_ImportFileSetup_Post(data_dir, save_data_in_this_folder, name_of_data, loc_or_name_of_Post_data);

% Pull out the metadata from the ND2 file
cd(data_dir(1).folder)
if isnumeric(loc_or_name_of_Post_data)
    meta_Post_data = KLS_ParseND2Metadata(data_dir(loc_or_name_of_Post_data).name);
elseif ischar(loc_or_name_of_Post_data)
    meta_Post_data = KLS_ParseND2Metadata([loc_or_name_of_Post_data '.nd2']);
else
    for i = 3:length(data_dir)
        if strcmp(data_dir(i).name(1:end-4), loc_or_name_of_Post_data)
            meta_Post_data = KLS_ParseND2Metadata(data_dir(i).name);
        end
    end
end

disp(' ')
disp(' ')
disp('Channel names')
for i = 1:length(meta_Post_data)
    disp(['Ch' num2str(i) ' --- ' meta_Post_data(i).Name])
end

KLS_Check_tif_Imwrite()
clear loc_of_data_in_dir curr_dir

%% Section_Post_02a: Img Process -- Raw Data and  --     
%---------------------------------------------------------%
% Stuff to Change
%---------------------------------------------------------%
num_ch_Post = length(meta_Post_data);

% Which channel denotes where the expression mark is? (i.e. T2A GFP marker)
Marker_Channel_Post = 3;

% Which channel denotes where the cell contact is?
Contact_Channel_Post = 2;

% Is this Impulse-Response Data? Which channel denotes the response?
    % set to [] if no response
Response_Channel_Post = []; % WIP currently not used KLS

% Make sure that all of the following are

% channel frequencies, list must be num_ch long 
channel_freq_Post = ones([1 num_ch_Post]); 
% channel labels, cell array must be num_ch long 
channel_labels_Post = cell([num_ch_Post 1]);
for i = 1:num_ch_Post
    channel_labels_Post{i} = meta_Post_data(i).Name;
end
% channel colors, cell array must be num_ch long, hex-code or RGB acceptable
channel_colors_Post = {}; % WIP currently not used KLS
    % '#007bff'; % Blue
    % '#ffad00'; % Orange
    % '#ff0000'; % Red
    % '#ffffff';   % White
% channel LUT, cell array must be num_ch long, individual LUTs are
    % 2-element row vectors
channel_LUTs_Post = {}; % WIP currently not used KLS

% channel numbers that you want later tracking data to be overlayed on
Impulse_track_annotation_channels_Post = []; % WIP currently not used KLS
Response_track_annotation_channels_Post = []; % WIP currently not used KLS

% binary flag indicating which channels get median image filtered, row vector must be num_ch long 
median_filter_ch_flag_Post = [1 1 0 0 0 0 0 0];
% Of the channels being median image filtered, which can serve as their own
    % background data? (i.e. a tiled data set with mutliple FOVs/timepoint), 
    % row vector must be num_ch long 
self_med_filter_ch_flag_Post = [1 1 0 0 0 0 0 0];


med_filter_lower_thresholds_Post = {[5 5]; [3 5]; []; []; []; []; []; [];};

% binary flag indicating which channels get shade corrected, row vector must be num_ch long 
shade_correct_ch_flag_Post = [0 0 0 1 1 1 1 1];
% Shade_NDFin is a 4 by 3 cell array
    % row 1-4, are t405, t488, t561, t647
    % column 1-2 are TIRF direction 0°, 45°, 135°
NDFin_or_NDFout_Post = 1; % in = 1, out = 0

% binary flag indicating which channels get converted from arbitrary 
    % fluorescence units to eletron volts (photons), row vector must be num_ch long
AU_to_Photon_flag_Post = [0 0 1 1 1 1 1 1];

% binary flag indicating which channels get photobleach corrected using
    % the mean image intensity decay as the data for a single exponential 
    % fit, row vector must be num_ch long
bleach_correct_flag_Post = [0 0 0 0 0 0 0 0];

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%
base_Post_data = cell([num_ch_Post 1]); % a cell array {n_channels by 1}
processed_Post_data = cell([num_ch_Post 1]); % a cell array {n_channels by 1}
median_Post_threshold = zeros([num_ch_Post 1]);

for i = 1:num_ch_Post
    base_Post_data{i,1} = Raw_Post_data(:,:,i:num_ch_Post*channel_freq_Post(i):end);
    [x, y, t] = size(base_Post_data{i,1});
    reshapedData = reshape(base_Post_data{i,1}, x*y, t);
    [uniqueColumns, ~, ~] = unique(reshapedData', 'rows', 'stable');
    uniqueData = reshape(uniqueColumns', x, y, size(uniqueColumns, 1));
    nonZeroFrames = any(any(uniqueData, 1), 2);
    base_Post_data{i,1} = uniqueData(:, :, nonZeroFrames);

    if channel_freq_Post(i) == 1
        channel_freq_Post(i) = round(t / size(base_Post_data{i,1},3));
    else
        ii = find(channel_freq_Post == 1);
        channel_freq_Post(i) = round(length(ii:ii*channel_freq_Post(i):size(Raw_Post_data,3)) / size(base_Post_data{i,1},3));
    end
    clear x y t reshapedData uniqueColumns uniqueData nonZeroFrames
end

%---------------------------------------------------------%
% If you need an external background reference image complete the details
% for where that data is below
%---------------------------------------------------------%
if any(median_filter_ch_flag_Post .* ~self_med_filter_ch_flag_Post)
    % If you have multiple channels in your background, let the program
    % know which slice is approperate for each data chennal you are median
    % image filtering.
    external_ch_num_Post = [1 2 0 0 0 0 0 0]; % Must be a length of num_ch
    %---------------------------------------------------------%
    % Manual Median Image Correct the IRM data
    %---------------------------------------------------------%
    
    % What folder is your IRM bkgd image in?
    cd 'G:\...\FolderWithTargetImageData' % <--- Change me
    data_dir = dir;
    loc_or_name_of_median_data = [25 26 27]; % <--- Change Me (What file do you want to import?) = 
    
    temp_cell = cell([length(loc_or_name_of_median_data) 1]);
    for i = 1:length(loc_or_name_of_median_data)
        temp_cell{i} = KLS_ND2ImportAll(data_dir(loc_or_name_of_median_data(i)).name);
    end
end

% Let's do median image correction, AU to photon convertion, shade
% correction, and blech correct as per requested above
for i = 1:num_ch_Post
    processed_Post_data{i,1} = base_Post_data{i,1};

    if median_filter_ch_flag_Post(i) == 1
        switch self_med_filter_ch_flag_Post(i)
            case 1
                [processed_Post_data{i,1}, ~, median_Post_threshold(i), ~] = KLS_RICM_bkgd_correction(processed_Post_data{i,1}, med_filter_lower_thresholds_Post{i});
            case 0
                temp_bkgd = zeros([size(temp_cell{1}, 1) size(temp_cell{1}, 2) length(loc_or_name_of_median_data)]);

                for ii = 1:length(loc_or_name_of_median_data)
                    temp_bkgd(:,:,ii) = temp_cell{ii}(:,:,external_ch_num_Post(i));
                end
                
                [~, median_img, median_Post_threshold(i), ~] = KLS_RICM_bkgd_correction(temp_bkgd, med_filter_lower_thresholds_Post{i});
                % Pull bkgd IRM image from the large data set.
                    % If no large image data available construct a background from different
                    % IRM image data sets taken that day.
                    % Check if processed_data is 512x512 in the first
                    % two dimensions, then resize median_img.
                rows = round(size(processed_Post_data{i,1},1)/512);
                cols = round(size(processed_Post_data{i,1},2)/512);
                if rows + cols > 2
                    median_img = repmat(median_img,[rows cols 1]);
                end
                processed_Post_data{i,1} = (processed_Post_data{i,1} - median_img) + round(mean(median_img,'all'));
        end
    end

    if AU_to_Photon_flag_Post(i) == 1
            gain = meta_Post_data(i).Multiplier;
            [conversion, offset] = KLS_gain_basic(gain);
    
            processed_Post_data{i,1} = (processed_Post_data{i,1}-offset) .* conversion;
    end
    
    if shade_correct_ch_flag_Post(i) == 1
        rows = round(size(base_Post_data{i,1}(:,:,1),1)/512);
        cols = round(size(base_Post_data{i,1}(:,:,1),2)/512);
        
        if floor(rows) ~= rows || floor(cols) ~= cols
            disp('Current script only handles image data that is divisible by 512.')
            return;
        end
    if isscalar(meta_Post_data(i).ExWavelength)
            switch meta_Post_data(i).ExWavelength
                case 405
                    WL_num = 1;
                case 488
                    WL_num = 2;
                case 561
                    WL_num = 3;
                case 640
                    WL_num = 4;
                otherwise
                    % Check meta_data(i).Name for wavelengths
                    if contains(meta_Post_data(i).Name, '405')
                        WL_num = 1;
                    elseif contains(meta_Post_data(i).Name, '488')
                        WL_num = 2;
                    elseif contains(meta_Post_data(i).Name, '561')
                        WL_num = 3;
                    elseif contains(meta_Post_data(i).Name, '640') || contains(meta_Post_data(i).Name, '647')
                        WL_num = 4;
                    else
                        WL_num = 2; % Default to 488
                        warning('ExWavelength not recognized. Setting WL_num to 2 (488).');
                    end
            end
        else
            % Check meta_data(i).Name for wavelengths
            if contains(meta_Post_data(i).Name, '405')
                WL_num = 1;
            elseif contains(meta_Post_data(i).Name, '488')
                WL_num = 2;
            elseif contains(meta_Post_data(i).Name, '561')
                WL_num = 3;
            elseif contains(meta_Post_data(i).Name, '640') || contains(meta_Post_data(i).Name, '647')
                WL_num = 4;
            else
                WL_num = 2; % Default to 488
                warning('ExWavelength not recognized. Setting WL_num to 2 (488).');
            end
        end
        switch meta_Post_data(i).TIRF_Direction
            case 0
                Direction_num = 1;
            case 45
                Direction_num = 2;
            case 135
                Direction_num = 3;
            otherwise
                Direction_num = 1;
        end

        if NDFin_or_NDFout_Post == 1
            shade_img = Shade_NDFin{WL_num,Direction_num};
        else
            shade_img = Shade_NDFout{WL_num,Direction_num};
        end

        Resized_shade_img = repmat(shade_img,[rows cols]);
    
        ii = 1;
        while ii <= size(processed_Post_data{i,1},3)
            processed_Post_data{i,1}(:,:,ii) = processed_Post_data{i,1}(:,:,ii)./Resized_shade_img;
            ii = ii+1;
        end
    end

    if bleach_correct_flag_Post(i) == 1 && size(processed_Post_data{i,1},3) > 3
        %---------------------------------------------------------%
        % Basic Bleach Correction using Mean Image Intensity Bleaching Rate
        %---------------------------------------------------------%
        figure()
        y = median(processed_Post_data{i,1},[1 2]);
        y = squeeze(unique(y,'stable'));
        x = 0:size(y,1)-1;

        scatter(x,y)
        [fitFunction, gof, fit_str, lifetime_tau] = KLS_Exponentialfit_and_plot(processed_Post_data{i,1}, 1);
        
        xlabel('Frame')
        if AU_to_Photon_flag_Post(i) == 1
            ylabel('Mean Intensity (Photons)')
        else
            ylabel('Mean Intensity (AU)')
        end

        legend('Data',fit_str,'location','best')
        title(['Ch_' num2str(i) ' Bleach Correction'])

        x = 0:size(processed_Post_data{i,1},3)-1;
        y = fitFunction(x);
        y = KLS_NormStack(y);

        for ii = 1:size(processed_Post_data{i,1},3)
            processed_Post_data{i,1}(:,:,ii) = processed_Post_data{i,1}(:,:,ii) ./ y(ii);
        end
    end
end

filtered_IRM_Post_data = processed_Post_data{Contact_Channel_Post,1};
IRM_thres_Post = median_Post_threshold(Contact_Channel_Post);

for i = 1:num_ch_Post
    figure('Position',Pos + (i-1)*[Pos(3) 0 0 0])
    histogram(processed_Post_data{i,1},'Normalization','PDF')
    title(['Histogram: ' channel_labels_Post{i}])

    box off
end

clear median_img Resized_shade_img temp_bkgd temp_cell conversion loc_or_name_of_median_data

%% Section_Post_03a: -- Select ROIs --   
close all
Add_more_ROIs_flag_Post = 0; % Want to add more ROIs?

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%
assay_data_FOV_rows = round(size(base_data{1,1}(:,:,1),1)/512);
assay_data_FOV_cols = round(size(base_data{1,1}(:,:,1),2)/512);

Post_data_FOV_rows = round(size(processed_Post_data{1,1}(:,:,1),1)/512);
Post_data_FOV_cols = round(size(processed_Post_data{1,1}(:,:,1),2)/512);

if floor(size(processed_Post_data{1,1}(:,:,1),1)/512) ~= Post_data_FOV_rows || floor(size(processed_Post_data{1,1}(:,:,1),2)/512) ~= Post_data_FOV_cols
    disp('Current script only handles image data that is divisible by 512.')
    return;
end
Data_FOVs_nums = [assay_data_FOV_rows assay_data_FOV_cols ...
    Post_data_FOV_rows Post_data_FOV_cols];
manual_offset = [90 -256];
[Post_roi_corners, previous_ROI_count_Post] = IRCE_CellROIs_Post(filtered_IRM_Post_data, Save_individual_acq_dir, Add_more_ROIs_flag_Post, IRM_thres_Post, roi_corners, Data_FOVs_nums, manual_offset);

%% Section_Post_03a(II): [As Needed] -- Edit -- Revise Specific ROIs
close all
Specific_ROIs = [];

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%
filtered_IRM_Post_data = processed_Post_data{Contact_Channel_Post,1};
IRM_thres_Post = median_Post_threshold(Contact_Channel_Post);

if ~isempty(Specific_ROIs)
    Post_roi_corners = IRCE_Edit_CellROIs_Post(filtered_IRM_Post_data, Save_individual_acq_dir, IRM_thres_Post, Post_roi_corners, Specific_ROIs);
end

close all

%% Section_Post_03a(III): [As Needed] -- Edit -- Remove Specific ROIs
close all
Specific_ROIs = [];

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%
filtered_IRM_Post_data = processed_Post_data{Contact_Channel_Post,1};

if ~isempty(Specific_ROIs)
    Post_roi_corners = IRCE_Remove_CellROIs(filtered_IRM_Post_data, Save_individual_acq_dir, Post_roi_corners, Specific_ROIs);
end

%% Section_Post_03b: -- Automated Mask Cells & Manually Seperate Connecting Cells --    
%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%

[Base_label_ROIs_Post] = IRCE_MaskSeperation_Post(Post_roi_corners, filtered_IRM_Post_data, IRM_thres_Post, channel_LUTs{Contact_Channel}, Save_individual_acq_dir);

%% Section_Post_03c(I): [Video Output] -- Label Video to Review Masking --   
Specific_ROIs = [];
multi_mask_list = IRCE_AutomatedMaskReview(Specific_ROIs, Base_label_ROIs_Post);

if all(isnan(multi_mask_list))
    disp('No ROIs with multiple labels, check for no mask ROIs')
else
    Specific_ROIs = multi_mask_list;
    IRCE_MaskVideo_Post(Save_individual_acq_dir, Specific_ROIs, Base_label_ROIs_Post);
end

IRCE_AutomatedMaskReview(Specific_ROIs, Base_label_ROIs_Post);

clear multi_mask_list

%% Section_Post_03c(II): [As Needed] -- Edit -- Manually Draw Specific ROI Masks -- 
cd(Save_individual_acq_dir)
Specific_ROIs = [];

if ~isempty(Specific_ROIs)
    Base_label_ROIs_Post = IRCE_RedoMaskSeperation_Post(Specific_ROIs, Post_roi_corners, filtered_IRM_Post_data, IRM_thres_Post, channel_LUTs{Contact_Channel}, Save_individual_acq_dir);

    IRCE_MaskVideo_Post(Save_individual_acq_dir, Specific_ROIs, Base_label_ROIs_Post);
else
    
end

IRCE_AutomatedMaskReview(Specific_ROIs, Base_label_ROIs_Post)

%% Section_Post_03c(III): -- Manually Connect Any Remaining Disconnected Masks --
%---------------------------------------------------------%
% Manually connect cell rois
%---------------------------------------------------------%
% KLS_generateBrackets_withOnes(n) <-- This command generates a blank set of brackets

% List out each Label # that represents the target cell
% i.e. [1 2] if the target cell is both label mask 1 and 2
% i.e. [4] if the target cell is label mask 4
% [] is equivalent to [1]
labelsForTheTargetCell_Post = {
    [1]; ...
    [1]; ...
    [1]; ...
    [1]; ...
    [1]; ...
        [1]; ...
        [1]; ...
        [1]; ...
        [1]; ...
        [1]; ... %010 
    [1]; ...
    [1]; ...
    [1]; ...
    [1]; ...
    [1]; ...
        [1]; ...
        [1]; ...
        [1]; ...
        [1]; ...
        [1]; ... %020 
    [1]; ...
    [1]; ...
    [1]; ...
    [1]; ...
    [1]; ...
        [1]; ...
        [1]; ...
        [1]; ...
        [1]; ...
        [1]; ... %030 
    [1]; ...
    [1]; ...
    [1]; ...
    [1]; ...
    [1]; ...
        [1]; ...
        [1]; ...
};

%---------------------------------------------------------%
% Combine ROIS that are one cell
%---------------------------------------------------------%
i = 1;
while i <= length(labelsForTheTargetCell_Post)
    %{
    % Check if this has section has already already been run before
    final_ROI = length(labelsForTheTargetCell_Post);
    uniqueVals = unique(Base_label_ROIs_Post{final_ROI,1}(:));
    expectedVals = [0; final_ROI];
    if isequal(sort(uniqueVals), sort(expectedVals))
        disp('Manual label assignment appears complete.');
        i = length(labelsForTheTargetCell_Post)+1;
        
        continue
    end
    %}
    
    % Extract the current set of labels
    currentSet = labelsForTheTargetCell_Post{i};
    if isempty(currentSet)
        currentSet(1) = 1;
    end

    
    % Iterate through each label in the current set, starting from the second label
    idx_target_baseLabel = (Base_label_ROIs_Post{i,1} == currentSet(1));
    
    for j = 2:length(currentSet)
        % combined idx of the target label
        idx_target_baseLabel = bitor(Base_label_ROIs_Post{i,1} == currentSet(j),idx_target_baseLabel);
    end
    % for all the pixels not associated with the target cells, set them to
        % 999
    Base_label_ROIs_Post{i,1}(bitand(~idx_target_baseLabel, Base_label_ROIs_Post{i,1}>0)) = 999;
    
    % for all the pixels that are associated with the target cells, set them to
        % ROI#, i
    Base_label_ROIs_Post{i,1}(idx_target_baseLabel) = i;
    
    i = i + 1;
end


clear i j labelsOnTheTargetCell targetLabel currentSet idx_target_dilatedLabel idx_target_baseLabel

%---------------------------------------------------------%
% Dilate the final label to yeild the dilated mask
%---------------------------------------------------------%
SE_3 = strel('disk', 6); % Flat Structuring Element for Image dilation, disk of size 9

Dilated_label_ROIs_Post = Base_label_ROIs_Post;

for  i = 1:length(Base_label_ROIs_Post)
    ii = 1;
    while ii <= size(Base_label_ROIs_Post{1},3)
        Dilated_label_ROIs_Post{i,1}(:,:,ii) = imdilate(Base_label_ROIs_Post{i,1}(:,:,ii),SE_3);
        ii = ii+1;
    end
end

%% Section_Post_03c(IV): [Video Output] -- Review Mask Assignments --    
Specific_ROIs = [];

IRCE_FinalMaskVideo_Post(Save_individual_acq_dir, Specific_ROIs, Dilated_label_ROIs_Post);


%% Section_Post_O3d: ROI channel data -- Gen channel data for each ROI --  
clc

% Pull out data for each channel corrisponding to the Cell ROIs
base_data_ROIs_Post = cell([num_ch_Post 1]); % a cell array {n_channels by 1}
processed_data_ROIs_Post = cell([num_ch_Post 1]); % a cell array {n_channels by 1}

for i = 1:num_ch_Post
    base_data_ROIs_Post{i,1} = IRCE_ROIchannelcrops(Dilated_label_ROIs_Post, base_Post_data{i,1}, Post_roi_corners);
end

for i = 1:num_ch_Post
    processed_data_ROIs_Post{i,1} = IRCE_ROIchannelcrops(Dilated_label_ROIs_Post, processed_Post_data{i,1}, Post_roi_corners);
end

[Base_label_ROIs_Post, ~] = IRCE_MaskAndROICorrnersCrop(Dilated_label_ROIs_Post, Base_label_ROIs_Post, Post_roi_corners);
[Dilated_label_ROIs_Post, Post_roi_corners] = IRCE_MaskAndROICorrnersCrop(Dilated_label_ROIs_Post, Dilated_label_ROIs_Post, Post_roi_corners);

cd(Save_individual_acq_dir)
save('Post_roi_corners.mat','Post_roi_corners','-v7.3')


%% Section_Post_05a: [Optional] [TIF Output] -- Entire FOV (Raw and Processed) Data --    
%---------------------------------------------------------%
% Save the seperate raw data from the channels
%---------------------------------------------------------%
clc

for i = 1:num_ch_Post
    file_name = ['Post_raw_' channel_labels_Post{i,1}]; % <--- Change Me as needed = 
    KLS_save_double2tif(base_Post_data{i,1}, file_name, Raw_dir);
end

%---------------------------------------------------------%
% Save the seperate processed data from the channels
%---------------------------------------------------------%
% Move to the folder
cd(Processed_dir);

for i = 1:num_ch_Post
    file_name = ['Post_processed_' channel_labels_Post{i,1}]; % <--- Change Me as needed = 
    KLS_save_double2tif(processed_Post_data{i,1}, file_name, Raw_dir);
end



%% Section_Post_05b: [Optional] [TIF Output] -- Save Cells ROI (Raw, Processed and Tracks under cell) Data  
%---------------------------------------------------------%
% Save individual Cell ROIs
%---------------------------------------------------------%
clc 

Specific_ROIs = [];

if isempty(Specific_ROIs)
    Specific_ROIs = 1:size(Base_label_ROIs_Post,1);
end

for n = 1:length(Specific_ROIs) % Loop over manually selected ROIs
    ROI_n = (Specific_ROIs(n));
    
    folderName = ['Cell_' num2str(ROI_n,'%03.f')];
    if ~exist(fullfile(ROIs_dir, folderName), 'dir')
        % If folder doesn't exist, create it
        error('Need to run timepoint data before doing endpoint.')
    end

    folderPath = fullfile(ROIs_dir, folderName);

    %---------------------------------------------------------%
    % Save Labeled Data as tif files in the cell ROIs
    %---------------------------------------------------------%    
    file_name = ['Cell_' num2str(ROI_n,'%03.f') '_Post_processed_label_dilated']; % <--- Change Me as needed = 
    KLS_save_double2tif(Dilated_label_ROIs_Post{ROI_n,1}, file_name, folderPath);
    file_name = ['Cell_' num2str(ROI_n,'%03.f') '_Post_processed_label_base']; % <--- Change Me as needed = 
    KLS_save_double2tif(Base_label_ROIs_Post{ROI_n,1}, file_name, folderPath);

    %---------------------------------------------------------%
    % Save Processed Data as tif files in the cell ROIs
    %---------------------------------------------------------%    
    for i = 1:num_ch_Post
        file_name = ['Cell_' num2str(ROI_n,'%03.f') '_Post_processed_' channel_labels_Post{i,1}]; % <--- Change Me as needed =
        KLS_save_double2tif(processed_data_ROIs_Post{i,1}{ROI_n,1}, file_name, folderPath);
    end
end
 
close all


pause(1)




%% Section_Post_08a: ROI Stats 
%---------------------------------------------------------%
% Save individual Cell ROIs
%---------------------------------------------------------%
max_ROI = size(Base_label_ROIs_Post, 1);
Invalid_ROIs = [];

cd(Save_individual_acq_dir)
folderName = 'Stats';
if ~exist(fullfile(Save_individual_acq_dir, folderName),'dir')
    % If folder doesn't exist, create it
    mkdir(fullfile(Save_individual_acq_dir, folderName));
end

folderPath = fullfile(Save_individual_acq_dir, folderName);

if isfile(fullfile(folderPath,'Stats_ROIs.mat'))
    disp('Old stats file found: loading...');
    load(fullfile(folderPath,'Stats_ROIs.mat'))
    disp('Done')
    have_old_stats_flag = 1;
else
    disp('No existing stats file found. Let''s make one shall we.')
    % Initialize Stats_ROIs cell array
    Stats_ROIs = cell(max_ROI, 1);
    have_old_stats_flag = 0;
end

ROI_n = 1;
if size(processed_data_ROIs_Post{1,1},1) == size(Stats_ROIs,1)
    while ROI_n <= max_ROI
        % Get current ROI data
        
        % Initialize stats structure
        stats = struct();
        stats.eGFP = nan(1); % GFP in units of the image
        
        % WIP
        GPF_img = processed_data_ROIs_Post{Marker_Channel_Post,1}{ROI_n}; 
        mask_img = Base_label_ROIs_Post{ROI_n,1};
        
        stats.eGFP = mean(GPF_img(logical(mask_img)),'all');
        %WIP

        % Save stats in the cell array
        if have_old_stats_flag == 1
            Stats_ROIs{ROI_n} = KLS_updateStruct(Stats_ROIs{ROI_n}, stats);
        else
            Stats_ROIs{ROI_n} = stats;
        end
    
        ROI_n = ROI_n + 1;
    end
else
    disp('Post ROIs do not match # of ROIs in the initial data')
end

for i = 1:length(Invalid_ROIs)
    if isempty(Invalid_ROIs)
        continue;
    end
    ROI_n = Invalid_ROIs(i);
    % Get current ROI data
    
    % Initialize stats structure
    stats = struct();
    stats.eGFP = nan(1); % GFP in units of the image

    % Save stats in the cell array
    if have_old_stats_flag == 1
        Stats_ROIs{ROI_n} = KLS_updateStruct(Stats_ROIs{ROI_n}, stats);
    else
        Stats_ROIs{ROI_n} = stats;
    end
end

save(fullfile(folderPath, 'Stats_ROIs.mat'),'Stats_ROIs','-v7.3')

close all



%% Section_Post_08b: [Optional] [.mat Output] -- Save to Cells ROI (individual stats files, and timing variable) Data  
%---------------------------------------------------------%
% Save individual Cell ROIs
%---------------------------------------------------------%
clc 

Specific_ROIs = [];

if isempty(Specific_ROIs)
    Specific_ROIs = 1:size(Base_label_ROIs_Post,1);
end

% Move to the folder
cd(ROIs_dir);

for n = 1:length(Specific_ROIs) % Loop over manually selected ROIs
    ROI_n = (Specific_ROIs(n));
    
    folderName = ['Cell_' num2str(ROI_n,'%03.f')];
    if ~exist(fullfile(ROIs_dir, folderName), 'dir')
        % If folder doesn't exist, create it
        error('Need to run kinetic data first.')
    end
    
    folderPath = fullfile(ROIs_dir, folderName);

    ROI_stats = Stats_ROIs{ROI_n};

    var_name = ['Cell_' num2str(ROI_n,'%03.f') '_Stats' '.mat'];

    save(fullfile(folderPath, var_name),'ROI_stats', '-v7.3')
end
 
close all

