%%
%{
Work in process sections for specific tasks
%}

%%

%%

%% Section_04b: [As Needed] -- IRM Deformation Mask Automated --    
% WIP as of 20240716 - KLS
[roi_IRM_interface] = IRCE_IRMDeformationAnnotation_Automated(Save_individual_acq_dir, roi_corners, filtered_IRM_data, Base_label_ROIs, channel_LUTs{Contact_Channel}, Dilated_label_ROIs,median_threshold(Contact_Channel));

%% Section_06a(i): [Optional - TrackMate Tracking] -- Aggregate files for batch TrackMate Tracking 
clc 
%file_name_pattern = '*_processed_LysoBrite Blue.tif';
%file_name_pattern = '*_processed_forTracking_LysoTracker DND-99.tif';

file_name_pattern = '*_processed_forTracking_Binding tirf 647.tif';
file_name_pattern = '*_processed_forTracking_LysoTracker DND-99.tif';
file_name_pattern = '*_processed_forTracking_FOLR1 Binding.tif';
file_name_pattern = '*_processed_forTracking_LysoBrite Blue.tif';

file_name_pattern = '*_processed_Synced_Binding tirf 647.tif';
file_name_pattern = '*_processed_Synced_LysoTracker DND-99.tif';
file_name_pattern = '*_processed_Synced_FOLR1 Binding.tif';
file_name_pattern = '*_processed_Synced_LysoBrite Blue.tif';

file_name_pattern = '*_processed_Synced_LysoTracker DeepRed.tif';


%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%
mask_name_pattern = '*_processed_label_dilated.tif';

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

    % Search for files with a specific pattern in the current folder
    files = dir(fullfile(folderPath, mask_name_pattern));
    
    % Copy each found file to the target directory and store its original path
    for j = 1:length(files)
        srcFile = fullfile(folderPath, files(j).name);
        dstFile = fullfile(aggregatedDir, files(j).name);
        copyfile(srcFile, dstFile);
        
        % Store the original path with the file name as the key
        originalPaths(files(j).name) = folderPath;
    end
end

disp('Base files copied successfully.');

% This section of the script processes paired data and label stacks
% in the aggregatedDir: applies the label mask to the image stack,
% saves the masked image over the original synced file, and deletes the label file.

% Define the directory where paired files were copied
aggregatedDir = fullfile(Save_individual_acq_dir, 'Aggregated');

% Find all synced files in the aggregated directory
imgFiles = dir(fullfile(aggregatedDir, file_name_pattern));

for k = 1:length(imgFiles)
    % Get synced file name and parse prefix (e.g., 'Cell_001')
    imgName = imgFiles(k).name;

    % Construct corresponding label file name
    nameParts = split(imgName, '_processed_');
    prefix = nameParts{1}; 
    
    labelName = sprintf('%s_processed_label_dilated.tif', prefix);
    
    % Full paths
    syncedPath = fullfile(aggregatedDir, imgName);
    labelPath  = fullfile(aggregatedDir, labelName);
    
    % Import stacks
    img_stack   = KLS_ND2ImportAll(syncedPath);
    label_stack = KLS_ND2ImportAll(labelPath);
    
    % Determine numeric ID from prefix (e.g., '001' -> 1)
    idStr = regexp(prefix, '\d+', 'match', 'once');
    idNum = str2double(idStr);
    
    % Create binary mask for this object
    mask_stack = label_stack == idNum;
    
    % Apply mask to image stack
    img_stack_out = img_stack .* mask_stack;
    
    % Save masked image, overwriting the original synced file
    % KLS_save_double2tif(data, baseName, outputDir)
    baseName = imgName(1:end-4);  % remove .tif
    KLS_save_double2tif(img_stack_out, baseName, aggregatedDir);
    
    % Delete the label file
    if exist(labelPath, 'file')
        delete(labelPath);
    end
end

disp('Masked image stacks saved and label files removed.');

clear folders i folderPath j srcFile dstFile copyfile file_name_pattern




%% Section_06a(ii): [Optional - TrackMate Tracking] -- Add empty trackfiles | Rename | Move TrackMate Tracking files
% Need to have run section_06a(i)
clc

% Define file_name_pattern and file_name_pattern_new as cell arrays
file_name_pattern = {
    '*_processed_forTracking_Binding tirf 647-all-spots.csv',  ...
    '*_processed_forTracking_LysoBrite Blue-all-spots.csv', ...
    '*_processed_forTracking_t647 Short Expo. Binding-all-spots.csv', ...
    '*_processed_forTracking_LysoTracker DND-99-all-spots.csv', ...
    '*_processed_forTracking_FOLR1 Binding-all-spots.csv', ...
    ...
    '*_processed_Synced_Binding tirf 647-all-spots.csv',  ...
    '*_processed_Synced_LysoBrite Blue-all-spots.csv', ...
    '*_processed_Synced_t647 Short Expo. Binding-all-spots.csv', ...
    '*_processed_Synced_LysoTracker DND-99-all-spots.csv', ...
    '*_processed_Synced_FOLR1 Binding-all-spots.csv', ...
    '*_processed_Synced_LysoTracker DR-all-spots.csv', ...
};

file_name_pattern_new = {
    '_All_Spots_Binding.csv', ...
    '_All_Spots_Response.csv', ...
    '_All_Spots_Binding.csv', ...
    '_All_Spots_Response.csv', ...
    '_All_Spots_Binding.csv', ...
    ...
    '_All_Spots_Binding.csv', ...
    '_All_Spots_Response.csv', ...
    '_All_Spots_Binding.csv', ...
    '_All_Spots_Response.csv', ...
    '_All_Spots_Binding.csv', ...
    '_All_Spots_Response.csv', ...
};

%---------------------------------------------------------%
% Do not change remaining code in this section
%---------------------------------------------------------%

% Define the template CSV file path
templateCSVPath = 'I:\010_SC_Data\Cell_001_emptyTrackMate.csv';

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

%% Section_06b: Tracking Import

maxT = length(Time_stamps);

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

% Convert data to total seconds
% data(1, :) are minutes, data(2, :) are seconds, data(3, :) are milliseconds
Tstamps = Time_stamps/60; % convert s to min

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
    
Timing_seconds = Time_stamps;

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
    
Timing_seconds = Time_stamps;

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
        Time_stamps, Video_Name, ...
        ROI_n, roi_IRM_interface(ROI_n,:,:), temp_ch_labels); 
end