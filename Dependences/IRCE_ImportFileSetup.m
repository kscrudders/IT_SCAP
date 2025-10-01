function [Raw_data, base_save_dir, Raw_dir, Processed_dir, ROIs_dir] = IRCE_ImportFileSetup(data_dir, save_data_in_this_folder, name_of_data)
% Last update: 20241203 KLS
%    - Check for import file existance and allow ND2 or tif format,
%       preference ND2
% 20240821 KLS
%   - Allow either idx in current folder or file name
% 20250516 KLS
%   - Also generate the stats folder

    base_save_dir = [save_data_in_this_folder '\' name_of_data];
    %---------------------------------------------------------%
    % Do not change this stuff
    %---------------------------------------------------------%

    %---------------------------------------------------------%
    % import data
    %---------------------------------------------------------%
    % Construct the full file path for the .nd2 file
    nd2_file = fullfile(data_dir, [name_of_data '.nd2']);
    
    % Construct the full file path for the .tif file
    tif_file = fullfile(data_dir, [name_of_data '.tif']);
    
    % Check if the .nd2 file exists
    if isfile(nd2_file)
        % If .nd2 file exists, import using KLS_ND2ImportAll
        Raw_data = KLS_ND2ImportAll(nd2_file);
    elseif isfile(tif_file)
        % If .nd2 file does not exist but .tif file does, import using KLS_TIFImportAll
        Raw_data = KLS_TifImportAll(tif_file);
    else
        % If neither file exists, throw an error
        error('That file does not exists as a .tif nor .nd2 in the specified directory.');
    end
    
    %---------------------------------------------------------%
    % This is where all the workflow's data will be saved
    %---------------------------------------------------------%
    % Check if the folder exists, and if not, create it
    if ~exist(base_save_dir, 'dir')
        % If the folder doesn't exist, create it
        mkdir(base_save_dir);
    end

    %---------------------------------------------------------%
    % Save the seperate raw data from the channels
    %---------------------------------------------------------%
    folderName = 'Raw';
    folder_path = fullfile(base_save_dir, folderName);
    if ~exist(folder_path, 'dir')
        % If folder doesn't exist, create it
        mkdir(folder_path);
    end
    Raw_dir = folder_path;

    %---------------------------------------------------------%
    % Save the seperate processed data from the channels
    %---------------------------------------------------------%
    folderName = 'Processed';
    folder_path = fullfile(base_save_dir, folderName);
    if ~exist(folder_path, 'dir')
        % If folder doesn't exist, create it
        mkdir(folder_path);
    end
    Processed_dir = folder_path;

    %---------------------------------------------------------%
    % Save the Individual ROIs
    %---------------------------------------------------------%
    folderName = 'Cell_ROIs';
    folder_path = fullfile(base_save_dir, folderName);
    if ~exist(folder_path, 'dir')
        % If folder doesn't exist, create it
        mkdir(folder_path);
    end
    ROIs_dir = folder_path;

    %---------------------------------------------------------%
    % Save the Stats
    %---------------------------------------------------------%
    folderName = 'Stats';
    folder_path = fullfile(base_save_dir, folderName);
    if ~exist(folder_path,'dir')
        % If folder doesn't exist, create it
        mkdir(folder_path);
    end
end