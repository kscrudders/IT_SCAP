function [Raw_data, Save_individual_acq_dir, Raw_dir, Processed_dir, ROIs_dir] = IRCE_ImportFileSetup_Post(data_dir, save_data_in_this_folder, loc_or_name_of_data, loc_or_name_of_Post_data)
% Last update: 20240821 KLS
%     - Allow either idx in current folder or file name

    %---------------------------------------------------------%
    % Do not change this stuff
    %---------------------------------------------------------%
    %---------------------------------------------------------%
    % import data
    %---------------------------------------------------------%
    cd(data_dir(1).folder)
    if isnumeric(loc_or_name_of_data)
        folderName = data_dir(loc_or_name_of_data).name(1:end-4); % Generate a folder name for this data set
    else
        folderName = loc_or_name_of_data; % Generate a folder name for this data set
    end
    if isnumeric(loc_or_name_of_Post_data)
        Raw_data = KLS_ND2ImportAll(data_dir(loc_or_name_of_Post_data).name); 
    else
        Raw_data = KLS_ND2ImportAll([loc_or_name_of_Post_data '.nd2']); 
    end

    %---------------------------------------------------------%
    % Move cd to data parent data save folder
    %---------------------------------------------------------%
    cd(save_data_in_this_folder)
    if ~exist(fullfile(cd, folderName), 'dir')
        % If folder doesn't exist, create it
        mkdir(folderName);
    end
    % Move to the folder
    cd(folderName);
    Save_individual_acq_dir = cd;

    %---------------------------------------------------------%
    % Save the seperate raw data from the channels
    %---------------------------------------------------------%
    cd(Save_individual_acq_dir)
    folderName = 'Raw';
    if ~exist(fullfile(cd, folderName), 'dir')
        % If folder doesn't exist, create it
        mkdir(folderName);
    end
    cd(folderName);
    Raw_dir = cd;

    %---------------------------------------------------------%
    % Save the seperate processed data from the channels
    %---------------------------------------------------------%
    cd(Save_individual_acq_dir)
    folderName = 'Processed';
    if ~exist(fullfile(cd, folderName), 'dir')
        % If folder doesn't exist, create it
        mkdir(folderName);
    end
    cd(folderName);
    Processed_dir = cd;

    %---------------------------------------------------------%
    % Save the Individual ROIs
    %---------------------------------------------------------%
    cd(Save_individual_acq_dir)
        % Move to Image Folder
    folderName = 'Cell_ROIs';
    if ~exist(fullfile(cd, folderName), 'dir')
        % If folder doesn't exist, create it
        mkdir(folderName);
    end
    cd(folderName);
    ROIs_dir = cd;

    clear folderName loc_of_data_in_dir save_data_in_this_folder
end