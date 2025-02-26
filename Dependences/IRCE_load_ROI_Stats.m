function loadedData = IRCE_load_ROI_Stats(rootDir)
    % Search for all "Stats_ROIs.mat" files in the directory and its subfolders
    files = dir(fullfile(rootDir, '**', 'Stats_ROIs.mat'));

    % Initialize a cell array to hold the loaded data
    loadedData = cell([length(files) 2]);

    % Loop through each found file and load its content into the cell array
    for k = 1:length(files)
        % Get the full path of the file
        filePath = fullfile(files(k).folder, files(k).name);

        % Load the file
        temp = load(filePath);
        data = temp.Stats_ROIs;

        % Extract the parent folder name
        relativePath = strrep(filePath, [rootDir filesep], '');
        subFolders = strsplit(relativePath, filesep);
        firstSubFolder = subFolders{1};

        % Store the loaded data and the parent folder name in the cell array
        loadedData{k, 1} = firstSubFolder;  % Store the parent folder name
        loadedData{k, 2} = data;  % Store the loaded data
    end

    clear k filePath files firstSubFolder relativePath rootDir temp subFolders data
end