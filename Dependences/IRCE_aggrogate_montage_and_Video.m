function [no_output_copy_files] = IRCE_aggrogate_montage_and_Video(rootDir,targetDirMontage,targetDirVideo)
    % Search for all folders in the root directory
    folders = dir(rootDir);
    folders = folders([folders.isdir]); % Filter out anything that's not a directory

    %---------------------------------------------------------%
    % Videos
    %---------------------------------------------------------%
    video_names = {'*Response_Video.mp4', ...
        '*_Impulse_Response_Video_*.mp4', ...
        '*_IRM_Binding_Video.mp4', ...
        '*_Channel_Video*', ...
        };
    % Loop through each folder to find and copy video files ending in Response_Video.avi
    for z = 1:length(video_names)
        for i = 1:length(folders)
            % Skip the '.' and '..' folders
            if strcmp(folders(i).name, '.') || strcmp(folders(i).name, '..')
                continue;
            end

            % Construct the folder path
            folderPath = fullfile(rootDir, folders(i).name);

            files = dir(fullfile(folderPath, video_names{1,z}));
            % Copy each found file to the target directory
            for j = 1:length(files)
                srcFile = fullfile(folderPath, files(j).name);
                dstFile = fullfile(targetDirVideo, files(j).name);
                copyfile(srcFile, dstFile);
            end
        end        
    end
    
    %---------------------------------------------------------%
    % Montages
    %---------------------------------------------------------%
    Montage_names = {'*_Impulse_Response_*.png', ...
        '*_Impulse__Response_*.png', ...
        '*Response_Montage.png', ...
        '*_IRM_Binding_Montage.png', ...
        };
    % Loop through each folder to find and copy video files ending in Response_Video.avi
    for z = 1:length(Montage_names)
        % Loop through each folder to find and copy montage images with
        for i = 1:length(folders)
            % Skip the '.' and '..' folders
            if strcmp(folders(i).name, '.') || strcmp(folders(i).name, '..')
                continue;
            end

            % Construct the folder path
            folderPath = fullfile(rootDir, folders(i).name);

            files = dir(fullfile(folderPath, Montage_names{1,z}));
            % Copy each found file to the target directory
            for j = 1:length(files)
                srcFile = fullfile(folderPath, files(j).name);
                dstFile = fullfile(targetDirMontage, files(j).name);
                copyfile(srcFile, dstFile);
            end
        end   
    end

    disp('Files copied successfully.');
end

