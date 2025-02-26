function meta_data = KLS_GetMetadata(base_save_dir)
    % Define the path to the ch_metadata file
    metadata_file = fullfile(base_save_dir, 'meta_data.mat');
    
    if exist(metadata_file, 'file')
        % If the metadata file exists, load ch_metadata from it
        loaded_data = load(metadata_file, 'meta_data');
        meta_data = loaded_data.meta_data;
        disp('Loaded existing meta_data from file.');
    else
        % If the file does not exist, run the manual input GUI
        meta_data = KLS_ManualMetadataInput();
        % Save ch_metadata to the specified directory
        save(metadata_file, 'meta_data');
        disp('Saved meta_data to file.');
    end
end