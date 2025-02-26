function KLS_save_double2tif(data_to_save, file_name, output_dir)
    if nargin < 3
        output_dir = cd;
    end

    % Convert data to uint16 for saving
    data_uint16 = uint16(data_to_save);
    
    % Calculate the approximate size of the data in bytes
    % Each uint16 element takes 2 bytes, so we multiply the total number of elements by 2
    dataSizeInBytes = numel(data_uint16) * 2;
    sizeLimit = 4 * 1024^3; % 4 GB limit on base tif files
    
    full_path = fullfile(output_dir, [file_name '.tif']);


    if dataSizeInBytes <= sizeLimit
        % If the data size is smaller than 4 GB, use the simple tiff code
        LF_smalldata_save(data_uint16, full_path)
        %LF_bigdata_save(data_uint16, file_name)
    else
        % If the data size is larger than 4 GB, use the BigTIFF version
        LF_bigdata_save(data_uint16, full_path)
    end
end

function LF_smalldata_save(data_uint16, full_path)
    % If the data size is smaller than 4 GB, use the simple tiff code
    imwrite(data_uint16(:,:,1), full_path);
    i = 2;
    while i <= size(data_uint16, 3)
        imwrite(data_uint16(:,:,i), full_path, 'WriteMode', 'append');
        i = i + 1;
    end
end

function LF_bigdata_save(data_uint16, full_path)
    % If the data size is larger than 4 GB, use the BigTIFF version
    tiffObj = Tiff(full_path, 'w8');
    
    % Define basic tags for the TIFF file
    tagstruct.ImageLength = size(data_uint16, 1);    % Image height
    tagstruct.ImageWidth = size(data_uint16, 2);     % Image width
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack; % Greyscale image
    tagstruct.BitsPerSample = 16;                    % Each pixel is 16 bits
    tagstruct.SamplesPerPixel = 1;                   % Single channel image
    tagstruct.RowsPerStrip = 16;                     % Number of rows per strip
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression = Tiff.Compression.None;   % No compression

    % Loop over all slices and save them to the BigTIFF file
    for i = 1:size(data_uint16, 3)
        tiffObj.setTag(tagstruct);               % Set the tags for each slice
        tiffObj.write(data_uint16(:,:,i));       % Write the slice
        
        if i < size(data_uint16, 3)
            tiffObj.writeDirectory();            % Create a new directory for each slice
        end
    end

    % Close the TIFF object to finalize the file
    tiffObj.close();
end