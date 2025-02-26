function ch_metadata = KLS_ParseND2Metadata(FileAddress)
    % Dependence: ND2Info
    % https://www.mathworks.com/matlabcentral/fileexchange/71345-nd2read
    [ImageInfo] = ND2Info(FileAddress);
    
    % Initialize a structure array for channel metadata
    ch_metadata = struct('Name', '', 'CameraType', '', 'Binning', {}, ...
        'Exposure_ms', 0, 'Multiplier', 0, 'ReadoutSpeed_MHz', 0, ...
        'ConversionGain', {}, 'VerticalShiftSpeed_s', {}, 'VerticalClockVoltageAmplitude', {}, ...
        'TriggerMode', '', 'Temperature_C', 0, 'TIRF_Direction', 0, ...
        'ExWavelength', 0,'Power_percentage', 0);
    
    metadata_text = ImageInfo.description;
    
    % Find all plane sections
    plane_indices = regexp(metadata_text, 'Plane #\d+:', 'start');
    if isempty(plane_indices)
        plane_indices = 1;
    end
    
    % Loop through each plane section to extract metadata
    for i = 1:ImageInfo.Component
        if i == ImageInfo.Component
            % Last plane, get until the end of the text
            plane_text = metadata_text(plane_indices(i):end);
        else
            % Get text from current plane to the start of the next plane
            plane_text = metadata_text(plane_indices(i):plane_indices(i+1)-1);
        end
        
        % Extract channel name
        ch_metadata(i).Name = ImageInfo.metadata.channels(i).channel.name;
        
        % Extract other fields
        temp_cell = strtrim(extractBetween(plane_text, 'Camera Type:', 'Binning:'));
        ch_metadata(i).CameraType = temp_cell{1};
        temp_cell = strtrim(extractBetween(plane_text, 'Binning:', 'Exposure:'));
        ch_metadata(i).Binning = temp_cell;
        ch_metadata(i).Exposure_ms = str2double(extractBetween(plane_text, 'Exposure:', ' ms'));
        ch_metadata(i).Multiplier = str2double(extractBetween(plane_text, 'Multiplier:', 'Readout Speed:'));
        ch_metadata(i).ReadoutSpeed_MHz = str2double(extractBetween(plane_text, 'Readout Speed:', ' MHz'));
        ch_metadata(i).ConversionGain = strtrim(extractBetween(plane_text, 'Conversion Gain:', 'Vertical Shift Speed:'));
        ch_metadata(i).VerticalShiftSpeed_s = str2double(extractBetween(plane_text, 'Vertical Shift Speed:', 's'));
        ch_metadata(i).VerticalClockVoltageAmplitude = strtrim(extractBetween(plane_text, 'Vertical Clock Voltage Amplitude:', 'Trigger Mode:'));
        ch_metadata(i).TriggerMode = strtrim(extractBetween(plane_text, 'Trigger Mode:', 'Temperature:'));
        ch_metadata(i).Temperature_C = str2double(extractBetween(plane_text, 'Temperature:', 'C'));
        ch_metadata(i).TIRF_Direction = str2double(extractBetween(plane_text, 'Direction: ', 'LUN-F,'));
        ch_metadata(i).ExWavelength = str2double(extractBetween(plane_text, 'ExW:', ';'));
        ch_metadata(i).Power_percentage = str2double(extractBetween(plane_text, 'Power:', ';'));
    end
    
    %{
    text = ImageInfo.capturing;

    % Split the text into samples
    chParts = regexp(text, 'Sample \d+:', 'split');
    numCh = length(chParts) - 1; % First part is before 'Sample 1:'

    % Initialize structure array
    ch_metadata(numCh) = struct();

    for i = 1:numCh
        % Extract sample text
        chText = chParts{i+1};

        % Parse metadata
        ch_metadata(i).Name = ImageInfo.metadata.channels(i).channel.name;
        ch_metadata(i).CameraType = strtrim(extractBetween(chText, 'Camera Type:', 'Binning:'));
        ch_metadata(i).Binning = strtrim(extractBetween(chText, 'Binning:', 'Exposure:'));
        ch_metadata(i).Exposure_ms = str2double(extractBetween(chText, 'Exposure:', ' ms'));
        ch_metadata(i).Multiplier = str2double(extractBetween(chText, 'Multiplier:', 'Readout Speed:'));
        ch_metadata(i).ReadoutSpeed_MHz = str2double(extractBetween(chText, 'Readout Speed:', ' MHz'));
        ch_metadata(i).ConversionGain = strtrim(extractBetween(chText, 'Conversion Gain:', 'Vertical Shift Speed:'));
        ch_metadata(i).VerticalShiftSpeed_s = str2double(extractBetween(chText, 'Vertical Shift Speed:', ' s'));
        ch_metadata(i).VerticalClockVoltageAmplitude = strtrim(extractBetween(chText, 'Vertical Clock Voltage Amplitude:', 'Trigger Mode:'));
        ch_metadata(i).TriggerMode = strtrim(extractBetween(chText, 'Trigger Mode:', 'Temperature:'));
        ch_metadata(i).Temperature_C = str2double(extractBetween(chText, 'Temperature:', 'C'));
    end
    %}
end

