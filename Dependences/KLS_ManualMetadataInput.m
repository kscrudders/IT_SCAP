function metadata = KLS_ManualMetadataInput()
    % Function to manually input ch_metadata using a GUI

    % Ask user for number of channels
    prompt = {'Enter number of channels:'};
    dlgtitle = 'Number of Channels';
    dims = [1 35];
    definput = {'1'};
    answer = inputdlg(prompt, dlgtitle, dims, definput);
    if isempty(answer)
        disp('No input provided. Exiting.');
        metadata = [];
        return;
    end
    numChannels = str2double(answer{1});
    if isnan(numChannels) || numChannels < 1
        error('Invalid number of channels.');
    end

    % Initialize ch_metadata structure array
    metadata = struct('Name', '', 'CameraType', '', 'Binning', '', ...
        'Exposure_ms', 0, 'Multiplier', 0, 'ReadoutSpeed_MHz', 0, ...
        'ConversionGain', '', 'VerticalShiftSpeed_us', 0, 'VerticalClockVoltageAmplitude', '', ...
        'TriggerMode', '', 'Temperature_C', 0, 'TIRF_Direction', 0, ...
        'ExWavelength', 0,'Power_percentage', 0);
    metadata = repmat(metadata, numChannels, 1);

    % Create GUI figure
    f = figure('Name', 'Manual Metadata Input', 'NumberTitle', 'off', ...
        'MenuBar', 'none', 'ToolBar', 'none', 'Resize', 'off');

    % Set figure size based on number of channels
    figWidth = 600;
    numFields = 14;
    fieldHeight = 35; % Height per field including padding
    tabHeight = numFields * fieldHeight + 100; % Additional space for buttons etc.
    figHeight = tabHeight;
    if numChannels > 1
        figHeight = tabHeight + 50; % Additional space for tabs
    end
    set(f, 'Position', [100, 100, figWidth, figHeight]);

    % Create tab group if more than one channel
    if numChannels > 1
        tgroup = uitabgroup('Parent', f);
        for ch = 1:numChannels
            tab(ch) = uitab('Parent', tgroup, 'Title', ['Channel ' num2str(ch)]);
            parent = tab(ch);
            % Create UI controls for each channel
            createChannelUI(parent, ch);
        end
    else
        % For single channel, use the figure as parent
        parent = f;
        createChannelUI(parent, 1);
    end

    % Create 'Submit' button
    submitBtn = uicontrol('Parent', f, 'Style', 'pushbutton', 'String', 'Submit', ...
        'Position', [figWidth/2 - 50, 20, 100, 40], 'Callback', @submitCallback);

    % Wait for user to close the GUI
    uiwait(f);

    % Callback function for 'Submit' button
    function submitCallback(~, ~)
        % Retrieve data from UI controls and populate ch_metadata
        for ch = 1:numChannels
            if numChannels > 1
                parent = tab(ch);
            else
                parent = f;
            end

            metadata(ch).Name = get(findobj(parent, 'Tag', ['Name' num2str(ch)]), 'String');
            metadata(ch).CameraType = getPopupValue(parent, 'CameraType', ch);
            metadata(ch).Binning = getPopupValue(parent, 'Binning', ch);
            metadata(ch).Exposure_ms = str2double(getPopupValue(parent, 'Exposure_ms', ch));
            metadata(ch).Multiplier = str2double(get(findobj(parent, 'Tag', ['Multiplier' num2str(ch)]), 'String'));
            metadata(ch).ReadoutSpeed_MHz = str2double(getPopupValue(parent, 'ReadoutSpeed_MHz', ch));
            metadata(ch).ConversionGain = getPopupValue(parent, 'ConversionGain', ch);
            metadata(ch).VerticalShiftSpeed_us = str2double(getPopupValue(parent, 'VerticalShiftSpeed_us', ch));
            metadata(ch).VerticalClockVoltageAmplitude = getPopupValue(parent, 'VerticalClockVoltageAmplitude', ch);
            metadata(ch).TriggerMode = getPopupValue(parent, 'TriggerMode', ch);
            metadata(ch).Temperature_C = str2double(get(findobj(parent, 'Tag', ['Temperature_C' num2str(ch)]), 'String'));
            metadata(ch).TIRF_Direction = str2double(getPopupValue(parent, 'TIRF_Direction', ch));
            metadata(ch).ExWavelength = str2double(getPopupValue(parent, 'ExWavelength', ch));
            metadata(ch).Power_percentage = str2double(get(findobj(parent, 'Tag', ['Power_percentage' num2str(ch)]), 'String'));
        end
        uiresume(f);
        close(f);
    end

    % Helper function to get value from popup menu
    function value = getPopupValue(parent, tagBase, ch)
        popupHandle = findobj(parent, 'Tag', [tagBase num2str(ch)]);
        contents = cellstr(get(popupHandle, 'String'));
        value = contents{get(popupHandle, 'Value')};
    end

    % Function to create UI controls for each channel
    function createChannelUI(parent, ch)
        % Define positions
        labelWidth = 200;
        editWidth = 200;
        height = 25;
        padding = 10;
        startY = figHeight - 80;
        y = startY;

        % Name
        uicontrol('Parent', parent, 'Style', 'text', 'String', 'Name:', ...
            'Position', [padding, y, labelWidth, height], 'HorizontalAlignment', 'right');
        uicontrol('Parent', parent, 'Style', 'edit', 'Tag', ['Name' num2str(ch)], ...
            'Position', [padding+labelWidth+padding, y, editWidth, height]);

        % CameraType
        y = y - height - padding;
        uicontrol('Parent', parent, 'Style', 'text', 'String', 'Camera Type:', ...
            'Position', [padding, y, labelWidth, height], 'HorizontalAlignment', 'right');
        cameraTypes = {'EMCCD', 'sCMOS', 'CCD'};
        uicontrol('Parent', parent, 'Style', 'popupmenu', 'Tag', ['CameraType' num2str(ch)], ...
            'Position', [padding+labelWidth+padding, y, editWidth, height], 'String', cameraTypes);

        % Binning
        y = y - height - padding;
        uicontrol('Parent', parent, 'Style', 'text', 'String', 'Binning:', ...
            'Position', [padding, y, labelWidth, height], 'HorizontalAlignment', 'right');
        binningOptions = {'1x1', '2x2', '4x4'};
        uicontrol('Parent', parent, 'Style', 'popupmenu', 'Tag', ['Binning' num2str(ch)], ...
            'Position', [padding+labelWidth+padding, y, editWidth, height], 'String', binningOptions);

        % Exposure_ms
        y = y - height - padding;
        uicontrol('Parent', parent, 'Style', 'text', 'String', 'Exposure (ms):', ...
            'Position', [padding, y, labelWidth, height], 'HorizontalAlignment', 'right');
        exposureOptions = {'18', '20', '50', '100', '200', '500'};
        uicontrol('Parent', parent, 'Style', 'popupmenu', 'Tag', ['Exposure_ms' num2str(ch)], ...
            'Position', [padding+labelWidth+padding, y, editWidth, height], 'String', exposureOptions);

        % Multiplier
        y = y - height - padding;
        uicontrol('Parent', parent, 'Style', 'text', 'String', 'Multiplier:', ...
            'Position', [padding, y, labelWidth, height], 'HorizontalAlignment', 'right');
        uicontrol('Parent', parent, 'Style', 'edit', 'Tag', ['Multiplier' num2str(ch)], ...
            'Position', [padding+labelWidth+padding, y, editWidth, height]);

        % ReadoutSpeed_MHz
        y = y - height - padding;
        uicontrol('Parent', parent, 'Style', 'text', 'String', 'Readout Speed (MHz):', ...
            'Position', [padding, y, labelWidth, height], 'HorizontalAlignment', 'right');
        readoutSpeeds = {'17', '15', '7'};
        uicontrol('Parent', parent, 'Style', 'popupmenu', 'Tag', ['ReadoutSpeed_MHz' num2str(ch)], ...
            'Position', [padding+labelWidth+padding, y, editWidth, height], 'String', readoutSpeeds);

        % ConversionGain
        y = y - height - padding;
        uicontrol('Parent', parent, 'Style', 'text', 'String', 'Conversion Gain:', ...
            'Position', [padding, y, labelWidth, height], 'HorizontalAlignment', 'right');
        conversionGains = {'4', '3', '2', '1'};
        uicontrol('Parent', parent, 'Style', 'popupmenu', 'Tag', ['ConversionGain' num2str(ch)], ...
            'Position', [padding+labelWidth+padding, y, editWidth, height], 'String', conversionGains);

        % VerticalShiftSpeed_us
        y = y - height - padding;
        uicontrol('Parent', parent, 'Style', 'text', 'String', 'Vertical Shift Speed (µs):', ...
            'Position', [padding, y, labelWidth, height], 'HorizontalAlignment', 'right');
        verticalShiftSpeeds = {'0.5', '0.9', '1.5', '2.0'};
        uicontrol('Parent', parent, 'Style', 'popupmenu', 'Tag', ['VerticalShiftSpeed_us' num2str(ch)], ...
            'Position', [padding+labelWidth+padding, y, editWidth, height], 'String', verticalShiftSpeeds);

        % VerticalClockVoltageAmplitude
        y = y - height - padding;
        uicontrol('Parent', parent, 'Style', 'text', 'String', 'Vertical Clock Voltage Amplitude:', ...
            'Position', [padding, y, labelWidth, height], 'HorizontalAlignment', 'right');
        vcvaOptions = {'Standard', 'Extended'};
        uicontrol('Parent', parent, 'Style', 'popupmenu', 'Tag', ['VerticalClockVoltageAmplitude' num2str(ch)], ...
            'Position', [padding+labelWidth+padding, y, editWidth, height], 'String', vcvaOptions);

        % TriggerMode
        y = y - height - padding;
        uicontrol('Parent', parent, 'Style', 'text', 'String', 'Trigger Mode:', ...
            'Position', [padding, y, labelWidth, height], 'HorizontalAlignment', 'right');
        triggerModes = {'Internal', 'External', 'Software'};
        uicontrol('Parent', parent, 'Style', 'popupmenu', 'Tag', ['TriggerMode' num2str(ch)], ...
            'Position', [padding+labelWidth+padding, y, editWidth, height], 'String', triggerModes);

        % Temperature_C
        y = y - height - padding;
        uicontrol('Parent', parent, 'Style', 'text', 'String', 'Temperature (C):', ...
            'Position', [padding, y, labelWidth, height], 'HorizontalAlignment', 'right');
        uicontrol('Parent', parent, 'Style', 'edit', 'Tag', ['Temperature_C' num2str(ch)], ...
            'Position', [padding+labelWidth+padding, y, editWidth, height]);

        % TIRF_Direction
        y = y - height - padding;
        uicontrol('Parent', parent, 'Style', 'text', 'String', 'TIRF Direction:', ...
            'Position', [padding, y, labelWidth, height], 'HorizontalAlignment', 'right');
        tirfDirections = {'0', '45', '90', '135', '180', '225', '270', '315', '360'};
        uicontrol('Parent', parent, 'Style', 'popupmenu', 'Tag', ['TIRF_Direction' num2str(ch)], ...
            'Position', [padding+labelWidth+padding, y, editWidth, height], 'String', tirfDirections);

        % ExWavelength
        y = y - height - padding;
        uicontrol('Parent', parent, 'Style', 'text', 'String', 'Excitation Wavelength (nm):', ...
            'Position', [padding, y, labelWidth, height], 'HorizontalAlignment', 'right');
        exWavelengths = {'405', '488', '561', '647'};
        uicontrol('Parent', parent, 'Style', 'popupmenu', 'Tag', ['ExWavelength' num2str(ch)], ...
            'Position', [padding+labelWidth+padding, y, editWidth, height], 'String', exWavelengths);

        % Power_percentage
        y = y - height - padding;
        uicontrol('Parent', parent, 'Style', 'text', 'String', 'Power (%):', ...
            'Position', [padding, y, labelWidth, height], 'HorizontalAlignment', 'right');
        uicontrol('Parent', parent, 'Style', 'edit', 'Tag', ['Power_percentage' num2str(ch)], ...
            'Position', [padding+labelWidth+padding, y, editWidth, height]);
    end

    
end
