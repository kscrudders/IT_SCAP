function Timing_seconds = IRCE_ND2_TimeStamps(Time_stamps_address)
    % Pull in the recorded timepoints
    fileID = fopen(Time_stamps_address, 'r');

    % Assuming time format in the file is like 12:34.567 (mm:ss.ms)
    % '%d' reads an integer, ':%d.%f' reads the seconds and milliseconds
    data = fscanf(fileID, '%d:%d.%d', [3, Inf]);

    % Read all lines into a cell array
    data_lines = textscan(fileID, '%s', 'Delimiter', '\n');
    data_lines = data_lines{1};
    if isempty(data_lines)
        line = [''];
    else
        line = data_lines{1};
    end

    if (contains(line, ':') || isempty(data_lines)) && any(data > 0,'all') % Time format is likely mm:ss.ms
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
    else % Time format is likely ss.ms
        %{
        % Open the file
        fileID = fopen(Time_stamps_address, 'r');
        
        % Read all lines into a cell array
        data_lines = textscan(fileID, '%s', 'Delimiter', '\n');
        data_lines = data_lines{1};
        
        % Close the file
        fclose(fileID);
        
        % Initialize the output array
        Timing_seconds = zeros(length(data_lines), 1);
        
        % Loop through each line to parse the time
        for i = 1:length(data_lines)
            line = data_lines{i};
            if contains(line, ':')
                % Time format is mm:ss.ms
                % Use sscanf to parse minutes and seconds
                time_parts = sscanf(line, '%d:%f');
                minutes = time_parts(1);
                seconds = time_parts(2);
                total_seconds = minutes * 60 + seconds;
            else
                % Time format is ss.ms
                % Directly convert the string to a number
                total_seconds = str2double(line);
            end
            Timing_seconds(i) = total_seconds;
        end
        %}
        
        % Open the file
        fileID = fopen(Time_stamps_address, 'r');
        
        % Read all lines into a cell array
        data_lines = textscan(fileID, '%s', 'Delimiter', '\n');
        data_lines = data_lines{1};
        
        % Close the file
        fclose(fileID);
        
        % Initialize the output array
        Timing_seconds = zeros(length(data_lines), 1);
        
        % Loop through each line to parse the time
        for i = 1:length(data_lines)
            line = data_lines{i};
            colon_count = count(line, ':');
            
            if colon_count == 2
                % Time format is hh:mm:ss.ms
                % Use sscanf to parse hours, minutes, and seconds
                time_parts = sscanf(line, '%d:%d:%f');
                hours = time_parts(1);
                minutes = time_parts(2);
                seconds = time_parts(3);
                total_seconds = hours * 3600 + minutes * 60 + seconds;
            elseif colon_count == 1
                % Time format is mm:ss.ms
                % Use sscanf to parse minutes and seconds
                time_parts = sscanf(line, '%d:%f');
                minutes = time_parts(1);
                seconds = time_parts(2);
                total_seconds = minutes * 60 + seconds;
            else
                % Time format is ss.ms
                % Directly convert the string to a number
                total_seconds = str2double(line);
            end
            Timing_seconds(i) = total_seconds;
        end
    end 
end