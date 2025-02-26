function time_str = KLS_format_seconds_to_time_string(Tstamp_final, Tstamp_max)
    % Input: Tstamp_final (a single number value in seconds)
    % Output: time_str (a string repsenting the input in 
        % HH:MM:ss.ms format)
    if Tstamp_final < 0
        % Get total seconds
        total_seconds = abs(Tstamp_final);
        
        % Calculate hours, minutes, and seconds
        hours = floor(total_seconds / 3600);
        total_seconds = mod(total_seconds, 3600);
        minutes = floor(total_seconds / 60);
        seconds = floor(mod(total_seconds, 60));
        
        % Get milliseconds
        num = (total_seconds - seconds) ;
        integ = fix(num);
        fract = abs(num - integ);
        milliseconds = round(fract * 1000);
        
        % Ensure milliseconds have leading zeros
        % Convert milliseconds to string
        ms_str = num2str(milliseconds);
        
        % Remove trailing zeros from ms_str
        if milliseconds ~= 0
            ms_str = regexprep(ms_str, '0+$', '');
        end
        
        
        if Tstamp_max > 3600
            % Combine into HH:MM:SS.ms format
            time_str = sprintf('%02d:%02d:%02d.%s', hours, minutes, seconds, ms_str);
        elseif Tstamp_max > 60
            % Combine into MM:SS.ms format
            time_str = sprintf('%02d:%02d.%s', minutes, seconds, ms_str);
        else
            % Combine into SS.ms format
            time_str = sprintf('%02d.%s', seconds, ms_str);
        end
        time_str = ['-' time_str];
        
    else
        % Get total seconds
        total_seconds = Tstamp_final;
        
        % Calculate hours, minutes, and seconds
        hours = floor(total_seconds / 3600);
        total_seconds = mod(total_seconds, 3600);
        minutes = floor(total_seconds / 60);
        seconds = floor(mod(total_seconds, 60));
        
        % Get milliseconds
        num = (total_seconds - seconds) ;
        integ = fix(num);
        fract = abs(num - integ);
        milliseconds = round(fract * 1000);
        
        % Ensure milliseconds have leading zeros
        % Convert milliseconds to string
        ms_str = num2str(milliseconds);
        
        % Remove trailing zeros from ms_str
        if milliseconds ~= 0
            ms_str = regexprep(ms_str, '0+$', '');
        end
        
        
        if Tstamp_max > 3600
            % Combine into HH:MM:SS.ms format
            time_str = sprintf('%02d:%02d:%02d.%s', hours, minutes, seconds, ms_str);
        elseif Tstamp_max > 60
            % Combine into MM:SS.ms format
            time_str = sprintf('%02d:%02d.%s', minutes, seconds, ms_str);
        else
            % Combine into SS.ms format
            time_str = sprintf('%02d.%s', seconds, ms_str);
        end
    end
end