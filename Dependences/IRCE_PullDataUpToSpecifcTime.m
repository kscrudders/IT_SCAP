function parsed_variable = IRCE_PullDataUpToSpecifcTime(stats, time_in_min)
    % Copy the input variable to the output
    parsed_variable = stats;

    % Convert Timing_sec to minutes
    Timing_min = stats.Timing_sec / 60;

    % Calculate relative timing to landing
    relative_timing = Timing_min - Timing_min(stats.LandingIdx);

    % Find the first index where relative timing >= time_in_min
    end_idx = find(round(relative_timing) >= time_in_min, 1);

    % If no valid end index is found, return nan
    if isempty(end_idx)
        parsed_variable = nan;
        disp(['Contact duraction = ' num2str(relative_timing(end)) ' min'])
        return;
    end

    % Define variables to index
    variable_names = fieldnames(stats);

    % Index the variables with the same length as Timing_sec
    for i = 1:numel(variable_names)
        var_name = variable_names{i};
        current_value = stats.(var_name);

        % Check if the variable has the same length as Timing_sec
        if isnumeric(current_value) && isequal(size(current_value, 1), numel(stats.Timing_sec))
            parsed_variable.(var_name) = current_value(stats.LandingIdx:end_idx, :);
        elseif isnumeric(current_value) && isequal(size(current_value, 2), numel(stats.Timing_sec))
            parsed_variable.(var_name) = current_value(:, stats.LandingIdx:end_idx);
        end
    end
end