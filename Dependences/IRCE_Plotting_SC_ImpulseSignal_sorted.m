function IRCE_Plotting_SC_ImpulseSignal_sorted(Stats_ROIs, fig_color, inclusion_idx, idx_map_to_time_idx)

% Updated 20250204: Impoved indexing into the timestamps based on frequency
% of impulse channel and landing point

    if isempty(inclusion_idx)
        inclusion_idx = 1:length(Stats_ROIs);
    end
    if isempty(idx_map_to_time_idx)
        idx_map_to_time_idx = 1:length(Stats_ROIs{1,1}.integrated_impulse_above_background);
    end

    if any(size(fig_color) ~= [1 3])
        try 
            fig_color = hex2rgb(fig_color); %R2024a and later
        catch
            if any(size(fig_color) ~= [1 3])
                disp('"fig_color" needs to be a RGB or hexcode value')
                return;
            end            
        end
    end
        
    inclusion_idx = find(inclusion_idx); % convert to a logical
    
    if isempty(inclusion_idx)
        disp('Must have some value to plot for "inclusion_idx"')
        return;
    end
    
    %-----%-----%
    % Find landing idx to time idx mapping
    %-----%-----%
    num_timepoints = length(Stats_ROIs{1,1}.Timing_sec);
    num_IRM_frames = length(Stats_ROIs{1,1}.Area);
    idx_mapContact_to_time_idx = 1:round(num_timepoints/num_IRM_frames):num_timepoints;


    num_ROIs = length(Stats_ROIs);

    num_frames = 0;
    for i = 1:length(Stats_ROIs)
        if length(Stats_ROIs{i,1}.integrated_impulse_above_background) > num_frames
            ROI_with_TimeStamps = i;
        end        
        num_frames = max([num_frames length(Stats_ROIs{i,1}.integrated_impulse_above_background)]);
    end
    
    % Define the base color
    base_color = fig_color;

    % Number of shades to generate
    num_shades = length(inclusion_idx);

    % Generate shades
    shades = [linspace(1, base_color(1), num_shades)', ...
              linspace(1, base_color(2), num_shades)', ...
              linspace(1, base_color(3), num_shades)'];
          
    for n = 1:length(inclusion_idx)
        i = inclusion_idx(n);

        landing_time_idx = idx_mapContact_to_time_idx(Stats_ROIs{i,1}.LandingIdx); % First frame, in contact time idxing, of contact
        first_frame = find(idx_map_to_time_idx >= landing_time_idx,1,'first'); % First frame, in impulse time idxing, after contact

        last_contact_time_idx = idx_mapContact_to_time_idx(find(~isnan(Stats_ROIs{i,1}.Area),1,'last')); % Last frame, in contact time idxing, of contact
        last_frame = find(idx_map_to_time_idx >= last_contact_time_idx,1,'first'); % Last frame, in impulse time idxing, of contact

        % pull out the impulse data for the contact
        y = Stats_ROIs{i,1}.integrated_impulse_above_background;
        y = y(first_frame:last_frame);

        if isempty(y)
            continue
        end

        x = Stats_ROIs{i,1}.Timing_sec;
        % convert first impulse data idx to time idx
        interval_x = round(mean(diff(idx_map_to_time_idx)));
        x = x(idx_map_to_time_idx(first_frame):interval_x:idx_map_to_time_idx(last_frame));
        x = x-x(1); % Reorder starting at 0
        x = x./60; % convert frames to minutes
       
        hold on
            plot(x,y, 'Color', shades(n,:));
        hold off
    end
end