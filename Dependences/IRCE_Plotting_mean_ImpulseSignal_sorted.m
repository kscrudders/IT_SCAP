function IRCE_Plotting_mean_ImpulseSignal_sorted(Stats_ROIs, fig_color, inclusion_idx, idx_map_to_time_idx)
%

% Updated 20240822: KLS zero values once a cell lands are valid signals,
    % include them in the definition for a valid timepoint (i.e a timepoint
    % with more than 2 cells with signal)
% Updated 20240919: KLS let's exclude zeros from the mean. This changes the
    % meaning of the plot. Instead of the mean intensity *of the population of
    % cells* its the mean intensity of cells with signal.
% Updated 20250204: Impoved indexing into the timestamps based on frequency
% of impulse channel

    if isempty(inclusion_idx)
        inclusion_idx = 1:length(Stats_ROIs);
    end
    if isempty(idx_map_to_time_idx)
        idx_map_to_time_idx = 1:length(Stats_ROIs{1,1}.integrated_impulse_above_background);
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

    num_ROIs = length(inclusion_idx);

    num_frames = 0;
    for i = 1:length(Stats_ROIs)
        if length(Stats_ROIs{i,1}.integrated_impulse_above_background) > num_frames
            ROI_with_TimeStamps = i;
        end        
        num_frames = max([num_frames length(Stats_ROIs{i,1}.integrated_impulse_above_background)]);
    end

    comb_y = nan([num_ROIs num_frames]);
    for n = 1:length(inclusion_idx)
        i = inclusion_idx(n);
        
        landing_time_idx = idx_mapContact_to_time_idx(Stats_ROIs{i,1}.LandingIdx); % First frame, in contact time idxing, of contact
        first_frame = find(idx_map_to_time_idx >= landing_time_idx,1,'first'); % First frame, in impulse time idxing, after contact

        last_contact_time_idx = idx_mapContact_to_time_idx(find(~isnan(Stats_ROIs{i,1}.Area),1,'last')); % Last frame, in contact time idxing, of contact
        last_frame = find(idx_map_to_time_idx >= last_contact_time_idx,1,'first'); % Last frame, in impulse time idxing, of contact
        
        % pull out the impulse data for the contact
        y = Stats_ROIs{i,1}.integrated_impulse_above_background;
        y = y(first_frame:last_frame);


        % Find the length of the cell trace
        data_width = last_frame-first_frame;

        %Insert the current cell data into a combined variable 1:length cell
            %trace
        comb_y(i,1:data_width+1) = y';
    end

    % remove data once only 2 cells trace remain in time
    % non_zero_comb_y = comb_y > 0;
    non_zero_comb_y = ~isnan(comb_y);

    comb_y(comb_y == 0) = nan; % 20240919 set values inside cell trace that are zero to nan to exclude them

    cells_at_diff_points = sum(non_zero_comb_y,1);
    valid_idx = cells_at_diff_points > 2;

    hold on
        x = Stats_ROIs{ROI_with_TimeStamps,1}.Timing_sec; % Start with the base timeline (for the most frequent channel)
        x = x(idx_map_to_time_idx); % Pull out impulse timepoints
        x = x(1:size(comb_y,2)) ./ 60; % trancate to size of combined y and convert to minutes
        
        y = mean(comb_y,1,'omitnan');
        y_sem = SEM_calc(comb_y,0.05); % Calculate 95% CI using SEM, assumes guassian population

        y_plus_error = y + y_sem;
        y_minus_error = y - y_sem;

        x_error = [x, fliplr(x)];
        AreainBetween_errors = [y_plus_error, fliplr(y_minus_error)];
    hold on    
        opacity_fig = 0.4;
        edge_opacity = 0.8;
        plot_color = fig_color;
        fill(x_error([valid_idx,fliplr(valid_idx)]), AreainBetween_errors([valid_idx,fliplr(valid_idx)]),'g','FaceAlpha',opacity_fig,'FaceColor',plot_color,'EdgeAlpha',edge_opacity);
        plot(x(valid_idx),y(valid_idx),'Color',plot_color,'LineWidth',4)
    hold off
end