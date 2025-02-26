function IRCE_Plotting_mean_Ch3Signal_sorted_withLocalizationFilter(Stats_ROIs, fig_color, inclusion_idx)
%

% Updated 20240822: KLS zero values once a cell lands are valid signals,
    % include them in the definition for a valid timepoint (i.e a timepoint
    % with more than 2 cells with signal)
% Updated 20240919: KLS let's exclude zeros from the mean. This changes the
% meaning of the plot. Instead of the mean intensity *of the population of
% cells* its the mean intensity of cells with signal.

    if nargin < 3
        inclusion_idx = 1:length(Stats_ROIs);
    end
    
    inclusion_idx = find(inclusion_idx); % convert to a logical
    
    if isempty(inclusion_idx)
        disp('Must have some value to plot for "inclusion_idx"')
        return;
    end

    num_ROIs = length(inclusion_idx);

    num_frames = 0;
    for i = 1:length(Stats_ROIs)
        if length(Stats_ROIs{i,1}.integrated_response_above_background_guass_withlocalization) > num_frames
            ROI_with_TimeStamps = i;
        end        
        num_frames = max([num_frames length(Stats_ROIs{i,1}.integrated_response_above_background_guass_withlocalization)]);
    end

    comb_y = nan([num_ROIs num_frames]);
    for n = 1:length(inclusion_idx)
        i = inclusion_idx(n);
        
        first_frame = find(~isnan(Stats_ROIs{i,1}.Area),1,'first');
        last_frame = find(~isnan(Stats_ROIs{i,1}.Area),1,'last');

        y = Stats_ROIs{i,1}.integrated_response_above_background_guass_withlocalization;
        y = y(first_frame:last_frame);
        %y(y==0) = nan;

        % Find the length of the cell trace
        data_width = last_frame-first_frame;

        %Insert the current cell data into a combined variable 1:length cell
            %trace
        comb_y(i,1:data_width+1) = y';
    end

    % remove data once only 2 cells trace remain in time
    % non_zero_comb_y = comb_y > 0;
    non_zero_comb_y = ~isnan(comb_y); % 20240822 zeros are valid data points
    cells_at_diff_points = sum(non_zero_comb_y,1);
    valid_idx = cells_at_diff_points > 2;

    hold on
        x = Stats_ROIs{ROI_with_TimeStamps,1}.Timing_sec(1:size(comb_y,2)) ./ 60;
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