function IRCE_Plotting_SC_ResponseSignal_sorted(Stats_ROIs, fig_color, inclusion_idx, IRM_idx_to_time_idx)
    if isempty(inclusion_idx)
        inclusion_idx = 1:length(Stats_ROIs);
    end
    if isempty(IRM_idx_to_time_idx)
        IRM_idx_to_time_idx = 1:length(Stats_ROIs{1,1}.integrated_response_above_background_guass);
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
    
    num_ROIs = length(Stats_ROIs);

    num_frames = 0;
    for i = 1:length(Stats_ROIs)
        if length(Stats_ROIs{i,1}.integrated_response_above_background_guass) > num_frames
            ROI_with_TimeStamps = i;
        end        
        num_frames = max([num_frames length(Stats_ROIs{i,1}.integrated_response_above_background_guass)]);
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
        
        first_frame = find(~isnan(Stats_ROIs{i,1}.Area),1,'first');
        last_frame = find(~isnan(Stats_ROIs{i,1}.Area),1,'last');

        first_frame = IRM_idx_to_time_idx(first_frame);
        last_frame = IRM_idx_to_time_idx(last_frame);
        
        x = Stats_ROIs{i,1}.Timing_sec;
        x = x(first_frame:last_frame);
        x = x-x(1); % Reorder starting at 0
        x = x./60; % convert frames to minutes
        
        y = Stats_ROIs{i,1}.integrated_response_above_background_guass;
        y = y(first_frame:last_frame);
        hold on
            plot(x,y, 'Color', shades(n,:));
        hold off
    end
end