function is_overlapping = IRCE_check_overlap(line1, line2)
    % Define a small tolerance value for overlap
    tolerance = 5/.157; % 2 micron tolerance
    
    % Check if the lines share any common points within the tolerance
    is_overlapping = any(pdist2(line1, line2) < tolerance, 'all');
end