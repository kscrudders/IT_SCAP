function out = KLS_Fill_TrackZeros(tracks)
    if sum(isnan(tracks(:))) > 0
        tracks(isnan(tracks)) = 0;
        back2nan = 1;
    else
        back2nan = 0;
    end
    
    
    out = tracks;
    
    i = 1;
    while i <= size(tracks,1)
        [~,t_i,~] = find(tracks(i,:,1),1,'first'); % First time point
        [~,t_l,~] = find(tracks(i,:,1),1,'last'); % Last time point
        
        t_all = find(tracks(i,t_i:t_l,1)); % idx of time points between start and end of track
        num_gaps = find(~tracks(i,t_i:t_l,1)); % idx of gaps between start and end of track
        
        if num_gaps > 0 % if there are gaps
            ii = 1;
            while ii <= size(num_gaps,2)
                last_time = (num_gaps(ii)-1)+t_i-1; % last non-zero time point
                next_time_idx = find(t_all > last_time,1,'first'); % idx of next timepoint
                next_time = t_all(next_time_idx)+t_i-1; % next timepoint
                
                gap_size = next_time - last_time -1;
                
                out(i,num_gaps(ii)+t_i-1,1) = out(i,last_time,1) + ...
                    (out(i,next_time,1) - out(i,last_time,1))/(gap_size+1);
                out(i,num_gaps(ii)+t_i-1,2) = out(i,last_time,2) + ...
                    (out(i,next_time,2) - out(i,last_time,2))/(gap_size+1);
                ii=ii+1;
            end
        end
        i = i+1;
    end
    
    if back2nan == 1
        out(out == 0) = NaN;
    end
end