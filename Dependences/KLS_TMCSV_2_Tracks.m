function out = KLS_TMCSV_2_Tracks(All_Spots_Stats_CSV, length_t)
    if nargin < 2
        t = max((All_Spots_Stats_CSV(:,9)))+1; % number of time points
    else
        t = length_t; % number of time points
    end
    
    % If the TrackMate column labels were imported, remove them
    if find(all(isnan(All_Spots_Stats_CSV(1,:))))
        All_Spots_Stats_CSV = All_Spots_Stats_CSV(2:end,:);
    end

    N = max((All_Spots_Stats_CSV(:,3)))+1; % number of tracks

    if isnan(N) % if there are only localizations, no tracks
        N = sum(~isnan(All_Spots_Stats_CSV(:,2))); % number of localizations
        Tracks = NaN([N t 2]); % empty Structure
    
        idx = find(~isnan(All_Spots_Stats_CSV(:,2)));
        for i = 1:length(idx)
            z = idx(i);
            Tracks(z, All_Spots_Stats_CSV(z,9)+1, 1) = All_Spots_Stats_CSV(z,5); % set x
            Tracks(z, All_Spots_Stats_CSV(z,9)+1, 2) = All_Spots_Stats_CSV(z,6); % set y
        end

    else % if there are tracks
        % if there are tracks and isolated localizations
        if any(isnan(All_Spots_Stats_CSV(:,3)))
            num_isolated_localizations = sum(isnan(All_Spots_Stats_CSV(:,3)));

            Tracks = NaN([N+num_isolated_localizations t 2]); % empty Structure
        
            % Bring over the localizations in tracks first
            i = 0;
            while i <= N
                 idx = find((All_Spots_Stats_CSV(:,3))==i);
                 L = size(idx,1);
                 ii = 1;
                 while ii <= L
                      Tracks(i+1, All_Spots_Stats_CSV(idx(ii),9)+1, 1) = (All_Spots_Stats_CSV(idx(ii),5)); % set x
                      Tracks(i+1, All_Spots_Stats_CSV(idx(ii),9)+1, 2) = (All_Spots_Stats_CSV(idx(ii),6)); % set y
                      ii = ii+1;
                 end
                 i = i+1;
            end

            % Bring over the localizations without tracks
            idx = find(isnan(All_Spots_Stats_CSV(:,3)));
            for i = 1:length(idx)
                z = idx(i);
                 Tracks(i+N, All_Spots_Stats_CSV(z,9)+1, 1) = All_Spots_Stats_CSV(z,5); % set x
                 Tracks(i+N, All_Spots_Stats_CSV(z,9)+1, 2) = All_Spots_Stats_CSV(z,6); % set y
            end
        else % if there are only tracks
            Tracks = NaN([N t 2]); % empty Structure
        
            i = 0;
            while i <= N
                 idx = find((All_Spots_Stats_CSV(:,3))==i);
                 L = size(idx,1);
                 ii = 1;
                 while ii <= L
                      Tracks(i+1, All_Spots_Stats_CSV(idx(ii),9)+1, 1) = (All_Spots_Stats_CSV(idx(ii),5)); % set x
                      Tracks(i+1, All_Spots_Stats_CSV(idx(ii),9)+1, 2) = (All_Spots_Stats_CSV(idx(ii),6)); % set y
                      ii = ii+1;
                 end
                 i = i+1;
            end
        end
    end

    out = Tracks ./ 0.157; % Convert from micron to px
end