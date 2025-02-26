function label_out = KLS_Label_fullStackTracking(label_in)
    % consideration from BoT-SORT: Robust Associations Multi-Pedestrian Tracking
    % Nir Aharon, Roy Orfaig, Ben-Zion Bobrovsky
    % https://doi.org/10.48550/arXiv.2206.14651

    % Convert input to binary image if necessary
    bw_in = label_in > 0;

    % Initialize output label image
    label_out = zeros(size(label_in));

    % Initialize tracks (as a structure array)
    tracks = initializeTracks();

    nextId = 1; % ID of the next track

    % Loop through frames
    for i = 1:size(bw_in, 3)
        frame = bw_in(:,:,i);

        % Detect regions in the current frame
        [current_labels, ~] = bwlabeln(frame, 4);

        % Remove labels touching the edge of the FOV
        edgeLabels = unique([current_labels(1,:), current_labels(end,:), current_labels(:,1)', current_labels(:,end)']);
        edgeLabels(edgeLabels == 0) = [];    

        for k = 1:length(edgeLabels)
            current_labels(current_labels == edgeLabels(k)) = 0;
        end
        
        [current_labels, ~] = bwlabeln(current_labels > 0, 4);

        % Extract properties of detected regions
        detections = regionprops(current_labels, 'Centroid', 'Area', 'BoundingBox', 'PixelIdxList', 'MajorAxisLength');

        % Predict new locations of existing tracks
        tracks = predictNewLocationsOfTracks(tracks);

        % Extract appearance features for detections
        detections = extractAppearanceFeatures(detections, label_in(:,:,i));

        % Assign detections to tracks
        [assignments, unassignedTracks, unassignedDetections] = detectionToTrackAssignment(tracks, detections);

        % Update assigned tracks
        tracks = updateAssignedTracks(tracks, detections, assignments);

        % Update unassigned tracks
        tracks = updateUnassignedTracks(tracks, unassignedTracks);

        % Delete lost tracks
        tracks = deleteLostTracks(tracks);

        % Create new tracks
        [tracks, nextId] = createNewTracks(tracks, detections, unassignedDetections, nextId);

        % Update label_out with current tracks
        label_out(:,:,i) = createLabelImageFromTracks(tracks, size(frame));
    end
end

% Initialize tracks structure
function tracks = initializeTracks()
    tracks = struct(...
        'id', {}, ...
        'bbox', {}, ...
        'kalmanFilter', {}, ...
        'age', {}, ...
        'totalVisibleCount', {}, ...
        'consecutiveInvisibleCount', {}, ...
        'centroid', {}, ...
        'appearance', {}, ...
        'PixelIdxList', {}, ...
        'maxLength', {});
end

% Predict new locations of tracks using Kalman filter
function tracks = predictNewLocationsOfTracks(tracks)
    for i = 1:length(tracks)
        if ~isempty(tracks(i).kalmanFilter)
            % Predict the current location of the track
            predictedCentroid = predict(tracks(i).kalmanFilter);

            % Check for NaN in the predicted centroid
            if any(isnan(predictedCentroid))
                % Handle the NaN case
                tracks(i).centroid = tracks(i).centroid; % Keep the last known centroid
            else
                % Update the bounding box
                bbox = tracks(i).bbox;
                predictedCentroid = predictedCentroid - bbox(3:4) / 2;
                tracks(i).bbox(1:2) = predictedCentroid;

                % Update the centroid
                tracks(i).centroid = predictedCentroid + bbox(3:4) / 2;
            end
        end
    end
end

% Extract appearance features for detections
function detections = extractAppearanceFeatures(detections, frame)
    for i = 1:length(detections)
        % Extract the pixel values within the detection
        mask = false(size(frame));
        mask(detections(i).PixelIdxList) = true;
        % Compute appearance features (e.g., histogram)
        detections(i).appearance = imhist(frame(mask));
    end
end

% Assign detections to tracks
function [assignments, unassignedTracks, unassignedDetections] = detectionToTrackAssignment(tracks, detections)
    nTracks = length(tracks);
    nDetections = length(detections);

    % Initialize outputs
    assignments = [];
    unassignedTracks = [];
    unassignedDetections = [];

    if nTracks == 0
        % If there are no tracks, all detections are unassigned
        unassignedDetections = 1:nDetections;
        return;
    end

    if nDetections == 0
        % If there are no detections, all tracks are unassigned
        unassignedTracks = 1:nTracks;
        return;
    end

    % Compute cost matrix
    cost = zeros(nTracks, nDetections);
    for i = 1:nTracks
        for j = 1:nDetections
            % Check if centroids are valid
            if any(isnan(tracks(i).centroid)) || any(isnan(detections(j).Centroid))
                cost(i,j) = Inf; % Assign a high cost to invalid pairs
                continue;
            end

            % Compute the Euclidean distance between centroids
            distance = norm(tracks(i).centroid - detections(j).Centroid);

            % Get maximum allowed distance
            maxDistance = 1.5 * tracks(i).maxLength;
            
            % Handle missing or zero maxDistance values. Likely unneeded:
            % In detectionToTrackAssignment, adjust maxDistance computation
            if isempty(maxDistance) || maxDistance == 0
                maxDistance = 20; % e.g., 20
            end

            if distance > maxDistance
                cost(i,j) = Inf; % Assign high cost to pairs that are too far apart
                continue;
            end

            % Compute the cost using Euclidean distance between centroids
            spatialCost = distance;

            % Compute appearance cost
            appearanceCost = compareAppearance(tracks(i).appearance, detections(j).appearance);

            % Combine costs (e.g., weighted sum)
            cost(i,j) = spatialCost + appearanceCost;
        end
    end

    % Replace any NaN or Inf in cost matrix with a high cost
    cost(isnan(cost) | isinf(cost)) = 1e5;

    % Solve the assignment problem
    % Values between 10 - 50 appear to work well.
    % Smaller values make linking less likely
    % Larger values make linking more likely
    costOfNonAssignment = 50; % Set a threshold
    [assignments, unassignedTracks, unassignedDetections] = assignDetectionsToTracks(cost, costOfNonAssignment);
end


% Update assigned tracks
function tracks = updateAssignedTracks(tracks, detections, assignments)
    numAssignedTracks = size(assignments, 1);
    for i = 1:numAssignedTracks
        trackIdx = assignments(i, 1);
        detectionIdx = assignments(i, 2);

        % Correct the estimate of the object's location using the new detection
        centroid = detections(detectionIdx).Centroid;

        if any(isnan(centroid))
            continue; % Skip updating if centroid is invalid
        end

        correct(tracks(trackIdx).kalmanFilter, centroid);

        % Update track's properties
        tracks(trackIdx).bbox = detections(detectionIdx).BoundingBox;
        tracks(trackIdx).centroid = centroid;
        tracks(trackIdx).appearance = detections(detectionIdx).appearance;
        tracks(trackIdx).PixelIdxList = detections(detectionIdx).PixelIdxList;
        tracks(trackIdx).maxLength = detections(detectionIdx).MajorAxisLength;

        % Update visibility
        tracks(trackIdx).age = tracks(trackIdx).age + 1;
        tracks(trackIdx).totalVisibleCount = tracks(trackIdx).totalVisibleCount + 1;
        tracks(trackIdx).consecutiveInvisibleCount = 0;
    end
end



% Update unassigned tracks
function tracks = updateUnassignedTracks(tracks, unassignedTracks)
    for i = 1:length(unassignedTracks)
        ind = unassignedTracks(i);
        tracks(ind).age = tracks(ind).age + 1;
        tracks(ind).consecutiveInvisibleCount = tracks(ind).consecutiveInvisibleCount + 1;
    end
end

% Delete lost tracks
function tracks = deleteLostTracks(tracks)
    if isempty(tracks)
        return;
    end

    invisibleForTooLong = 5;
    ageThreshold = 8;

    % Compute the fraction of the track's age for which it was visible
    ages = [tracks(:).age];
    totalVisibleCounts = [tracks(:).totalVisibleCount];
    visibility = totalVisibleCounts ./ ages;

    % Find the indices of 'lost' tracks
    lostInds = (ages < ageThreshold & visibility < 0.6) | ...
        [tracks(:).consecutiveInvisibleCount] >= invisibleForTooLong;

    % Delete lost tracks
    tracks = tracks(~lostInds);
end

% Create new tracks
function [tracks, nextId] = createNewTracks(tracks, detections, unassignedDetections, nextId)
    for i = 1:length(unassignedDetections)
        detectionIdx = unassignedDetections(i);
        centroid = detections(detectionIdx).Centroid;
        bbox = detections(detectionIdx).BoundingBox;

        % Create a Kalman filter object
        kalmanFilter = configureKalmanFilter('ConstantVelocity', centroid, [200, 50], [100, 25], 100);

        % Create new track
        newTrack = struct(...
            'id', nextId, ...
            'bbox', bbox, ...
            'kalmanFilter', kalmanFilter, ...
            'age', 1, ...
            'totalVisibleCount', 1, ...
            'consecutiveInvisibleCount', 0, ...
            'centroid', centroid, ...
            'appearance', detections(detectionIdx).appearance, ...
            'PixelIdxList', detections(detectionIdx).PixelIdxList, ...
            'maxLength', detections(detectionIdx).MajorAxisLength);

        % Add to tracks
        tracks(end + 1) = newTrack;

        % Increment the nextId
        nextId = nextId + 1;
    end
end


% Create label image from tracks
function labelImg = createLabelImageFromTracks(tracks, frameSize)
    labelImg = zeros(frameSize);
    for i = 1:length(tracks)
        if tracks(i).consecutiveInvisibleCount == 0  % Only include visible tracks
            % Assign the track ID to the pixels in PixelIdxList
            labelImg(tracks(i).PixelIdxList) = tracks(i).id;
        end
    end
end

% Compare appearance features (e.g., histogram intersection)
function cost = compareAppearance(appearance1, appearance2)
    % Check for empty or invalid appearances
    if isempty(appearance1) || isempty(appearance2) || any(isnan(appearance1)) || any(isnan(appearance2))
        cost = 1; % Maximum cost
        return;
    end

    total = sum(appearance1);
    if total == 0
        cost = 1; % Assign maximum cost to avoid division by zero
        return;
    end

    intersection = sum(min(appearance1, appearance2));
    cost = 1 - (intersection / total); % Normalize to [0,1]
end
