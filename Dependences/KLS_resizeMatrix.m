function newMatrix = KLS_resizeMatrix(oldMatrix, maxT)
    x = size(oldMatrix, 1);
    y = size(oldMatrix, 2);
    t = size(oldMatrix, 3);
    newMatrix = zeros(x, y, maxT); % Pre-allocate new matrix with the maximum size
    
    for i = 1:t
        % Repeat each slice to fill the new matrix
        repeatFactor = ceil(maxT / t); % Calculate how many times each slice should be repeated
        startIdx = (i-1)*repeatFactor + 1;
        endIdx = min(i*repeatFactor, maxT);
        newMatrix(:, :, startIdx:endIdx) = repmat(oldMatrix(:, :, i), [1, 1, endIdx - startIdx + 1]);
    end
end