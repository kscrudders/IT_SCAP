function Image = KLS_TifImportAll(FileName)
    % Imports an tif file into a Matlab Matrix
    bfcell = bfopen(FileName);
    H = size(bfcell{1,1}{1,1},1);
    W = size(bfcell{1,1}{1,1},2);
    Z = size(bfcell{1,1},1);
    
    Image = zeros(H,W,Z);
    i = 1;
    while i <= Z
        Image(:,:,i) = bfcell{1,1}{i,1};
        i = i+1;
    end
end