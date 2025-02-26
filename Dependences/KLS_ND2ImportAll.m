function Image = KLS_ND2ImportAll(FileAddress)
% Imports an ND2 file into a Matlab Matrix
% Inputs (1): Full File Address, or file name if cd is the file's folder
% Outputs (1): 2D or 3D stack of images as organized in the nd2 file

    bfcell = bfopen(FileAddress);
    H = size(bfcell{1,1}{1,1},1);
    W = size(bfcell{1,1}{1,1},2);
    
    Z_All= size(bfcell,1); % How many data cells are there?
    Z_01 = size(bfcell{1,1},1); % How big is z in the first cell
    
    Z = Z_01*(Z_01>=Z_All)+Z_All*(Z_All>Z_01); 
    
    Image = zeros(H,W,Z);
    i = 1;
    if Z_01>Z_All
        while i <= Z
            Image(:,:,i) = bfcell{1,1}{i,1};
            i = i+1;
        end
    else
        while i <= Z
            Image(:,:,i) = bfcell{i,1}{1,1};
            i = i+1;
        end
    end
end