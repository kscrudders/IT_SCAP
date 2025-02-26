function out = KLS_NormStack(raw_data)
    Z = size(raw_data,3); % # of frames
    i = 1;
    while i <= Z
       raw_data(:,:,i) = raw_data(:,:,i)./max(squeeze(raw_data(:,:,i)),[],'all');
       i = i+1;
    end
    out = raw_data;
end