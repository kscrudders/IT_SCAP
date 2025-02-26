function KLS_Check_tif_Imwrite()
    try
        imwrite(uint8(randi([0, 255], 256, 256)),'sample.tiff','tiff');
        
        if isfile('sample.tiff')
            delete('sample.tiff');
        end
    catch
        disp('Imwrite is not working')
    end
    
end