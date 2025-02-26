function updatedStruct = KLS_updateStruct(originalStruct, newStruct)
    % ChatGPT 4o code from 20240910, check for funcationality 20240910 KLS
    % 

    % Get the field names of the new structure
    newFields = fieldnames(newStruct);
    
    % Loop over the fields of the new structure and update the original structure
    for i = 1:length(newFields)
        field = newFields{i};

        originalStruct.(field) = newStruct.(field);
    end
    
    % Return the updated structure
    updatedStruct = originalStruct;
end