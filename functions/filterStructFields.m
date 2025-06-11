function filteredStruct = filterStructFields(inputStruct, validIndices)
    % FILTERSTRUCTFIELDS Filters fields of a structure based on valid indices.
    %
    %   Description: 
    %   filters the fields of the input structure `inputStruct` using the specified
    %   `validIndices`. It updates the 'ID' field and all fields within
    %   'timelock.HEP' to only include the elements at the positions specified
    %   by `validIndices`.
    %
    %   Inputs:
    %       inputStruct  - A structure containing an 'ID' field and a nested
    %                      structure 'timelock.HEP' with fields to be filtered.
    %       validIndices - A logical or numerical array indicating which indices
    %                      to retain in the filtering process.
    %
    %   Outputs:
    %       filteredStruct - The resulting structure with fields filtered
    %                        according to `validIndices`.
    %

    % Filter the 'ID' field
    inputStruct.ID  = inputStruct.ID(validIndices);

    % Get all field names within 'timelock.HEP'
    fieldNames = fieldnames(inputStruct.timelock.HEP); 
    for iField = 1:length(fieldNames)
        currentField = fieldNames{iField};
        % Filter each field in 'timelock.HEP' using validIndices
        inputStruct.timelock.HEP.(currentField) = inputStruct.timelock.HEP.(currentField)(validIndices);
    end
    % Return the filtered structure
    filteredStruct = inputStruct;
end
