function data_in = rem_empty(data_in, condition, IDs)
    % REM_EMPTY Remove empty cells from a specified field in a structure and update IDs.
    %
    %   Removes empty cells from the field specified by 'condition' 
    %   in the input structure 'data_in'. It also updates the
    %   corresponding IDs to maintain alignment with the non-empty entries.
    %
    %   Inputs:
    %       data_in   - A structure containing data fields.
    %       condition - A string or cell array specifying the field name to check for empty cells.
    %       IDs       - An array of IDs corresponding to the entries in the specified field.
    %
    %   Outputs:
    %       data_in   - The updated structure with empty cells removed from the specified field
    %                   and corresponding IDs updated.

    % if condition is a cell, convert it to a character array
    if iscell(condition)
        condition_name = char(condition);
    elseif ischar(condition)
        condition_name = condition;
    end
    
    %% remove empty cells
    empty_cells = cell2mat(cellfun(@(x) isempty(x), data_in.(condition_name), 'UniformOutput', false));
    data_in.(condition_name) = data_in.(condition_name)(~empty_cells);
    
    %remove correspoding IDs
    if ~isempty(IDs)
        IDs = IDs(~empty_cells);
    end
    
    %save updated IDs
    data_in.IDs = IDs;

end