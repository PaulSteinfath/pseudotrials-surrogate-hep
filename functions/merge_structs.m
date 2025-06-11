function merged = merge_structs(defaults, custom)
    % MERGE_STRUCTS Helper function to merge two structures, with custom values overriding defaults.
    %
    %   merged = MERGE_STRUCTS(defaults, custom) takes two input structures,
    %   'defaults' and 'custom', and returns a new structure 'merged'. The 
    %   'custom' structure's fields will override those in 'defaults'. If a 
    %   field in 'custom' is also a structure, the function will recursively 
    %   merge these substructures.
    %
    %   Inputs:
    %       defaults - A structure containing default values.
    %       custom   - A structure containing custom values that override the defaults.
    %
    %   Output:
    %       merged   - A structure containing the merged result of 'defaults' 
    %                  and 'custom'.

    merged = defaults;
    if ~isempty(custom)
        fields = fieldnames(custom);
        for i = 1:length(fields)
            if isstruct(custom.(fields{i})) && isfield(merged, fields{i})
                merged.(fields{i}) = merge_structs(merged.(fields{i}), custom.(fields{i}));
            else
                merged.(fields{i}) = custom.(fields{i});
            end
        end
    end
end