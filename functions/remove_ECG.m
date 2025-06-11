function [ECG_struct, data_struct] = remove_ECG(data_struct, conditions, condition_type)
    % REMOVE_ECG Extracts and removes ECG data from a data structure.
    % 
    %   [ECG_STRUCT, DATA_STRUCT] = REMOVE_ECG(DATA_STRUCT, CONDITIONS, CONDITION_TYPE)
    %   extracts the ECG data from the input DATA_STRUCT for the specified
    %   CONDITIONS and CONDITION_TYPE, and returns it in ECG_STRUCT. The ECG
    %   data is removed from DATA_STRUCT.
    %
    %   Inputs:
    %     DATA_STRUCT    - A structure containing timelock and optionally
    %                      grand average data with ECG data included.
    %     CONDITIONS     - A cell array of condition names to process.
    %     CONDITION_TYPE - A string specifying the type of condition to
    %                      process (e.g., 'pre', 'post').
    %
    %   Outputs:
    %     ECG_STRUCT     - A structure containing only the ECG data extracted
    %                      from DATA_STRUCT.
    %     DATA_STRUCT    - The input data structure with ECG data removed.

    ECG_struct = data_struct;
    
    for cond = 1:length(conditions)
        num_cells = size(data_struct.timelock.(condition_type).(conditions{cond}), 2);
    
        % get ECG index
        ECG_idx =  find(strcmp(data_struct.timelock.(condition_type).(conditions{cond}){1}.label, 'ECG'));
    
        for i = 1:num_cells
            ECG_struct.timelock.(condition_type).(conditions{cond}){i}.avg = data_struct.timelock.(condition_type).(conditions{cond}){i}.avg(ECG_idx, :);
            ECG_struct.timelock.(condition_type).(conditions{cond}){i}.var = data_struct.timelock.(condition_type).(conditions{cond}){i}.var(ECG_idx, :);
            ECG_struct.timelock.(condition_type).(conditions{cond}){i}.dof = data_struct.timelock.(condition_type).(conditions{cond}){i}.dof(ECG_idx, :);
            ECG_struct.timelock.(condition_type).(conditions{cond}){i}.label = {'ECG'};
    
            data_struct.timelock.(condition_type).(conditions{cond}){i}.avg(ECG_idx, :) = [];
            data_struct.timelock.(condition_type).(conditions{cond}){i}.var(ECG_idx, :) = [];
            data_struct.timelock.(condition_type).(conditions{cond}){i}.dof(ECG_idx, :) = [];
            data_struct.timelock.(condition_type).(conditions{cond}){i }.label(ECG_idx, :) = [];
        end
    
        % copy ECG from grandavg
        % ECG_struct.grandAvg.(condition_type).(conditions{cond}).avg = []
        if isfield( ECG_struct, 'grandAvg')
            ECG_struct.grandAvg.(condition_type).(conditions{cond}).avg = data_struct.grandAvg.(condition_type).(conditions{cond}).avg(ECG_idx, :);
            ECG_struct.grandAvg.(condition_type).(conditions{cond}).var = data_struct.grandAvg.(condition_type).(conditions{cond}).var(ECG_idx, :);
            ECG_struct.grandAvg.(condition_type).(conditions{cond}).dof = data_struct.grandAvg.(condition_type).(conditions{cond}).dof(ECG_idx, :);
            ECG_struct.grandAvg.(condition_type).(conditions{cond}).label = {'ECG'};
    
            %remove ECG from grandavg
            data_struct.grandAvg.(condition_type).(conditions{cond}).avg(ECG_idx, :) =  [];
            data_struct.grandAvg.(condition_type).(conditions{cond}).var(ECG_idx, :) =  [];
            data_struct.grandAvg.(condition_type).(conditions{cond}).dof(ECG_idx, :) =  [];
        end
    end
end
