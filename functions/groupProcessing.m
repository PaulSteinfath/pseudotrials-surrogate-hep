function allsubj = groupProcessing(allsubj)
% GROUPPROCESSING Corrects pre-target heartbeat data by subtracting random pseudo trials.
%
%   Description:
%   takes a structure `allsubj` containing
%   subject data and performs a correction on the pre-target heartbeat data.
%   The correction involves subtracting the average of random pseudo heartbeat
%   trials from the actual pre-target heartbeat trials for each subject.
%
%   Input:
%     allsubj - A structure with fields:
%       - timelock.HEP.pre_target_heartbeat: Cell array of structures with
%         average pre-target heartbeat data for each subject.
%       - timelock.HEP.pre_target_random_pseudo_heartbeat: Cell array of
%         structures with average random pseudo heartbeat data for each subject.
%
%   Output:
%     allsubj - The input structure with an additional field:
%       - timelock.HEP.pre_target_heartbeat_correct: Cell array of structures
%         containing the corrected pre-target heartbeat data.
%
%   Note:
%     - If the subtraction fails for a subject, the corresponding entry in
%       `pre_target_heartbeat_correct` will be left empty.
%     - Ensure that the input structure contains the necessary fields and data
%       before calling this function.
%

    %% Individual subject subtract random pseudo trials
    allsubj_corr = allsubj;
    
    % Pre-target heartbeat corrected with random pseudo heartbeat
    if isfield(allsubj.timelock.HEP, 'pre_target_heartbeat')
        % Create dummy structure for the correction output
        allsubj_corr.timelock.HEP.pre_target_heartbeat_correct = allsubj_corr.timelock.HEP.pre_target_heartbeat;
        
        for i = 1:length(allsubj.timelock.HEP.pre_target_heartbeat)
            try
                allsubj_corr.timelock.HEP.pre_target_heartbeat_correct{i}.avg = ...
                    allsubj.timelock.HEP.pre_target_heartbeat{i}.avg - ...
                    allsubj.timelock.HEP.pre_target_random_pseudo_heartbeat{i}.avg;
            catch
                % Keep empty if subtraction fails
            end
        end
        
        % Add corrected data to output structure
        allsubj.timelock.HEP.pre_target_heartbeat_correct = allsubj_corr.timelock.HEP.pre_target_heartbeat_correct;
    end

end