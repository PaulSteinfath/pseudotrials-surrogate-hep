function [EEG, cV, rejV] = ecg_ica_corr(EEG, ECG_template, time_window, R_marker, SDcrit)
% ECG_ICA_CORR - Identifies ECG-related ICA components by correlation.
%
% Syntax:  [EEG, cV, rejV] = ecg_ica_corr(EEG, ECG_template, time_window, R_marker, SDcrit)
%
% Inputs:
%    EEG - EEG structure containing ICA activations.
%    ECG_template - ECG template for correlation (optional).
%    time_window - Time window for epoching around R-peak markers.
%    R_marker - Marker for R-peaks in the ECG signal.
%    SDcrit - Standard deviation criterion for identifying ECG components.
%
% Outputs:
%    EEG - EEG structure with identified ECG components removed.
%    cV - Correlation values of ICA components with ECG template.
%    rejV - Logical vector indicating rejected components.
%
% Description:
%    This function identifies ECG-related Independent Component Analysis (ICA)
%    components by correlating them with an ECG template. It performs the following steps:
%    1) Epochs the ICA activations based on R-peak markers.
%    2) Calculates the correlation of each average ICA waveform to the ECG template.
%    3) Identifies components above a certain standard deviation threshold as ECG artifacts.
%    4) Checks for oscillatory peaks to preserve potential brain-like activity.
%
% Example:
%    [EEG, cV, rejV] = ecg_ica_corr(EEG, [], [-0.2 0.6], 'R', 1.5);
%
% See also: pop_epoch, pop_rmbase, Compute_RP

% Create a temporary EEG structure for ICA data
EEGtemp = EEG;
EEGtemp.data = EEG.icaact;

% Check if the data is 2D or less
if ndims(EEGtemp.data) <=2
    % Calculate the number of additional channels needed
    add_chans =  size(EEG.data,1) -size(EEGtemp.data,1);

    % exchange EEG data with ICA data for easier epoching
    if size(EEGtemp.data,1) ~= size(EEG.data,1)
        EEGtemp.data = [EEGtemp.data; zeros(add_chans, length(EEGtemp.data))];
    end

    % Epoch ICA data based on R-peak markers
    EEGtemp = pop_epoch(EEGtemp, {R_marker} , time_window);

    % Remove the added channels
    EEGtemp.data(end-add_chans:end,:,:) = [];

    % Epoch EEG according to ECG marker
    EEG = pop_epoch(EEG, {R_marker} , time_window);
end

% Remove baseline from EEG data
EEG = pop_rmbase(EEG, [time_window(1) 0]); 

% Calculate the average ICA component
avgIC = mean(EEGtemp.data,3);

% If no ECG template is provided, try to create one from the ECG channel
if isempty(ECG_template)
    try
        ECG_template = mean(EEG.data(strcmp({EEG.chanlocs.labels},'ECG'),:,:),3);
    catch
    end
end

% Calculate correlation values with the ECG template
cV  = abs(corr(avgIC',ECG_template'));
% Determine the correlation threshold
corthreshV = mean(cV)+ SDcrit * std(cV); 
% Identify components that exceed the threshold
rejV = cV > corthreshV; 
% Find the indices of rejected components
comps = find(rejV);

% If no components are identified, display a message
if sum(comps) == 0
    disp('No ECG components identified')
else
    
    % Check for alpha relative power to preserve potential brain-like activity
    for k=1:sum(rejV)
        
        comp = comps(k);
        alpha_rel_pow =  Compute_RP(EEGtemp.data(comp,:), 500,[1,45], [7 15]);

        % if peak greater than 50% - keep the component
        if alpha_rel_pow > 0.5  
            rejV(comps(k)) = 0;
        else
        end
        
    end
    
end



end