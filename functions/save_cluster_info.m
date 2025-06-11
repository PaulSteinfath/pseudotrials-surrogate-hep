function save_cluster_info(stat, name, save_data)
    % SAVE_CLUSTER_INFO - Identifies and saves significant clusters from statistical data.
    %
    % Inputs:
    %   stat - A structure containing statistical data with fields such as:
    %          - label: Cell array of channel labels.
    %          - posclusters: Array of structures for positive clusters, each with fields:
    %            - prob: Probability value of the cluster.
    %            - clusterstat: T-value of the cluster.
    %          - posclusterslabelmat: Matrix indicating channel and time indices for positive clusters.
    %          - negclusters: Array of structures for negative clusters, similar to posclusters.
    %          - negclusterslabelmat: Matrix indicating channel and time indices for negative clusters.
    %          - time: Array of time points.
    %   name - A string representing the name of the comparison or dataset.
    %   save_data - A string specifying the directory path where the CSV file should be saved.
    %
    % Outputs:
    %   A CSV file named '[name]_significant_clusters.csv' is saved in the specified directory
    %   if significant clusters are found. The file contains columns:
    %   - Comparison: Name of the comparison or dataset.
    %   - Direction: 'positive' or 'negative' indicating the cluster type.
    %   - T_Value: T-value of the cluster.
    %   - P_Value: P-value of the cluster.
    %   - Channels: Comma-separated list of significant channels.
    %   - Time_Range: Time range of the significant cluster.

    significant_clusters = {};
    
    % Get the list of channel labels
    channel_labels = stat.label;
    
    % Process positive clusters
    if isfield(stat, 'posclusters') && ~isempty(stat.posclusters)
        for i = 1:length(stat.posclusters)
            if stat.posclusters(i).prob < 0.05
                % Find the significant channels
                sig_channels_idx = find(any(stat.posclusterslabelmat == i, 2));
                sig_channels = channel_labels(sig_channels_idx);
                
                % Get the time range
                time_range = [stat.time(find(any(stat.posclusterslabelmat == i, 1), 1, 'first')), ...
                              stat.time(find(any(stat.posclusterslabelmat == i, 1), 1, 'last'))];
                
                % Format time range as a single cell
                time_range_str = sprintf('[%.4f, %.4f]', time_range(1), time_range(2));
                
                % Format T_Value and P_Value
                t_value = sprintf('%.2f', stat.posclusters(i).clusterstat);
                if stat.posclusters(i).prob < 0.001
                    p_value = '<= 0.001';
                else
                    p_value = sprintf('%.3f', stat.posclusters(i).prob);
                end
                
                % Add to the significant clusters list
                significant_clusters = [significant_clusters; ...
                    {name, 'positive', t_value, p_value, ...
                    strjoin(sig_channels', ', '), time_range_str}];
            end
        end
    end

    % Process negative clusters
    if isfield(stat, 'negclusters') && ~isempty(stat.negclusters)
        for i = 1:length(stat.negclusters)
            if stat.negclusters(i).prob < 0.05
                % Find the significant channels
                sig_channels_idx = find(any(stat.negclusterslabelmat == i, 2));
                sig_channels = channel_labels(sig_channels_idx);
                
                % Get the time range
                time_range = [stat.time(find(any(stat.negclusterslabelmat == i, 1), 1, 'first')), ...
                              stat.time(find(any(stat.negclusterslabelmat == i, 1), 1, 'last'))];
                
                % Format time range as a single cell
                time_range_str = sprintf('[%.4f, %.4f]', time_range(1), time_range(2));
                
                % Format T_Value and P_Value
                t_value = sprintf('%.2f', stat.negclusters(i).clusterstat);
                if stat.negclusters(i).prob < 0.001
                    p_value = '<= 0.001';
                else
                    p_value = sprintf('%.3f', stat.negclusters(i).prob);
                end
                
                % Add to the significant clusters list
                significant_clusters = [significant_clusters; ...
                    {name, 'negative', t_value, p_value, ...
                    strjoin(sig_channels', ', '), time_range_str}];
            end
        end
    end

    % Save the significant clusters to a CSV file
    if ~isempty(significant_clusters)
        T = cell2table(significant_clusters, 'VariableNames', {'Comparison', 'Direction', 'T_Value', 'P_Value', 'Channels', 'Time_Range'});
        writetable(T, fullfile(save_data, [name, '_significant_clusters.csv']));
    end
end
