function [cfg, yminmax] = plot_sig_topo(stat, cfg, data, time_array, clust_num, direction, inc_title)
    % PLOT_SIG_TOPO - Plot significant topographical data from EEG/MEG analysis.
    %
    % Description:
    %   This function visualizes significant clusters from EEG/MEG statistical analysis.
    %   It adjusts the data labels to match the statistical labels, configures the plot
    %   appearance, and highlights significant channels and time points for both positive
    %   and negative clusters. The function also customizes the colorbar and optionally
    %   includes a title with cluster information.
    %
    % Inputs:
    %   stat       - Statistical data structure containing cluster information.
    %   cfg        - Configuration structure for plotting.
    %   data       - Data structure containing EEG/MEG data.
    %   time_array - Array of time points.
    %   clust_num  - Cluster number to be plotted.
    %   direction  - Direction of the cluster ('pos' for positive, 'neg' for negative).
    %   inc_title  - Boolean to include a title in the plot (1 to include, 0 to exclude).
    %
    % Outputs:
    %   cfg        - Updated configuration structure.
    %   yminmax    - Minimum and maximum values of the data in the cluster.

    %

    cfg.comment = 'no';

    % Match size if ECG was still present
    if size(stat.prob,1) ~= size(data.label)
        matching_labels = intersect(stat.label, data.label);
        matching_indices = find(ismember(data.label, matching_labels));
        data.label = data.label(matching_indices);
        data.elec.chanpos = data.elec.chanpos(matching_indices, :);
        data.elec.chantype = data.elec.chantype(matching_indices);
        data.elec.chanunit = data.elec.chanunit(matching_indices);
        data.elec.elecpos = data.elec.elecpos(matching_indices, :);
        data.elec.label = data.elec.label(matching_indices);
    end

    % Set the figure size to 20x15 mm (convert mm to inches, 1 inch = 25.4 mm)
    figure('Units', 'inches', 'Position', [0, 0, 0.95, 0.8]); % 20x15 mm

    % Adjust marker and line sizes
    cfg.markersize = 0.1; % Smaller marker size
    cfg.highlightfontsize = 0.3; % Font size for highlighted channels
    cfg.highlightsize = 0.1; % Smaller highlight marker size
    cfg.linewidth = 0.3; % Set the outline width
    cfg.highlightsymbol  = []; %'*'
    
    % Change marker to star    %% Plot positive cluster
    if strcmp(direction, 'pos')
        pos_cluster_pvals = [stat.posclusters(:).prob];
        pos_clust = find(pos_cluster_pvals < 0.05);
        pos = stat.posclusterslabelmat == clust_num; % or == 2, or 3, etc.

        [pos_a, pos_b] = find(pos);
        pos_chan = unique(pos_a);

        % If cluster does not extend in time it was probably averaged so take the latency from the old cfg
        if size(pos,2) == 1
            cfg.xlim = stat.cfg.latency;
        else
            cfg.xlim = [time_array(min(pos_b)) time_array(max(pos_b))];
        end

        if isempty(pos_clust)
            cfg.highlight = 'off';
        else
            cfg.highlight = 'on';
            cfg.highlightchannel = pos_chan;
        end

        ft_topoplotER(cfg, data);
        c = colorbar;

        %"round" colorbar ticks  
        min_val = ceil(c.Limits(1) * 10) / 10;
        max_val = floor(c.Limits(2) * 10) / 10;
        c.Ticks = unique([min_val, 0, max_val]);
        

        c.FontSize = 5;
        c.Position = [0.87, 0.38, 0.055, 0.3]; % Position: [x, y, width, height] in normalized units
        
        % Adjust colorbar outline thickness
        c.LineWidth = 0.2;

        if inc_title == 1
            title(['positive cluster number: ', num2str(clust_num), '; between (s)', num2str(cfg.xlim)], 'Interpreter', 'none', 'FontSize', 5);
        end
        yminmax = [min(min(data.avg(pos_chan, min(pos_b):max(pos_b)))), max(max(data.avg(pos_chan, min(pos_b):max(pos_b))))];

    elseif strcmp(direction, 'neg')
        %% Plot negative cluster
        neg_cluster_pvals = [stat.negclusters(:).prob];
        neg_clust = find(neg_cluster_pvals < 0.05);
        neg = stat.negclusterslabelmat == clust_num;

        [neg_a, neg_b] = find(neg);
        neg_chan = unique(neg_a);

        if size(neg,2) == 1
            cfg.xlim = stat.cfg.latency;
        else
            cfg.xlim = [time_array(min(neg_b)) time_array(max(neg_b))];
        end

        if isempty(neg_clust)
            cfg.highlight = 'off';
        else
            cfg.highlight = 'on';
            cfg.highlightchannel = neg_chan;
        end

        ft_topoplotER(cfg, data);
        c = colorbar;
        
        %round colorbar ticks  
        min_val = ceil(c.Limits(1) * 10) / 10;
        max_val = floor(c.Limits(2) * 10) / 10;
        c.Ticks = unique([min_val, 0, max_val]);
        
        c.FontSize = 5;
        c.Position = [0.87, 0.38, 0.055, 0.3]; % Position: [x, y, width, height] in normalized units
        
        % Adjust colorbar outline thickness
        c.LineWidth = 0.2;

        if inc_title == 1
            title(['negative cluster number: ', num2str(clust_num), '; between (s): ', num2str(cfg.xlim)], 'Interpreter', 'none', 'FontSize', 5);
        end

        yminmax = [min(min(data.avg(neg_chan, min(neg_b):max(neg_b)))), max(max(data.avg(neg_chan, min(neg_b):max(neg_b))))];
    end

    % Adjust the outline thickness
    h = findobj(gca, 'Type', 'contour');
    set(h, 'LineWidth', 0.1);


    h = findobj(gca, 'Type', 'line');
    set(h, 'LineWidth', 0.3);

    % Adjust the contour line thickness
    hContour = findobj(gca, 'Type', 'contour', '-and', 'LineStyle', '-');
    for i = 1:length(hContour)
        set(hContour(i), 'LineWidth', 0.1);
    end

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 5); % Set font size to 5 pt for all elements
end
