function sig_area = plot_erp_var(data_cond_1, data_cond_2, time_roi, chan_names, bsl_window, std_sem, mask, cluster_direction, legend_str, legend_location, y_ticks, savepath, savename)
    % PLOT_ERP_VAR Plots event-related potentials (ERPs) with variability and highlights significant areas.
    %
    % Description:
    %   This function plots ERPs for two conditions with shaded areas representing variability
    %   (standard deviation or standard error). It highlights significant time regions based on
    %   a statistical mask and allows for baseline correction. The plot can be customized with
    %   legends, axis labels, and saved as an SVG file.
    %
    % Inputs:
    %   data_cond_1       - Structure containing ERP data for condition 1 with fields 'avg', 'var', 'time', and optionally 'mask'.
    %   data_cond_2       - Structure containing ERP data for condition 2 with fields 'avg', 'var', 'time'.
    %   time_roi          - Vector specifying the time range of interest [start, end] in seconds.
    %   chan_names        - Cell array of channel names or vector of channel indices to average.
    %   bsl_window        - Vector specifying the baseline window [start, end] in seconds for correction.
    %   std_sem           - String specifying 'std' for standard deviation or 'sem' for standard error of the mean.
    %   mask              - Statistical mask structure or matrix for highlighting significant areas.
    %   cluster_direction - String specifying 'pos' for positive or 'neg' for negative clusters.
    %   legend_str        - Cell array of strings for legend entries.
    %   legend_location   - String specifying the location of the legend (e.g., 'northeast').
    %   y_ticks           - Vector specifying custom y-axis ticks.
    %   savepath          - String specifying the path to save the plot.
    %   savename          - String specifying the name of the saved plot file.
    %
    % Outputs:
    %   sig_area          - Matrix containing the start and end times of significant areas in milliseconds.
    %
    % Note:
    %   Ensure that the 'boundedline' and 'ft_plot_box' functions are available in your MATLAB environment.

    %% average across channels for each condition

    %get indices of channels to average
    if iscell(chan_names)
        ch_idx = find(contains(data_cond_1.label, chan_names));
    else
        ch_idx = chan_names;
    end

    avg_cond_1 = mean(data_cond_1.avg(ch_idx,:), 1);
    avg_cond_2 = mean(data_cond_2.avg(ch_idx,:), 1);

    %% get std or sem 
    if strcmp(std_sem, 'std')
        avg_sem_cond_1 = mean(sqrt(data_cond_1.var(ch_idx,:)), 1);
        avg_sem_cond_2 = mean(sqrt(data_cond_2.var(ch_idx,:)), 1);
    elseif strcmp(std_sem, 'sem')
        avg_sem_cond_1 = mean(sqrt(data_cond_1.var(ch_idx,:)), 1) / sqrt(size(data_cond_1.cfg.previous, 2));
        avg_sem_cond_2 = mean(sqrt(data_cond_2.var(ch_idx,:)), 1) / sqrt(size(data_cond_1.cfg.previous, 2));
    elseif isempty(std_sem)
        avg_sem_cond_1 = zeros(size(data_cond_1.avg));
        avg_sem_cond_2 = zeros(size(data_cond_2.avg));
    end

    %% timerange of interest
    if ~isempty(time_roi)
        time_roi_idx = data_cond_1.time > time_roi(1) & data_cond_1.time < time_roi(2);
    else
        time_roi_idx = ones(length(data_cond_1.time), 1)';
    end

    plot_time = data_cond_1.time(time_roi_idx);

    %% baseline correction
    if ~isempty(bsl_window)
        bsl_window_idx = data_cond_1.time > bsl_window(1) & data_cond_1.time < bsl_window(2);
        avg_cond_1 = avg_cond_1 - mean(avg_cond_1(bsl_window_idx));
        avg_cond_2 = avg_cond_2 - mean(avg_cond_2(bsl_window_idx));
    end

    %% plot the results
    % % rest vs.task
    % figureWidth_mm = 37; % width in mm
    % figureHeight_mm = 26; % height in mm

    % only rest or task
    figureWidth_mm = 80; % width in mm
    figureHeight_mm = 63.5; % height in mm

    figureWidth_pixels = figureWidth_mm * 3.77953; % Convert mm to pixels (1mm = 3.77953 pixels for 300 dpi)
    figureHeight_pixels = figureHeight_mm * 3.77953;

    % Create the figure with the specified dimensions
    fig = figure('Position', [100, 100, figureWidth_pixels, figureHeight_pixels]); % [left, bottom, width, height]

    h1 = boundedline(data_cond_1.time(time_roi_idx) * 1000, avg_cond_1(time_roi_idx), avg_sem_cond_1(time_roi_idx), 'cmap', [0.8500, 0.3250, 0.0980], 'alpha', 'LineWidth', 1);
    hold on

    h2 = boundedline(data_cond_2.time(time_roi_idx) * 1000, avg_cond_2(time_roi_idx), avg_sem_cond_2(time_roi_idx), 'cmap', [0, 0.4470, 0.7410], 'alpha', 'LineWidth', 1);
    ylabel('Amplitude (µV)', 'FontSize', 5);
    xlabel('Time (ms)', 'FontSize', 5);

    %% add significance shade
    if ~isempty(mask)
        % Check if we have the stats struct as mask input, if so we can choose either positive or negative clusters only
        if isstruct(mask)
            if isfield(mask, 'posclusters')
                if strcmp(cluster_direction, 'pos')
                    pos_clust = find([mask.posclusters.prob] < 0.05);
                    pos = ismember(mask.posclusterslabelmat, pos_clust);
                    if size(mask.label, 1) == 1 % changed from size(x,2) to size(x,1) because it was not plotting the shade - very weird.
                        highlight = mean(pos(1, time_roi_idx), 1) > 0; %if analysis was performed only on single channel, select that mask
                    else
                        highlight = mean(pos(ch_idx, time_roi_idx), 1) > 0;
                    end
                elseif strcmp(cluster_direction, 'neg')
                    if isfield(mask, 'negclusters')
                        neg_clust = find([mask.negclusters.prob] < 0.05);
                        neg = ismember(mask.negclusterslabelmat, neg_clust);
                        if size(mask.label, 1) == 1
                            highlight = mean(neg(1, time_roi_idx), 1) > 0;
                        else
                            highlight = mean(neg(ch_idx, time_roi_idx), 1) > 0;
                        end
                    end
                end
            end
        elseif isfield(data_cond_1, 'mask') % if no mask struct was added, take the mask matrix from the data struct
            highlight = mean(data_cond_1.mask(ch_idx, time_roi_idx), 1) > 0; %if any channel was significant in the ROI, make it 1
        end
    end

    % If we have a highlight, apply it
    if exist('highlight', 'var')
        % find the sample number where the highlight begins and ends
        begsample = find(diff([0 highlight 0]) == 1);
        endsample = find(diff([0 highlight 0]) == -1) - 1;

        sig_area = [plot_time(begsample), plot_time(endsample)];

        % Somehow ylim changes in the loop, so we keep the original
        ylim_orig = ylim;

        % Loop over sig clusters
        for i = 1:length(begsample)
            begx = plot_time(begsample(i)) * 1000;
            endx = plot_time(endsample(i)) * 1000;
            ft_plot_box([begx endx ylim], 'facecolor', [.9 .9 .9], 'edgecolor', 'none');
            set(gca, 'YLim', ylim_orig);
        end

        % Send patches to back
        set(gca, 'children', flipud(get(gca, 'children')));
    end

    set(gca, 'Layer', 'top');
    axis tight;

    % Put legend
    leg = legend([h1, h2], legend_str, 'Location', legend_location);
    leg = legend('boxoff');
    leg.ItemTokenSize = [10, 10]; % Adjust the legend box size

    % Change font size
    ax = gca; % Get handle to the current axes
    ax.XAxis.FontSize = 5; % Set the font size of the x-axis ticks
    ax.YAxis.FontSize = 5; % Set the font size of the y-axis ticks

    % Add y ticks
    if ~isempty(y_ticks)
        yticks(y_ticks);
    end

    %% Save output
    if ischar(savepath)
        % Use fullfile to properly construct the path
        savefullpath = fullfile(savepath, savename);
        % Change extension from .png to .svg if needed
        [~, name, ~] = fileparts(savename);
        savefullpath = fullfile(savepath, [name, '.svg']);
        print('-dsvg', savefullpath);
    end
end
