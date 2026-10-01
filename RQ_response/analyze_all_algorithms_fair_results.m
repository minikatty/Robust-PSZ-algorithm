function results = analyze_all_algorithms_fair_results(run_dir)
%ANALYZE_ALL_ALGORITHMS_FAIR_RESULTS Summarize the frozen final comparison.

    response_dir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(run_dir)
        run_dir = local_latest_complete_run(fullfile(response_dir, ...
            'all_algorithms_fair_results'));
    end
    run_dir = char(string(run_dir));
    if ~exist(fullfile(run_dir, 'RUN_COMPLETE.txt'), 'file')
        error('FairAnalysis:IncompleteRun', ...
            'RUN_COMPLETE.txt is missing from %s.', run_dir);
    end

    result_files = dir(fullfile(run_dir, 'performance_matrices_*.mat'));
    if isempty(result_files)
        error('FairAnalysis:MissingResult', ...
            'No performance_matrices_*.mat exists in %s.', run_dir);
    end
    [~, newest] = max([result_files.datenum]);
    result_file = fullfile(result_files(newest).folder, ...
        result_files(newest).name);
    loaded = load(result_file, 'val_results', 'para', 'freq_params', ...
        'algorithms_to_test', 'metrics_to_evaluate', ...
        'uncertainty_config', 'algorithm_wall_seconds');

    temperatures = loaded.para.temperature_vector_celsius(:);
    frequencies = loaded.freq_params.target_freqs(:).';
    requested_frequencies = [500, 1000, 4000];
    frequency_indices = local_frequency_indices( ...
        frequencies, requested_frequencies);
    nominal_index = find(abs(temperatures - 22.5) < 1e-10, 1);
    if isempty(nominal_index)
        error('FairAnalysis:MissingNominalTemperature', ...
            'The result does not contain 22.5 deg C.');
    end

    analysis_dir = fullfile(run_dir, 'analysis');
    if ~exist(analysis_dir, 'dir'), mkdir(analysis_dir); end
    summary = local_summary(loaded.val_results, ...
        loaded.algorithms_to_test, requested_frequencies, ...
        frequency_indices, nominal_index);
    writetable(summary, fullfile(analysis_dir, ...
        'selected_frequency_summary.csv'));

    advantage = local_racc_pm_advantage(summary);
    writetable(advantage, fullfile(analysis_dir, ...
        'racc_pm_pairwise_differences.csv'));

    band_mask = frequencies >= 100 & frequencies <= 4000;
    band_summary = local_band_summary(loaded.val_results, ...
        loaded.algorithms_to_test, band_mask, nominal_index);
    writetable(band_summary, fullfile(analysis_dir, ...
        'band_100_4000_summary.csv'));

    near_nominal_summary = local_near_nominal_summary(loaded.val_results, ...
        loaded.algorithms_to_test, temperatures, band_mask, ...
        [22.2, 22.5, 22.8]);
    writetable(near_nominal_summary, fullfile(analysis_dir, ...
        'near_nominal_band_summary.csv'));

    figure_paths = local_plot_all_metrics(loaded.val_results, ...
        loaded.algorithms_to_test, temperatures, requested_frequencies, ...
        frequency_indices, analysis_dir);
    band_figure_paths = local_plot_band_metrics(loaded.val_results, ...
        loaded.algorithms_to_test, temperatures, band_mask, analysis_dir);

    timing = local_timing_table(loaded.algorithm_wall_seconds, ...
        loaded.algorithms_to_test);
    writetable(timing, fullfile(analysis_dir, ...
        'algorithm_design_wall_time.csv'));

    results = struct('run_dir', run_dir, 'result_file', result_file, ...
        'analysis_dir', analysis_dir, 'summary', summary, ...
        'racc_pm_advantage', advantage, 'band_summary', band_summary, ...
        'near_nominal_summary', near_nominal_summary, 'timing', timing, ...
        'figure_paths', {figure_paths}, ...
        'band_figure_paths', {band_figure_paths}, ...
        'uncertainty_config', loaded.uncertainty_config);
    save(fullfile(analysis_dir, 'analysis_results.mat'), 'results', '-v7.3');
    disp(summary);
    disp(advantage);
    disp(band_summary);
    disp(near_nominal_summary);
    fprintf('Completed fair-result analysis: %s\n', analysis_dir);
end

function summary = local_band_summary(performance, algorithms, ...
        band_mask, nominal_index)
    Algorithm = string(algorithms(:));
    rows = numel(algorithms);
    NominalMeanACdB = nan(rows, 1);
    SweepMeanACdB = nan(rows, 1);
    WorstTemperatureMeanACdB = nan(rows, 1);
    TemperatureMeanACRangedB = nan(rows, 1);
    NominalMeanNSREdB = nan(rows, 1);
    SweepMeanNSREdB = nan(rows, 1);
    WorstTemperatureMeanNSREdB = nan(rows, 1);
    NominalMeanAEdB = nan(rows, 1);
    SweepMeanAEdB = nan(rows, 1);
    WorstTemperatureMeanAEdB = nan(rows, 1);
    NominalMeanPlanarityPercent = nan(rows, 1);
    SweepMeanPlanarityPercent = nan(rows, 1);
    WorstTemperatureMeanPlanarityPercent = nan(rows, 1);

    for a = 1:rows
        algorithm = algorithms{a};
        ac_by_temperature = mean(performance.(algorithm).AC(:, band_mask), ...
            2, 'omitnan');
        nsre_by_temperature = mean( ...
            performance.(algorithm).NSRE(:, band_mask), 2, 'omitnan');
        ae_by_temperature = mean(performance.(algorithm).AE(:, band_mask), ...
            2, 'omitnan');
        planarity_by_temperature = mean( ...
            performance.(algorithm).Planarity(:, band_mask), 2, 'omitnan');
        NominalMeanACdB(a) = ac_by_temperature(nominal_index);
        SweepMeanACdB(a) = mean(ac_by_temperature, 'omitnan');
        WorstTemperatureMeanACdB(a) = min(ac_by_temperature, [], 'omitnan');
        TemperatureMeanACRangedB(a) = ...
            max(ac_by_temperature, [], 'omitnan') - ...
            min(ac_by_temperature, [], 'omitnan');
        NominalMeanNSREdB(a) = nsre_by_temperature(nominal_index);
        SweepMeanNSREdB(a) = mean(nsre_by_temperature, 'omitnan');
        WorstTemperatureMeanNSREdB(a) = ...
            max(nsre_by_temperature, [], 'omitnan');
        NominalMeanAEdB(a) = ae_by_temperature(nominal_index);
        SweepMeanAEdB(a) = mean(ae_by_temperature, 'omitnan');
        WorstTemperatureMeanAEdB(a) = ...
            max(ae_by_temperature, [], 'omitnan');
        NominalMeanPlanarityPercent(a) = ...
            planarity_by_temperature(nominal_index);
        SweepMeanPlanarityPercent(a) = ...
            mean(planarity_by_temperature, 'omitnan');
        WorstTemperatureMeanPlanarityPercent(a) = ...
            min(planarity_by_temperature, [], 'omitnan');
    end
    summary = table(Algorithm, NominalMeanACdB, SweepMeanACdB, ...
        WorstTemperatureMeanACdB, TemperatureMeanACRangedB, ...
        NominalMeanNSREdB, SweepMeanNSREdB, ...
        WorstTemperatureMeanNSREdB, NominalMeanAEdB, SweepMeanAEdB, ...
        WorstTemperatureMeanAEdB, NominalMeanPlanarityPercent, ...
        SweepMeanPlanarityPercent, WorstTemperatureMeanPlanarityPercent);
end

function summary = local_near_nominal_summary(performance, algorithms, ...
        temperatures, band_mask, requested_temperatures)
    rows = numel(algorithms) * numel(requested_temperatures);
    Algorithm = strings(rows, 1);
    TemperatureC = nan(rows, 1);
    MeanACdB = nan(rows, 1);
    MeanNSREdB = nan(rows, 1);
    MeanAEdB = nan(rows, 1);
    MeanPlanarityPercent = nan(rows, 1);
    row = 0;
    for a = 1:numel(algorithms)
        for t = 1:numel(requested_temperatures)
            [distance, temperature_index] = min(abs(temperatures - ...
                requested_temperatures(t)));
            if distance > 1e-10
                error('FairAnalysis:TemperatureMissing', ...
                    'Requested temperature %.1f C is unavailable.', ...
                    requested_temperatures(t));
            end
            row = row + 1;
            algorithm = algorithms{a};
            Algorithm(row) = string(algorithm);
            TemperatureC(row) = requested_temperatures(t);
            MeanACdB(row) = mean(performance.(algorithm).AC( ...
                temperature_index, band_mask), 'omitnan');
            MeanNSREdB(row) = mean(performance.(algorithm).NSRE( ...
                temperature_index, band_mask), 'omitnan');
            MeanAEdB(row) = mean(performance.(algorithm).AE( ...
                temperature_index, band_mask), 'omitnan');
            MeanPlanarityPercent(row) = mean( ...
                performance.(algorithm).Planarity( ...
                temperature_index, band_mask), 'omitnan');
        end
    end
    summary = table(Algorithm, TemperatureC, MeanACdB, MeanNSREdB, ...
        MeanAEdB, MeanPlanarityPercent);
end

function run_dir = local_latest_complete_run(root)
    directories = dir(fullfile(root, 'run_*'));
    directories = directories([directories.isdir]);
    complete = false(size(directories));
    for k = 1:numel(directories)
        complete(k) = exist(fullfile(directories(k).folder, ...
            directories(k).name, 'RUN_COMPLETE.txt'), 'file') == 2;
    end
    directories = directories(complete);
    if isempty(directories)
        error('FairAnalysis:NoCompleteRun', ...
            'No completed fair run exists under %s.', root);
    end
    [~, newest] = max([directories.datenum]);
    run_dir = fullfile(directories(newest).folder, ...
        directories(newest).name);
end

function indices = local_frequency_indices(all_frequencies, requested)
    indices = zeros(size(requested));
    for k = 1:numel(requested)
        [distance, indices(k)] = min(abs(all_frequencies - requested(k)));
        if distance > 1e-8
            error('FairAnalysis:FrequencyMissing', ...
                'Requested frequency %.1f Hz is unavailable.', requested(k));
        end
    end
end

function summary = local_summary(performance, algorithms, ...
        requested_frequencies, frequency_indices, nominal_index)
    rows = numel(algorithms) * numel(requested_frequencies);
    Algorithm = strings(rows, 1);
    FrequencyHz = nan(rows, 1);
    NominalACdB = nan(rows, 1);
    MeanACdB = nan(rows, 1);
    MinimumACdB = nan(rows, 1);
    ACRangedB = nan(rows, 1);
    NominalNSREdB = nan(rows, 1);
    MeanNSREdB = nan(rows, 1);
    MaximumNSREdB = nan(rows, 1);
    NominalAEdB = nan(rows, 1);
    MeanAEdB = nan(rows, 1);
    MaximumAEdB = nan(rows, 1);
    NominalPlanarityPercent = nan(rows, 1);
    MeanPlanarityPercent = nan(rows, 1);
    MinimumPlanarityPercent = nan(rows, 1);

    row = 0;
    for a = 1:numel(algorithms)
        algorithm = algorithms{a};
        for f = 1:numel(requested_frequencies)
            row = row + 1;
            index = frequency_indices(f);
            ac = performance.(algorithm).AC(:, index);
            nsre = performance.(algorithm).NSRE(:, index);
            ae = performance.(algorithm).AE(:, index);
            planarity = performance.(algorithm).Planarity(:, index);
            Algorithm(row) = string(algorithm);
            FrequencyHz(row) = requested_frequencies(f);
            NominalACdB(row) = ac(nominal_index);
            MeanACdB(row) = mean(ac, 'omitnan');
            MinimumACdB(row) = min(ac, [], 'omitnan');
            ACRangedB(row) = max(ac, [], 'omitnan') - ...
                min(ac, [], 'omitnan');
            NominalNSREdB(row) = nsre(nominal_index);
            MeanNSREdB(row) = mean(nsre, 'omitnan');
            MaximumNSREdB(row) = max(nsre, [], 'omitnan');
            NominalAEdB(row) = ae(nominal_index);
            MeanAEdB(row) = mean(ae, 'omitnan');
            MaximumAEdB(row) = max(ae, [], 'omitnan');
            NominalPlanarityPercent(row) = planarity(nominal_index);
            MeanPlanarityPercent(row) = mean(planarity, 'omitnan');
            MinimumPlanarityPercent(row) = min(planarity, [], 'omitnan');
        end
    end
    summary = table(Algorithm, FrequencyHz, NominalACdB, MeanACdB, ...
        MinimumACdB, ACRangedB, NominalNSREdB, MeanNSREdB, ...
        MaximumNSREdB, NominalAEdB, MeanAEdB, MaximumAEdB, ...
        NominalPlanarityPercent, MeanPlanarityPercent, ...
        MinimumPlanarityPercent);
end

function advantage = local_racc_pm_advantage(summary)
    proposed_name = "RACC_PM_Subpro";
    baselines = unique(summary.Algorithm(summary.Algorithm ~= proposed_name), ...
        'stable');
    frequencies = unique(summary.FrequencyHz, 'stable');
    rows = numel(baselines) * numel(frequencies);
    Baseline = strings(rows, 1);
    FrequencyHz = nan(rows, 1);
    DeltaMeanACdB = nan(rows, 1);
    DeltaMinimumACdB = nan(rows, 1);
    DeltaACRangedB = nan(rows, 1);
    DeltaMeanNSREdB = nan(rows, 1);
    DeltaMaximumNSREdB = nan(rows, 1);
    DeltaMeanAEdB = nan(rows, 1);
    DeltaMeanPlanarityPercent = nan(rows, 1);

    row = 0;
    for b = 1:numel(baselines)
        for f = 1:numel(frequencies)
            row = row + 1;
            proposed = summary(summary.Algorithm == proposed_name & ...
                summary.FrequencyHz == frequencies(f), :);
            baseline = summary(summary.Algorithm == baselines(b) & ...
                summary.FrequencyHz == frequencies(f), :);
            Baseline(row) = baselines(b);
            FrequencyHz(row) = frequencies(f);
            DeltaMeanACdB(row) = proposed.MeanACdB - baseline.MeanACdB;
            DeltaMinimumACdB(row) = ...
                proposed.MinimumACdB - baseline.MinimumACdB;
            DeltaACRangedB(row) = proposed.ACRangedB - baseline.ACRangedB;
            DeltaMeanNSREdB(row) = ...
                proposed.MeanNSREdB - baseline.MeanNSREdB;
            DeltaMaximumNSREdB(row) = ...
                proposed.MaximumNSREdB - baseline.MaximumNSREdB;
            DeltaMeanAEdB(row) = proposed.MeanAEdB - baseline.MeanAEdB;
            DeltaMeanPlanarityPercent(row) = ...
                proposed.MeanPlanarityPercent - ...
                baseline.MeanPlanarityPercent;
        end
    end
    advantage = table(Baseline, FrequencyHz, DeltaMeanACdB, ...
        DeltaMinimumACdB, DeltaACRangedB, DeltaMeanNSREdB, ...
        DeltaMaximumNSREdB, DeltaMeanAEdB, ...
        DeltaMeanPlanarityPercent);
end

function paths = local_plot_all_metrics(performance, algorithms, ...
        temperatures, requested_frequencies, frequency_indices, output_dir)
    specifications = { ...
        'AC', 'Acoustic contrast (dB)', 'ac_vs_temperature'; ...
        'NSRE', 'NSRE (dB)', 'nsre_vs_temperature'; ...
        'AE', 'Array effort (dB)', 'ae_vs_temperature'; ...
        'Planarity', 'Planarity (%)', 'planarity_vs_temperature'};
    colors = lines(numel(algorithms));
    line_styles = {'-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--'};
    markers = {'o', 's', '^', 'v', 'd', 'p', '<', '>', 'h', 'x'};
    paths = cell(size(specifications, 1), 3);

    for metric_index = 1:size(specifications, 1)
        metric = specifications{metric_index, 1};
        fig = figure('Color', 'w', 'Position', [100, 100, 1450, 430]);
        layout = tiledlayout(fig, 1, numel(requested_frequencies), ...
            'TileSpacing', 'compact', 'Padding', 'compact');
        legend_handles = gobjects(0);
        legend_labels = strings(0);
        first_ax = gobjects(1);
        for f = 1:numel(requested_frequencies)
            ax = nexttile(layout);
            if f == 1, first_ax = ax; end
            hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
            disableDefaultInteractivity(ax);
            if ~isempty(ax.Toolbar), ax.Toolbar.Visible = 'off'; end
            for a = 1:numel(algorithms)
                values = performance.(algorithms{a}).(metric)(:, ...
                    frequency_indices(f));
                if all(isnan(values)), continue; end
                handle = plot(ax, temperatures, values, ...
                    'LineWidth', 1.25, 'Color', colors(a, :), ...
                    'LineStyle', line_styles{a}, ...
                    'Marker', markers{a}, 'MarkerIndices', 1:5:numel(values), ...
                    'MarkerSize', 4);
                if f == 1
                    legend_handles(end + 1) = handle; %#ok<AGROW>
                    legend_labels(end + 1) = ...
                        local_display_name(algorithms{a}); %#ok<AGROW>
                end
            end
            xline(ax, 22.5, ':k', 'HandleVisibility', 'off');
            xlabel(ax, 'Temperature (^{\circ}C)');
            ylabel(ax, specifications{metric_index, 2});
            title(ax, sprintf('%d Hz', requested_frequencies(f)));
            set(ax, 'FontName', 'Times New Roman', 'FontSize', 10);
        end
        legend_handle = legend(first_ax, legend_handles, ...
            cellstr(legend_labels), ...
            'Orientation', 'horizontal', 'NumColumns', 5, ...
            'Location', 'southoutside');
        if isprop(legend_handle, 'Layout')
            legend_handle.Layout.Tile = 'south';
        end
        base = fullfile(output_dir, specifications{metric_index, 3});
        paths{metric_index, 1} = [base '.png'];
        paths{metric_index, 2} = [base '.pdf'];
        paths{metric_index, 3} = [base '.fig'];
        exportgraphics(fig, paths{metric_index, 1}, 'Resolution', 300);
        exportgraphics(fig, paths{metric_index, 2}, 'ContentType', 'vector');
        savefig(fig, paths{metric_index, 3});
        close(fig);
    end
end

function paths = local_plot_band_metrics(performance, algorithms, ...
        temperatures, band_mask, output_dir)
    specifications = { ...
        'AC', 'Mean acoustic contrast (dB)'; ...
        'NSRE', 'Mean NSRE (dB)'; ...
        'AE', 'Mean array effort (dB)'; ...
        'Planarity', 'Mean planarity (%)'};
    colors = lines(numel(algorithms));
    line_styles = {'-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--'};
    markers = {'o', 's', '^', 'v', 'd', 'p', '<', '>', 'h', 'x'};
    fig = figure('Color', 'w', 'Position', [100, 100, 1180, 760]);
    layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
        'Padding', 'compact');
    legend_handles = gobjects(0);
    legend_labels = strings(0);
    first_ax = gobjects(1);
    for metric_index = 1:size(specifications, 1)
        metric = specifications{metric_index, 1};
        ax = nexttile(layout);
        if metric_index == 1, first_ax = ax; end
        hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
        disableDefaultInteractivity(ax);
        if ~isempty(ax.Toolbar), ax.Toolbar.Visible = 'off'; end
        for a = 1:numel(algorithms)
            values = mean(performance.(algorithms{a}).(metric)(:, band_mask), ...
                2, 'omitnan');
            if all(isnan(values)), continue; end
            handle = plot(ax, temperatures, values, ...
                'LineWidth', 1.25, 'Color', colors(a, :), ...
                'LineStyle', line_styles{a}, 'Marker', markers{a}, ...
                'MarkerIndices', 1:5:numel(values), 'MarkerSize', 4);
            if metric_index == 1
                legend_handles(end + 1) = handle; %#ok<AGROW>
                legend_labels(end + 1) = ...
                    local_display_name(algorithms{a}); %#ok<AGROW>
            end
        end
        xline(ax, 22.5, ':k', 'HandleVisibility', 'off');
        xlabel(ax, 'Temperature (^{\circ}C)');
        ylabel(ax, specifications{metric_index, 2});
        title(ax, sprintf('%s, 100--4000 Hz', ...
            local_metric_title(metric)));
        set(ax, 'FontName', 'Times New Roman', 'FontSize', 10);
    end
    legend_handle = legend(first_ax, legend_handles, ...
        cellstr(legend_labels), 'Orientation', 'horizontal', ...
        'NumColumns', 5, 'Location', 'southoutside');
    if isprop(legend_handle, 'Layout')
        legend_handle.Layout.Tile = 'south';
    end
    base = fullfile(output_dir, 'band_100_4000_metrics_vs_temperature');
    paths = {[base '.png'], [base '.pdf'], [base '.fig']};
    exportgraphics(fig, paths{1}, 'Resolution', 300);
    exportgraphics(fig, paths{2}, 'ContentType', 'vector');
    savefig(fig, paths{3});
    close(fig);
end

function title_text = local_metric_title(metric)
    switch metric
        case 'AC', title_text = 'Acoustic contrast';
        case 'NSRE', title_text = 'Reproduction error';
        case 'AE', title_text = 'Array effort';
        case 'Planarity', title_text = 'Planarity';
        otherwise, title_text = metric;
    end
end

function name = local_display_name(identifier)
    switch identifier
        case 'ACC_Reg', name = "ACC-Reg";
        case 'ACC_PM', name = "ACC-PM";
        case 'NoCT_WCRACC', name = "NoCT-WCRACC";
        case 'Full_WCRACC', name = "Full-WCRACC";
        case 'POTDC_RACC', name = "POTDC-RACC";
        case 'RACC_PM_Subpro', name = "RACC-PM";
        otherwise, name = string(strrep(identifier, '_', '-'));
    end
end

function timing = local_timing_table(wall_times, algorithms)
    Algorithm = string(algorithms(:));
    WallTimeSeconds = nan(numel(algorithms), 1);
    for k = 1:numel(algorithms)
        if isfield(wall_times, algorithms{k})
            WallTimeSeconds(k) = wall_times.(algorithms{k});
        end
    end
    WallTimeMinutes = WallTimeSeconds / 60;
    timing = table(Algorithm, WallTimeSeconds, WallTimeMinutes);
end
