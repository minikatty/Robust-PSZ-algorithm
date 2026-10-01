function output = postprocess_racc_pm_nu_temperature_severity(run_dir)
%POSTPROCESS_RACC_PM_NU_TEMPERATURE_SEVERITY Reviewer-facing nu sweep plot.
%   Uses only the RACC-PM filters from RUN_DIR. The frozen filters are
%   evaluated on independent BZ/DZ evaluation points and separated into
%   mild (|T-T0| <= 0.3 C) and severe (|T-T0| > 0.3 C) temperature sets.
%   Native NSRE is computed against the evaluation-grid plane-wave target;
%   no post-hoc complex gain or phase alignment is applied.

    if nargin < 1 || strlength(string(run_dir)) == 0
        error('RACCPMNuSeverity:MissingRunDirectory', ...
            'A completed uncertainty-sweep run directory is required.');
    end
    run_dir = char(string(run_dir));
    result_file = fullfile(run_dir, 'experiment_results.mat');
    if ~isfile(result_file)
        error('RACCPMNuSeverity:MissingResults', ...
            'Missing experiment result: %s', result_file);
    end

    response_dir = fileparts(mfilename('fullpath'));
    project_root = fileparts(response_dir);
    addpath(fullfile(project_root, 'src', 'utils'));
    loaded = load(result_file, 'results');
    source = loaded.results;

    expected_frequencies = [500, 1000, 1500, 2000, 4000];
    if ~isequal(source.frequencies_hz(:).', expected_frequencies)
        error('RACCPMNuSeverity:UnexpectedFrequencies', ...
            'Expected frequencies %s Hz, but the run contains %s Hz.', ...
            mat2str(expected_frequencies), ...
            mat2str(source.frequencies_hz(:).'));
    end
    settings = source.settings;
    is_constant = strcmp({settings.type}, 'common-relative');
    settings = settings(is_constant);
    designs = source.designs(is_constant);
    nu_values = [settings.nominal_level];
    [nu_values, order] = sort(nu_values);
    settings = settings(order);
    designs = designs(order);

    frequencies_hz = source.frequencies_hz(:).';
    frequency_indices = source.frequency_indices(:).';
    temperatures_c = source.config.test_temperatures_c(:).';
    nominal_temperature_c = source.config.nominal_temperature_c;
    mild_limit_c = 0.3;

    data_dir = fullfile(project_root, 'data', 'SimulateRIR', ...
        'temperature');
    metadata = load(fullfile(data_dir, 'para.mat'), 'para');
    geometry = load(fullfile(project_root, 'data', 'arrayGeometry', ...
        'array_layout.mat'), 'roomArray');
    desired_saved = load(fullfile(project_root, 'data', ...
        'ATF_desired_plane.mat'), 'ATF_desired_plane');

    position = struct('src', geometry.roomArray.s, ...
        'pos_center', [2.6, 2.0, 1.6], ...
        'virtual_src_idx', source.config.virtual_source_idx, ...
        'mics_pos', geometry.roomArray.BZ_ctrl);
    desired_control_all = local_plane_wave_generator(position, ...
        metadata.para.freq_params);
    desired_control = desired_control_all(:, frequency_indices);
    saved_control = desired_saved.ATF_desired_plane(:, frequency_indices);
    desired_error = max(abs(desired_control(:) - saved_control(:)));
    if desired_error > 1e-10
        error('RACCPMNuSeverity:DesiredFieldMismatch', ...
            'Generated and saved control targets differ by %.3g.', ...
            desired_error);
    end
    position.mics_pos = geometry.roomArray.BZ_eval;
    desired_evaluation_all = local_plane_wave_generator(position, ...
        metadata.para.freq_params);
    desired_evaluation = desired_evaluation_all(:, frequency_indices);

    num_settings = numel(settings);
    num_frequencies = numel(frequencies_hz);
    num_temperatures = numel(temperatures_c);
    num_rows = num_settings * num_frequencies * num_temperatures;
    Nu = nan(num_rows, 1);
    FrequencyHz = nan(num_rows, 1);
    TemperatureC = nan(num_rows, 1);
    Severity = strings(num_rows, 1);
    EvaluationACdB = nan(num_rows, 1);
    EvaluationNSREdB = nan(num_rows, 1);
    row = 0;

    for temperature_index = 1:num_temperatures
        temperature_c = temperatures_c(temperature_index);
        test_file = get_data_filename(data_dir, 'temperature', ...
            temperature_c);
        test = load(test_file, 'ATF_BZ', 'ATF_DZ');
        HB_evaluation = test.ATF_BZ.eval(:, :, frequency_indices);
        HD_evaluation = test.ATF_DZ.eval(:, :, frequency_indices);
        if abs(temperature_c - nominal_temperature_c) <= ...
                mild_limit_c + 1e-9
            severity = "Mild";
        else
            severity = "Severe";
        end

        for setting_index = 1:num_settings
            filters = designs{setting_index}.filters.RACC_PM;
            for frequency_index = 1:num_frequencies
                row = row + 1;
                w = filters(:, frequency_index);
                if any(~isfinite(w)) || norm(w) <= realmin
                    error('RACCPMNuSeverity:InvalidFilter', ...
                        ['Invalid RACC-PM filter at nu=%.4g and ' ...
                        'frequency %.0f Hz.'], nu_values(setting_index), ...
                        frequencies_hz(frequency_index));
                end
                HBe = HB_evaluation(:, :, frequency_index);
                HDe = HD_evaluation(:, :, frequency_index);
                desired = desired_evaluation(:, frequency_index);
                Nu(row) = nu_values(setting_index);
                FrequencyHz(row) = frequencies_hz(frequency_index);
                TemperatureC(row) = temperature_c;
                Severity(row) = severity;
                EvaluationACdB(row) = local_ac_db(w, HBe, HDe);
                EvaluationNSREdB(row) = local_nsre_db(w, HBe, desired);
            end
        end
        fprintf('Severity postprocess %.1f C (%d/%d).\n', ...
            temperature_c, temperature_index, num_temperatures);
    end

    metrics = table(Nu, FrequencyHz, TemperatureC, Severity, ...
        EvaluationACdB, EvaluationNSREdB);
    metrics_file = fullfile(run_dir, ...
        'racc_pm_nu_temperature_severity_metrics.csv');
    writetable(metrics, metrics_file);
    summary = local_summary(metrics, nu_values, frequencies_hz);
    summary_file = fullfile(run_dir, ...
        'racc_pm_nu_temperature_severity_summary.csv');
    writetable(summary, summary_file);
    figure_paths = local_plot(summary, nu_values, frequencies_hz, run_dir);

    output = struct('source_result_file', result_file, ...
        'run_dir', run_dir, 'frequencies_hz', frequencies_hz, ...
        'nu_values', nu_values, 'nominal_temperature_c', ...
        nominal_temperature_c, 'mild_limit_c', mild_limit_c, ...
        'desired_control_validation_error', desired_error, ...
        'metrics', metrics, 'summary', summary, ...
        'metrics_file', metrics_file, 'summary_file', summary_file, ...
        'figure_paths', figure_paths);
    save(fullfile(run_dir, ...
        'racc_pm_nu_temperature_severity_analysis.mat'), ...
        'output', '-v7.3');
    fprintf('Completed temperature-severity analysis: %s\n', run_dir);
end

function desired = local_plane_wave_generator(position, freq_params)
    source_position = position.src(position.virtual_src_idx, :);
    direction = position.pos_center - source_position;
    direction = direction / norm(direction);
    wave_numbers = 2 * pi * freq_params.target_freqs(:) / 343;
    phase_delay = position.mics_pos * (wave_numbers .* direction).';
    desired = exp(-1i * phase_delay);
end

function value = local_ac_db(w, HB, HD)
    bright_energy = norm(HB * w)^2;
    dark_energy = norm(HD * w)^2;
    contrast = size(HD, 1) * bright_energy / max( ...
        size(HB, 1) * dark_energy, realmin);
    value = 10 * log10(max(real(contrast), realmin));
end

function value = local_nsre_db(w, HB, desired)
    error_energy = norm(HB * w - desired)^2;
    target_energy = norm(desired)^2;
    value = 10 * log10(max(error_energy / max( ...
        target_energy, realmin), realmin));
end

function summary = local_summary(metrics, nu_values, frequencies_hz)
    severity_names = ["Mild", "Severe"];
    rows = numel(nu_values) * numel(frequencies_hz) * ...
        numel(severity_names);
    Nu = nan(rows, 1);
    FrequencyHz = nan(rows, 1);
    Severity = strings(rows, 1);
    NumTemperatures = nan(rows, 1);
    MeanEvaluationACdB = nan(rows, 1);
    MinimumEvaluationACdB = nan(rows, 1);
    MeanEvaluationNSREdB = nan(rows, 1);
    MaximumEvaluationNSREdB = nan(rows, 1);
    row = 0;
    for nu_index = 1:numel(nu_values)
        for frequency_index = 1:numel(frequencies_hz)
            for severity_index = 1:numel(severity_names)
                row = row + 1;
                mask = abs(metrics.Nu - nu_values(nu_index)) < 1e-12 & ...
                    metrics.FrequencyHz == frequencies_hz(frequency_index) & ...
                    metrics.Severity == severity_names(severity_index);
                values = metrics(mask, :);
                Nu(row) = nu_values(nu_index);
                FrequencyHz(row) = frequencies_hz(frequency_index);
                Severity(row) = severity_names(severity_index);
                NumTemperatures(row) = height(values);
                MeanEvaluationACdB(row) = ...
                    mean(values.EvaluationACdB, 'omitnan');
                MinimumEvaluationACdB(row) = ...
                    min(values.EvaluationACdB, [], 'omitnan');
                MeanEvaluationNSREdB(row) = ...
                    mean(values.EvaluationNSREdB, 'omitnan');
                MaximumEvaluationNSREdB(row) = ...
                    max(values.EvaluationNSREdB, [], 'omitnan');
            end
        end
    end
    summary = table(Nu, FrequencyHz, Severity, NumTemperatures, ...
        MeanEvaluationACdB, MinimumEvaluationACdB, ...
        MeanEvaluationNSREdB, MaximumEvaluationNSREdB);
end

function paths = local_plot(summary, nu_values, frequencies_hz, output_dir)
    palette = [0.0000, 0.4470, 0.7410; ...
               0.8500, 0.3250, 0.0980; ...
               0.9290, 0.6940, 0.1250; ...
               0.4940, 0.1840, 0.5560; ...
               0.4660, 0.6740, 0.1880];
    markers = {'o', 's', '^', 'd', 'v'};
    severity_names = ["Mild", "Severe"];
    severity_styles = {'-', '--'};
    metric_fields = {'MeanEvaluationNSREdB', 'MeanEvaluationACdB'};
    y_labels = {'NSRE: [dB]', 'AC: [dB]'};
    titles = {'(a) Normalized reproduction error', ...
        '(b) Acoustic contrast'};

    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [80, 80, 1320, 500], ...
        'ToolBar', 'figure', 'MenuBar', 'figure');
    layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', ...
        'Padding', 'compact');
    axes_handles = gobjects(1, 2);
    frequency_handles = gobjects(1, numel(frequencies_hz));
    severity_handles = gobjects(1, numel(severity_names));

    for metric_index = 1:2
        ax = nexttile(layout);
        axes_handles(metric_index) = ax;
        hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
        disableDefaultInteractivity(ax);
        if ~isempty(ax.Toolbar), ax.Toolbar.Visible = 'off'; end
        for frequency_index = 1:numel(frequencies_hz)
            for severity_index = 1:numel(severity_names)
                mask = summary.FrequencyHz == ...
                    frequencies_hz(frequency_index) & ...
                    summary.Severity == severity_names(severity_index);
                subset = sortrows(summary(mask, :), 'Nu');
                semilogx(ax, subset.Nu, ...
                    subset.(metric_fields{metric_index}), ...
                    'LineStyle', severity_styles{severity_index}, ...
                    'Color', palette(frequency_index, :), ...
                    'LineWidth', 1.55, ...
                    'Marker', markers{frequency_index}, ...
                    'MarkerSize', 5.0, 'MarkerFaceColor', 'w', ...
                    'MarkerEdgeColor', palette(frequency_index, :), ...
                    'HandleVisibility', 'off');
            end
        end
        xlabel(ax, 'Normalized ATF radius \nu', ...
            'FontName', 'Times New Roman', 'FontSize', 12.5, ...
            'FontWeight', 'bold');
        ylabel(ax, y_labels{metric_index}, ...
            'FontName', 'Times New Roman', 'FontSize', 12.5, ...
            'FontWeight', 'bold');
        title(ax, titles{metric_index}, ...
            'FontName', 'Times New Roman', 'FontSize', 13.5, ...
            'FontWeight', 'bold');
        set(ax, 'FontName', 'Times New Roman', 'FontSize', 10.5, ...
            'LineWidth', 1, 'TickDir', 'in', 'Layer', 'top', ...
            'XScale', 'log', 'XLim', [min(nu_values), max(nu_values)], ...
            'XTick', nu_values, 'XTickLabel', compose('%.4g', nu_values), ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'XMinorGrid', 'on', 'YMinorGrid', 'on', ...
            'GridColor', [0.30, 0.30, 0.30], ...
            'MinorGridColor', [0.38, 0.38, 0.38], ...
            'GridAlpha', 0.31, 'MinorGridAlpha', 0.18, ...
            'GridLineStyle', '--', 'MinorGridLineStyle', ':', ...
            'XColor', [0.12, 0.12, 0.12], ...
            'YColor', [0.12, 0.12, 0.12]);
        drawnow;
        y_limits = ylim(ax);
        stable_patch = patch(ax, [0.005, 0.02, 0.02, 0.005], ...
            [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
            [0.78, 0.93, 0.82], 'EdgeColor', [0.25, 0.62, 0.35], ...
            'LineStyle', ':', 'LineWidth', 0.8, ...
            'FaceAlpha', 0.52, 'HandleVisibility', 'off');
        uistack(stable_patch, 'bottom');
        xline(ax, 0.01, ':', ...
            'Color', [0.22, 0.22, 0.22], 'LineWidth', 1.0, ...
            'HandleVisibility', 'off');
        text(ax, 0.0104, y_limits(1) + 0.055 * diff(y_limits), ...
            '\nu=0.01', 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', ...
            'FontName', 'Times New Roman', 'FontSize', 10.5, ...
            'FontWeight', 'bold', 'Color', [0.18, 0.18, 0.18], ...
            'HandleVisibility', 'off');
        ylim(ax, y_limits);
    end

    legend_ax = axes_handles(1);
    for severity_index = 1:numel(severity_names)
        severity_handles(severity_index) = plot(legend_ax, nan, nan, ...
            severity_styles{severity_index}, 'Color', [0.1, 0.1, 0.1], ...
            'LineWidth', 1.7);
    end
    for frequency_index = 1:numel(frequencies_hz)
        frequency_handles(frequency_index) = plot(legend_ax, nan, nan, ...
            '-', 'Color', palette(frequency_index, :), ...
            'LineWidth', 1.55, 'Marker', markers{frequency_index}, ...
            'MarkerSize', 5.0, 'MarkerFaceColor', 'w');
    end
    legend_handles = [severity_handles(1), severity_handles(2), ...
        frequency_handles];
    legend_labels = {'Mild', 'Severe', '500 Hz', '1000 Hz', ...
        '1500 Hz', '2000 Hz', '4000 Hz'};
    legend(legend_ax, legend_handles, ...
        legend_labels, 'Orientation', 'vertical', 'NumColumns', 1, ...
        'Location', 'northeast', 'FontName', 'Times New Roman', ...
        'FontSize', 10, 'Box', 'on');

    base_name = 'racc_pm_nu_mild_severe_temperature_ac_nsre';
    paths = struct();
    paths.png = fullfile(output_dir, [base_name '.png']);
    paths.pdf = fullfile(output_dir, [base_name '.pdf']);
    paths.fig = fullfile(output_dir, [base_name '.fig']);
    exportgraphics(fig, paths.png, 'Resolution', 300);
    exportgraphics(fig, paths.pdf, 'ContentType', 'vector');

    % Store the editable FIG as a normal interactive MATLAB window.  The
    % raster/vector exports above remain unaffected by these UI properties.
    set(fig, 'Visible', 'on', 'WindowStyle', 'normal', ...
        'MenuBar', 'figure', 'ToolBar', 'figure', 'DockControls', 'on');
    for axes_index = 1:numel(axes_handles)
        enableDefaultInteractivity(axes_handles(axes_index));
        if ~isempty(axes_handles(axes_index).Toolbar)
            axes_handles(axes_index).Toolbar.Visible = 'on';
        end
    end
    drawnow;
    savefig(fig, paths.fig);
    close(fig);
end
