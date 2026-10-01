function output = generate_nominal_frequency_metrics_figure()
%GENERATE_NOMINAL_FREQUENCY_METRICS_FIGURE Plot final nominal metrics.
% The plotted filters are the authoritative rho=10, nu=0.01 filters.  AC,
% AE, and planarity are read from the common monitor-grid evaluation.  To
% remain consistent with the manuscript's nominal audio-oriented protocol,
% NSRE uses the single complex gain estimated at the nominal BZ control
% points and frozen before evaluation on the independent monitor grid.

response_dir = fileparts(mfilename('fullpath'));
run_dir = fullfile(response_dir, 'racc_pm_parameter_comparison_results', ...
    'rho10_nu0p01_all_run_20260827_002614');
result_file = fullfile(run_dir, 'performance_matrices_20260827_021428.mat');
nsre_file = fullfile(run_dir, 'analysis', 'monitor_nsre_postprocess.mat');
out_dir = fullfile(response_dir, 'manuscript_figures');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

loaded = load(result_file, 'val_results', 'algorithms_to_test', ...
    'para', 'freq_params');
nsre_loaded = load(nsre_file, 'results');

algorithms = {'ACC', 'ACC_Reg', 'PM', 'ACC_PM', 'wcRACC', ...
    'NoCT_WCRACC', 'Full_WCRACC', 'POTDC_RACC', 'RPM', ...
    'RACC_PM_Subpro'};
assert(all(ismember(algorithms, loaded.algorithms_to_test)), ...
    'The authoritative result file does not contain all ten algorithms.');

temperatures = loaded.para.temperature_vector_celsius(:);
[temperature_error, nominal_index] = min(abs(temperatures - 22.5));
assert(temperature_error < 1e-9, 'The 22.5-degC realization is missing.');
frequencies = loaded.freq_params.target_freqs(:).';
band_mask = frequencies >= 100 & frequencies <= 4000;
frequencies = frequencies(band_mask);

metric_fields = {'AC', 'NSRE', 'AE', 'Planarity'};
metric_data = cell(1, numel(metric_fields));
for metric_no = 1:numel(metric_fields)
    metric_data{metric_no} = nan(numel(algorithms), numel(frequencies));
end

for algorithm_no = 1:numel(algorithms)
    algorithm = algorithms{algorithm_no};
    metric_data{1}(algorithm_no, :) = ...
        loaded.val_results.(algorithm).AC(nominal_index, band_mask);
    metric_data{2}(algorithm_no, :) = ...
        nsre_loaded.results.aligned_monitor_nsre_db.(algorithm)( ...
        nominal_index, band_mask);
    metric_data{3}(algorithm_no, :) = ...
        loaded.val_results.(algorithm).AE(nominal_index, band_mask);
    metric_data{4}(algorithm_no, :) = ...
        loaded.val_results.(algorithm).Planarity(nominal_index, band_mask);
end

display_names = {'ACC', 'ACC-Reg', 'PM', 'ACC-PM', 'WCRACC', ...
    'NoCT-WCRACC', 'Full-WCRACC', 'POTDC-RACC', 'WCRPM', 'RACC-PM'};
styles = local_styles();

fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [0.5, 0.5, 11.6666, 6.65]);
cleanup = onCleanup(@() close_if_valid(fig));
axes_handles = gobjects(4, 1);
axes_positions = [ ...
    0.065, 0.535, 0.430, 0.390; ...
    0.545, 0.535, 0.430, 0.390; ...
    0.065, 0.075, 0.430, 0.335; ...
    0.545, 0.075, 0.430, 0.335];
titles = {'(a) Acoustic contrast', '(b) Normalized reproduction error', ...
    '(c) Array effort', '(d) Planarity'};
ylabels = {'AC: [dB]', 'NSRE: [dB]', 'AE: [dB]', 'Planarity: [%]'};
ylimits = [0, 105; -65, 5; -30, 65; -10, 105];
fmax_hz = 1630;
legend_handles = gobjects(numel(algorithms), 1);

for metric_no = 1:4
    ax = axes('Parent', fig, 'Position', axes_positions(metric_no, :));
    axes_handles(metric_no) = ax;
    hold(ax, 'on');
    for algorithm_no = 1:numel(algorithms)
        style = styles(algorithm_no);
        line_handle = plot(ax, frequencies, ...
            metric_data{metric_no}(algorithm_no, :), ...
            'Color', style.Color, 'LineStyle', style.LineStyle, ...
            'LineWidth', style.LineWidth, 'Marker', style.Marker, ...
            'MarkerIndices', 1:8:numel(frequencies), ...
            'MarkerSize', style.MarkerSize, ...
            'MarkerFaceColor', style.MarkerFaceColor, ...
            'DisplayName', display_names{algorithm_no});
        if metric_no == 3
            legend_handles(algorithm_no) = line_handle;
        end
    end
    set(ax, 'XScale', 'log');
    xlim(ax, [100, 4000]);
    ylim(ax, ylimits(metric_no, :));
    xticks(ax, [100, 200, 500, 1000, 2000, 4000]);
    xticklabels(ax, {'100', '200', '500', '1000', '2000', '4000'});
    xline(ax, fmax_hz, '--', 'Color', [1.00, 0.00, 0.00], ...
        'LineWidth', 1.80, 'HandleVisibility', 'off');
    if metric_no == 1
        text(ax, fmax_hz * 0.94, ylimits(metric_no, 1) + 4.0, ...
            '$f_{\max}$', 'Interpreter', 'latex', ...
            'Color', [1.00, 0.00, 0.00], 'FontName', ...
            'Times New Roman', 'FontSize', 12.5, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    end
    xlabel(ax, 'Frequency: [Hz]');
    ylabel(ax, ylabels{metric_no});
    title(ax, titles{metric_no});
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 12, ...
        'FontWeight', 'bold', 'LineWidth', 1.0, 'Box', 'on', ...
        'Layer', 'top', 'XGrid', 'on', 'YGrid', 'on', ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'XMinorGrid', 'on', 'YMinorGrid', 'on', ...
        'GridLineStyle', '--', 'MinorGridLineStyle', ':', ...
        'GridAlpha', 0.24, 'MinorGridAlpha', 0.20);
    ax.Title.FontName = 'Times New Roman';
    ax.Title.FontSize = 13;
    ax.Title.FontWeight = 'bold';
    ax.XLabel.FontName = 'Times New Roman';
    ax.XLabel.FontSize = 12.5;
    ax.XLabel.FontWeight = 'bold';
    ax.YLabel.FontName = 'Times New Roman';
    ax.YLabel.FontSize = 12.5;
    ax.YLabel.FontWeight = 'bold';
end

legend_handle = legend(axes_handles(3), legend_handles, display_names, ...
    'NumColumns', 2, 'Orientation', 'vertical', ...
    'FontName', 'Times New Roman', 'FontSize', 7.8, ...
    'Box', 'on', 'Units', 'normalized', 'Location', 'none');
legend_handle.ItemTokenSize = [13, 6];
legend_handle.Position = [0.298, 0.218, 0.168, 0.165];

base = fullfile(out_dir, 'nominal_frequency_metrics_2x2_reference');
drawnow;
savefig(fig, [base '.fig']);
exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector', ...
    'BackgroundColor', 'white');
exportgraphics(fig, [base '.png'], 'Resolution', 300, ...
    'BackgroundColor', 'white');

Algorithm = string(display_names(:));
MeanACdB = mean(metric_data{1}, 2, 'omitnan');
MeanNSREdB = mean(metric_data{2}, 2, 'omitnan');
MeanAEdB = mean(metric_data{3}, 2, 'omitnan');
MeanPlanarityPercent = mean(metric_data{4}, 2, 'omitnan');
summary = table(Algorithm, MeanACdB, MeanNSREdB, MeanAEdB, ...
    MeanPlanarityPercent);
writetable(summary, [base '_band_summary.csv']);
disp(summary);

output = struct('figure_pdf', [base '.pdf'], ...
    'figure_png', [base '.png'], 'figure_fig', [base '.fig'], ...
    'summary_csv', [base '_band_summary.csv'], 'summary', summary);
close(fig);
clear cleanup;
end

function styles = local_styles()
styles = repmat(struct('Color', [0, 0, 0], 'LineStyle', '-', ...
    'LineWidth', 1.35, 'Marker', 'none', 'MarkerSize', 4.0, ...
    'MarkerFaceColor', 'w'), 10, 1);
colors = [ ...
    0.00, 0.45, 0.74; ... % ACC
    0.00, 0.45, 0.74; ... % ACC-Reg
    0.93, 0.69, 0.13; ... % PM
    0.93, 0.69, 0.13; ... % ACC-PM
    0.47, 0.67, 0.19; ... % WCRACC
    0.47, 0.67, 0.19; ... % NoCT-WCRACC
    0.00, 0.55, 0.55; ... % Full-WCRACC
    0.49, 0.18, 0.56; ... % POTDC-RACC
    0.85, 0.33, 0.10; ... % WCRPM
    0.80, 0.10, 0.10];    % RACC-PM
line_styles = {'-', '--', '-', '--', '-', '--', '-.', ':', '-.', '-'};
markers = {'o', 's', '^', 'v', 'd', '>', '<', 'p', 'h', 'o'};
for index = 1:10
    styles(index).Color = colors(index, :);
    styles(index).LineStyle = line_styles{index};
    styles(index).Marker = markers{index};
end
for index = 1:10
    styles(index).LineWidth = 2.10;
    styles(index).MarkerSize = 6.0;
end
styles(10).LineWidth = 2.80;
styles(10).MarkerSize = 7.2;
end

function close_if_valid(fig)
if isgraphics(fig)
    close(fig);
end
end
