%% AC robustness versus temperature with matched relative ATF-deviation axes
% This response-oriented figure combines the authoritative temperature
% evaluation (including SICER-VAST) with the normalized BZ control-ATF
% mismatch used in the earlier robustness figure. No filters are redesigned
% or re-evaluated here; the script only visualizes saved results.

clear; clc; close all;

this_file = mfilename('fullpath');
response_dir = fileparts(this_file);
project_root = fileparts(response_dir);

result_dir = fullfile(response_dir, 'sicer_vast_temperature_results', ...
    'run_20260828_221103');
result_file = fullfile(result_dir, 'sicer_vast_temperature_results.mat');
deviation_file = fullfile(response_dir, 'reviewer_temperature_results', ...
    'reviewer_temperature_results.mat');
output_dir = fullfile(result_dir, 'figures');

assert(exist(result_file, 'file') == 2, ...
    'Missing SICER-VAST result file: %s', result_file);
assert(exist(deviation_file, 'file') == 2, ...
    'Missing relative-ATF result file: %s', deviation_file);

S = load(result_file, 'experiment');
R = load(deviation_file);
experiment = S.experiment;

temperatures = experiment.temperatures(:);
frequencies = experiment.frequencies(:).';
results = experiment.combined_results;
nominal_temperature = R.nominal_temperature;

target_frequencies = [500, 1000, 3000];
temperature_ticks = 20:25;

% Preserve every algorithm in the authoritative comparison. The final two
% curves are the recent sensor-informed SICER-VAST endpoints.
algorithms = { ...
    'ACC', 'ACC_Reg', 'PM', 'ACC_PM', ...
    'wcRACC', 'NoCT_WCRACC', 'Full_WCRACC', 'POTDC_RACC', ...
    'RPM', 'RACC_PM_Subpro', 'SICER_VAST_V1', 'SICER_VAST_VL'};
display_names = { ...
    'ACC', 'ACC-Reg', 'PM', 'ACC-PM', ...
    'WCRACC', 'NoCT-WCRACC', 'Full-WCRACC', 'POTDC-RACC', ...
    'WCRPM', 'RACC-PM', 'SICER-VAST (V=1)', 'SICER-VAST (V=L)'};

colors = [ ...
    0.0000 0.4470 0.7410; ... % ACC
    0.3500 0.3500 0.3500; ... % ACC-Reg
    0.8510 0.3255 0.0980; ... % PM
    0.3010 0.7450 0.9330; ... % ACC-PM
    0.9294 0.6941 0.1255; ... % WCRACC
    0.1500 0.6000 0.7200; ... % NoCT-WCRACC
    0.6350 0.0780 0.1840; ... % Full-WCRACC
    0.4941 0.1843 0.5569; ... % POTDC-RACC
    0.4667 0.6745 0.1882; ... % WCRPM
    0.8863 0.0667 0.4235; ... % RACC-PM
    0.0000 0.0000 0.0000; ... % SICER-VAST V=1
    0.1000 0.5500 0.5000];    % SICER-VAST V=L
line_styles = {'--','--','--',':','-',':','-.','-','-','-','-','-'};
markers = {'o','x','s','>','^','<','h','v','d','p','o','s'};
line_widths = [1.8 1.8 1.8 1.8 1.8 1.8 1.8 1.8 1.8 2.0 2.1 2.1];
marker_sizes = [5.0 5.5 5.0 5.0 5.0 5.0 5.0 5.0 5.0 7.0 5.0 5.0];

font_name = 'Times New Roman';
tick_font_size = 12;
label_font_size = 12;
legend_font_size = 10.0;
marker_indices = 1:5:numel(temperatures);

fig = figure('Color', 'w', 'Position', [60, 90, 1320, 500], ...
    'MenuBar', 'figure', 'ToolBar', 'figure', ...
    'DockControls', 'on', 'Visible', 'on');
layout = tiledlayout(fig, 1, numel(target_frequencies), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

legend_handles = gobjects(numel(algorithms), 1);
panel_labels = {'(a)', '(b)', '(c)'};
first_axis = gobjects(1);
main_axes = gobjects(1, numel(target_frequencies));
top_axes = gobjects(1, numel(target_frequencies));

for panel = 1:numel(target_frequencies)
    target_frequency = target_frequencies(panel);
    [~, frequency_index] = min(abs(frequencies - target_frequency));
    actual_frequency = frequencies(frequency_index);

    ax = nexttile(layout, panel);
    main_axes(panel) = ax;
    hold(ax, 'on');
    if panel == 1
        first_axis = ax;
    end

    for algorithm_index = 1:numel(algorithms)
        algorithm = algorithms{algorithm_index};
        ac_values = results.(algorithm).AC(:, frequency_index);
        h = plot(ax, temperatures, ac_values, ...
            'Color', colors(algorithm_index, :), ...
            'LineStyle', line_styles{algorithm_index}, ...
            'LineWidth', line_widths(algorithm_index), ...
            'Marker', markers{algorithm_index}, ...
            'MarkerIndices', marker_indices, ...
            'MarkerSize', marker_sizes(algorithm_index), ...
            'MarkerFaceColor', colors(algorithm_index, :));
        if panel == 1
            legend_handles(algorithm_index) = h;
        end
    end

    xline(ax, nominal_temperature, ':', 'Color', [0.12 0.12 0.12], ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
    local_style_axis(ax, font_name, tick_font_size);
    xlim(ax, [temperatures(1), temperatures(end)]);
    xticks(ax, temperature_ticks);
    xlabel(ax, {sprintf('Temperature [^{\\circ}C] (bottom) / \\eta_H (top)'), ...
        sprintf('%s %.0f Hz', panel_labels{panel}, actual_frequency)}, ...
        'Interpreter', 'tex', ...
        'FontName', font_name, 'FontSize', label_font_size, ...
        'FontWeight', 'bold');
    if panel == 1
        ylabel(ax, 'AC [dB]', 'FontName', font_name, ...
            'FontSize', label_font_size, 'FontWeight', 'bold');
    end

    all_ac = nan(numel(temperatures), numel(algorithms));
    for algorithm_index = 1:numel(algorithms)
        all_ac(:, algorithm_index) = ...
            results.(algorithms{algorithm_index}).AC(:, frequency_index);
    end
    lower_limit = 5 * floor(min(all_ac, [], 'all', 'omitnan') / 5);
    upper_limit = 5 * ceil((max(all_ac, [], 'all', 'omitnan') + 2.5) / 5);
    if lower_limit == upper_limit
        upper_limit = lower_limit + 5;
    end
    ylim(ax, [lower_limit, upper_limit]);

    % Overlay a top horizontal axis whose ticks report the normalized BZ
    % control-ATF deviation at the same physical temperatures.
    deviation_frequency_index = find(abs(R.actual_curve_freqs_hz - ...
        actual_frequency) < 1e-9, 1);
    assert(~isempty(deviation_frequency_index), ...
        'No ATF-deviation data for %.0f Hz.', actual_frequency);
    deviation_temperatures = [R.analysis_temperatures(:); nominal_temperature];
    deviation_values = [R.relative_atf_deviation(:, deviation_frequency_index); 0];
    [deviation_temperatures, order] = sort(deviation_temperatures);
    deviation_values = deviation_values(order);
    top_values = interp1(deviation_temperatures, deviation_values, ...
        temperature_ticks, 'linear');

    ax_top = axes(fig, 'Position', ax.Position, 'Color', 'none', ...
        'XAxisLocation', 'top', 'YAxisLocation', 'right', ...
        'XLim', ax.XLim, 'XTick', temperature_ticks, ...
        'XTickLabel', compose('%.2f', top_values), ...
        'YTick', [], 'Box', 'off', 'FontName', font_name, ...
        'FontSize', tick_font_size, 'FontWeight', 'normal', ...
        'LineWidth', 1.0, 'HitTest', 'off');
    ax_top.XMinorTick = 'off';
    ax_top.XGrid = 'on';
    ax_top.GridLineStyle = '-';
    ax_top.GridAlpha = 0.15;
    ax_top.YColor = 'none';
    ax_top.XColor = [0.15 0.15 0.15];
    ax_top.HandleVisibility = 'off';
    top_axes(panel) = ax_top;
    linkaxes([ax, ax_top], 'x');
    local_enable_axis_interactions(ax);
end

legend_handle = legend(first_axis, legend_handles, display_names, ...
    'Orientation', 'horizontal', 'NumColumns', 12, ...
    'Box', 'on', 'FontName', font_name, ...
    'FontSize', legend_font_size);
if isprop(legend_handle, 'ItemTokenSize')
    legend_handle.ItemTokenSize = [12, 7];
end
legend_handle.Layout.Tile = 'north';

% The layout reflows after the shared legend is assigned to the north tile.
% Realign the overlaid top axes so the dual horizontal axes remain exact.
drawnow;
for panel = 1:numel(target_frequencies)
    top_axes(panel).Position = main_axes(panel).Position;
end

output_stem = fullfile(output_dir, ...
    'ac_temperature_relative_atf_sicer_comparison');
exportgraphics(fig, [output_stem '.png'], 'Resolution', 300);
exportgraphics(fig, [output_stem '.pdf'], 'ContentType', 'vector');
savefig(fig, [output_stem '.fig']);

fprintf('Saved editable figure: %s.fig\n', output_stem);
fprintf('Saved vector PDF:     %s.pdf\n', output_stem);
fprintf('Saved PNG preview:    %s.png\n', output_stem);

function local_style_axis(ax, font_name, font_size)
    ax.FontName = font_name;
    ax.FontSize = font_size;
    ax.FontWeight = 'normal';
    ax.LineWidth = 1.0;
    ax.Box = 'on';
    ax.Layer = 'top';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';
    ax.GridLineStyle = '-.';
    ax.MinorGridLineStyle = '-.';
    ax.GridAlpha = 0.20;
    ax.MinorGridAlpha = 0.25;
end

function local_enable_axis_interactions(ax)
    try
        toolbar_handle = axtoolbar(ax, 'default');
        toolbar_handle.Visible = 'on';
    catch
    end
    try
        enableDefaultInteractivity(ax);
    catch
    end
end
