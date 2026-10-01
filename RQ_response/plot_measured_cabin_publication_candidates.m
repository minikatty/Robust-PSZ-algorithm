function outputs = plot_measured_cabin_publication_candidates(run_dir, source_file)
%PLOT_MEASURED_CABIN_PUBLICATION_CANDIDATES Publication-style cabin plots.
%   Creates an AC-only candidate and a four-metric candidate without
%   replacing the preserved automatic diagnostic figures.

response_dir = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(run_dir)
    run_dir = fullfile(response_dir, 'final_measured_cabin_results', ...
        'run_20260829_152136');
end
if ~isfolder(run_dir)
    error('MeasuredCabinPlot:MissingRunDirectory', ...
        'Run directory does not exist: %s', run_dir);
end
if nargin < 2 || isempty(source_file)
    source_file = fullfile(run_dir, ...
        'measured_cabin_per_frequency_summary.csv');
end
if ~isfile(source_file)
    error('MeasuredCabinPlot:MissingSourceFile', ...
        'Per-frequency source file does not exist: %s', source_file);
end
output_dir = fullfile(run_dir, 'figures', 'publication_candidates');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

data = readtable(source_file, 'TextType', 'string');
data = data(data.Frequency_Hz >= 100 & data.Frequency_Hz <= 4000, :);

fields = ["ACC", "ACC_Reg", "PM", "ACC_PM", "wcRACC", ...
    "NoCT_WCRACC", "Full_WCRACC", "POTDC_RACC", "RPM", ...
    "RACC_PM_Subpro", "M_ACC"];
names = ["ACC", "ACC-Reg", "PM", "ACC-PM", "WCRACC", ...
    "NoCT-WCRACC", "Full-WCRACC", "POTDC-RACC", "WCRPM", ...
    "RACC-PM", "M-ACC"];

styles = local_styles();
font_name = 'Times New Roman';

% Candidate 1: engineering-focused AC-only figure.
fig_ac = figure('Color', 'w', 'Visible', 'off', 'Units', 'centimeters', ...
    'Position', [2, 2, 19, 11]);
ax = axes(fig_ac);
handles = local_plot_metric(ax, data, fields, names, styles, 'AC_Mean');
title(ax, 'Measured-cabin acoustic contrast', 'FontWeight', 'bold');
xlabel(ax, 'Frequency [Hz]', 'FontWeight', 'bold');
ylabel(ax, 'AC: [dB]', 'FontWeight', 'bold');
xlim(ax, [100, 4000]);
ylim(ax, [-10, 40]);
local_style_axes(ax, font_name);
lgd = legend(ax, handles, names, 'NumColumns', 4, ...
    'Location', 'southoutside', 'FontName', font_name, ...
    'FontSize', 8.5, 'Box', 'on');
lgd.ItemTokenSize = [23, 11];
ac_base = fullfile(output_dir, 'measured_cabin_ac_publication_candidate');
local_save(fig_ac, ac_base);
close(fig_ac);

% Candidate 2: compact engineering/profile figure.
fig_profile = figure('Color', 'w', 'Visible', 'off', ...
    'Units', 'centimeters', 'Position', [2, 2, 26, 16.5]);
layout = tiledlayout(fig_profile, 2, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
metric_fields = {'AC_Mean', 'NSRE_Mean', 'AE_Mean', 'BZ_SPL_Std_Mean'};
titles = {'(a) Acoustic contrast', ...
    '(b) Normalized reproduction error', ...
    '(c) Array effort', ...
    '(d) BZ SPL standard deviation'};
ylabels = {'AC: [dB]', 'NSRE: [dB]', 'AE: [dB]', ...
    'BZ SPL std.: [dB]'};
ylimits = {[-5, 38], [-32, 14], [-35, 48], [0.5, 8.0]};
handles = gobjects(numel(fields), 1);
for panel = 1:4
    ax = nexttile(layout, panel);
    current = local_plot_metric(ax, data, fields, names, styles, ...
        metric_fields{panel});
    if panel == 1
        handles = current;
    end
    title(ax, titles{panel}, 'FontWeight', 'bold');
    xlabel(ax, 'Frequency [Hz]', 'FontWeight', 'bold');
    ylabel(ax, ylabels{panel}, 'FontWeight', 'bold');
    xlim(ax, [100, 4000]);
    ylim(ax, ylimits{panel});
    local_style_axes(ax, font_name);
end
lgd = legend(handles, names, 'Orientation', 'horizontal', ...
    'NumColumns', 6, 'FontName', font_name, 'FontSize', 8.5, ...
    'Box', 'on');
lgd.Layout.Tile = 'north';
lgd.ItemTokenSize = [22, 10];
profile_base = fullfile(output_dir, ...
    'measured_cabin_ac_nsre_ae_uniformity_publication_candidate');
local_save(fig_profile, profile_base);
close(fig_profile);

% Candidate 3: full-width 1-by-4 layout with the legend inside panel (a).
fig_wide = figure('Color', 'w', 'Visible', 'off', ...
    'Units', 'centimeters', 'Position', [1, 2, 40, 11.8]);
layout = tiledlayout(fig_wide, 1, 4, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
handles = gobjects(numel(fields), 1);
wide_ylimits = {[-5, 38], [-32, 14], [-35, 48], [0.5, 8.0]};
for panel = 1:4
    ax = nexttile(layout, panel);
    current = local_plot_metric(ax, data, fields, names, styles, ...
        metric_fields{panel});
    if panel == 1
        handles = current;
    end
    title(ax, titles{panel}, 'FontWeight', 'bold');
    xlabel(ax, 'Frequency [Hz]', 'FontWeight', 'bold');
    ylabel(ax, ylabels{panel}, 'FontWeight', 'bold');
    xlim(ax, [100, 4000]);
    ylim(ax, wide_ylimits{panel});
    local_style_axes(ax, font_name);
end
first_axis = nexttile(layout, 1);
lgd = legend(first_axis, handles, names, 'NumColumns', 3, ...
    'Location', 'northeast', 'FontName', font_name, ...
    'FontSize', 6.2, 'Box', 'on');
lgd.ItemTokenSize = [10, 6];
drawnow;
lgd.Units = 'normalized';
lgd.Position = [0.102, 0.688, 0.135, 0.190];
wide_base = fullfile(output_dir, ...
    'measured_cabin_ac_nsre_ae_uniformity_publication_candidate_1x4');
local_save(fig_wide, wide_base);
close(fig_wide);

outputs = struct();
outputs.output_dir = output_dir;
outputs.ac_only = ac_base;
outputs.profile = profile_base;
outputs.wide_profile = wide_base;
end

function handles = local_plot_metric(ax, data, fields, names, styles, metric)
hold(ax, 'on');
handles = gobjects(numel(fields), 1);
for algorithm_no = 1:numel(fields)
    subset = data(data.AlgorithmField == fields(algorithm_no), :);
    subset = sortrows(subset, 'Frequency_Hz');
    style = styles(algorithm_no);
    marker_start = 1 + mod(algorithm_no - 1, 5);
    marker_indices = marker_start:7:height(subset);
    handles(algorithm_no) = plot(ax, subset.Frequency_Hz, ...
        subset.(metric), 'Color', style.Color, ...
        'LineStyle', style.LineStyle, 'LineWidth', style.LineWidth, ...
        'Marker', style.Marker, 'MarkerIndices', marker_indices, ...
        'MarkerSize', style.MarkerSize, 'MarkerFaceColor', ...
        style.MarkerFaceColor, 'DisplayName', names(algorithm_no));
end
end

function styles = local_styles()
styles = repmat(struct('Color', [0, 0, 0], 'LineStyle', '-', ...
    'LineWidth', 1.35, 'Marker', 'none', 'MarkerSize', 4.0, ...
    'MarkerFaceColor', 'none'), 11, 1);
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
    0.80, 0.10, 0.10; ... % RACC-PM
    0.30, 0.30, 0.30];    % M-ACC
line_styles = {'-', '--', '-', '--', '-', '--', '-.', ':', '-.', '-', '--'};
markers = {'o', 's', '^', 'v', 'd', '>', '<', 'p', 'h', 'o', 'x'};
for index = 1:11
    styles(index).Color = colors(index, :);
    styles(index).LineStyle = line_styles{index};
    styles(index).Marker = markers{index};
    if strcmp(markers{index}, 'x')
        styles(index).MarkerFaceColor = 'none';
    else
        styles(index).MarkerFaceColor = 'w';
    end
end
styles(10).LineWidth = 2.2;
styles(10).MarkerSize = 5.2;
end

function local_style_axes(ax, font_name)
set(ax, 'FontName', font_name, 'FontSize', 10.5, ...
    'FontWeight', 'bold', 'LineWidth', 1.0, 'Box', 'on', ...
    'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'on', ...
    'YMinorGrid', 'on', 'GridLineStyle', '--', ...
    'MinorGridLineStyle', ':', 'GridAlpha', 0.25, ...
    'MinorGridAlpha', 0.14, 'TickDir', 'in', 'Layer', 'top');
ax.Toolbar.Visible = 'off';
end

function local_save(fig, base)
exportgraphics(fig, [base, '.png'], 'Resolution', 400, ...
    'BackgroundColor', 'white');
exportgraphics(fig, [base, '.pdf'], 'ContentType', 'vector', ...
    'BackgroundColor', 'white');
set(fig, 'Visible', 'on');
savefig(fig, [base, '.fig']);
end
