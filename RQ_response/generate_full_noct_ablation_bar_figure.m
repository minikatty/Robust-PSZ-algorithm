function output = generate_full_noct_ablation_bar_figure()
%GENERATE_FULL_NOCT_ABLATION_BAR_FIGURE Focused cross-term ablation figure.
%   Reuses the authoritative common-scale frequency-band summary without
%   redesigning filters.  No existing figure or result asset is overwritten.

    response_dir = fileparts(mfilename('fullpath'));
    project_root = fileparts(response_dir);
    source_dir = fullfile(response_dir, ...
        'racc_pm_parameter_comparison_results', ...
        'all_algorithms_frequency_band_analysis');
    source_csv = fullfile(source_dir, ...
        'all_algorithms_frequency_band_summary.csv');
    if ~isfile(source_csv)
        error('CrossTermFigure:MissingSummary', ...
            'Required summary file not found: %s', source_csv);
    end

    summary = readtable(source_csv, 'TextType', 'string');
    setup_id = "rho10_nu0p01";
    algorithms = ["NoCT_WCRACC", "Full_WCRACC"];
    display_names = ["NoCT-WCRACC", "Full-WCRACC"];
    bands = ["band1", "band2", "band3"];
    band_titles = ["Band I: 100--500 Hz", ...
        "Band II: 550--1600 Hz", "Band III: 1650--4000 Hz"];
    scope_fields = ["NominalACdB", "NearMeanACdB", ...
        "SweepMeanACdB", "SweepWorstACdB"];
    scope_labels = ["Nominal", "Near mean", ...
        "Sweep mean", "Sweep worst"];

    values = nan(numel(bands), numel(scope_fields), numel(algorithms));
    rows = cell(numel(bands) * numel(scope_fields) * ...
        numel(algorithms), 5);
    row_index = 0;
    for b = 1:numel(bands)
        for s = 1:numel(scope_fields)
            for a = 1:numel(algorithms)
                selected = summary.Setup == setup_id & ...
                    summary.Band == bands(b) & ...
                    summary.Algorithm == algorithms(a);
                if nnz(selected) ~= 1
                    error('CrossTermFigure:Selection', ...
                        'Expected one row for %s/%s/%s.', ...
                        setup_id, bands(b), algorithms(a));
                end
                value = summary.(scope_fields(s))(selected);
                values(b, s, a) = value;
                row_index = row_index + 1;
                rows(row_index, :) = {char(bands(b)), char(band_titles(b)), ...
                    char(scope_labels(s)), char(display_names(a)), value};
            end
        end
    end

    run_stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    output_dir = fullfile(response_dir, ...
        'cross_term_ablation_figures', ['run_' run_stamp]);
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    plotted_values = cell2table(rows, 'VariableNames', ...
        {'Band', 'BandLabel', 'PerturbationScope', 'Algorithm', 'ACdB'});
    plotted_values.ACdB = round(plotted_values.ACdB, 2);
    writetable(plotted_values, fullfile(output_dir, ...
        'full_noct_ablation_values.csv'));

    font_name = 'Times New Roman';
    colors = [68, 119, 170; 238, 102, 119] / 255; % colorblind-safe
    fig = figure('Color', 'w', 'ToolBar', 'none', ...
        'Position', [80, 120, 1500, 500]);
    layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', ...
        'Padding', 'compact');
    layout.OuterPosition = [0.02, 0.02, 0.96, 0.84];

    for b = 1:numel(bands)
        ax = nexttile(layout);
        disableDefaultInteractivity(ax);
        ax.Toolbar.Visible = 'off';
        panel_values = squeeze(values(b, :, :));
        bars = bar(ax, panel_values, 'grouped', 'BarWidth', 0.82);
        hold(ax, 'on');
        for a = 1:numel(bars)
            bars(a).FaceColor = colors(a, :);
            bars(a).EdgeColor = [0.20, 0.20, 0.20];
            bars(a).LineWidth = 0.65;
        end
        ax.XTick = 1:numel(scope_labels);
        ax.XTickLabel = scope_labels;
        ax.XTickLabelRotation = 20;
        ax.YGrid = 'on';
        ax.YMinorGrid = 'on';
        ax.XGrid = 'off';
        ax.GridAlpha = 0.18;
        ax.MinorGridAlpha = 0.08;
        box(ax, 'on');
        title(ax, band_titles(b), 'FontName', font_name, ...
            'FontWeight', 'bold', 'FontSize', 11);
        if b == 1
            ylabel(ax, 'Acoustic contrast [dB]', 'FontName', font_name, ...
                'FontWeight', 'bold');
        end
        set(ax, 'FontName', font_name, 'FontSize', 9, ...
            'LineWidth', 0.75);

        y_max = max(panel_values, [], 'all');
        y_pad = max(1.5, 0.12 * y_max);
        ylim(ax, [0, y_max + y_pad]);
        label_offset = max(0.35, 0.012 * y_max);
        for a = 1:numel(bars)
            x_positions = bars(a).XEndPoints;
            y_positions = bars(a).YEndPoints;
            labels = compose('%.2f', bars(a).YData);
            text(ax, x_positions, y_positions + label_offset, labels, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'FontName', font_name, 'FontSize', 8, ...
                'Color', [0.10, 0.10, 0.10]);
        end

        hold(ax, 'off');
    end

    annotation(fig, 'textbox', [0.19, 0.945, 0.62, 0.04], ...
        'String', ...
        'Full versus truncated ATF-uncertainty expansion (\nu = 0.01)', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontName', font_name, 'FontWeight', 'bold', 'FontSize', 12, ...
        'Interpreter', 'tex', 'EdgeColor', 'none');
    annotation(fig, 'rectangle', [0.385, 0.895, 0.015, 0.022], ...
        'FaceColor', colors(1, :), 'EdgeColor', [0.20, 0.20, 0.20]);
    annotation(fig, 'textbox', [0.404, 0.888, 0.105, 0.035], ...
        'String', display_names(1), 'VerticalAlignment', 'middle', ...
        'FontName', font_name, 'FontSize', 9, 'EdgeColor', 'none');
    annotation(fig, 'rectangle', [0.515, 0.895, 0.015, 0.022], ...
        'FaceColor', colors(2, :), 'EdgeColor', [0.20, 0.20, 0.20]);
    annotation(fig, 'textbox', [0.534, 0.888, 0.105, 0.035], ...
        'String', display_names(2), 'VerticalAlignment', 'middle', ...
        'FontName', font_name, 'FontSize', 9, 'EdgeColor', 'none');

    base = fullfile(output_dir, 'full_noct_band_scope_ablation');
    png_path = [base '.png'];
    pdf_path = [base '.pdf'];
    fig_path = [base '.fig'];
    exportgraphics(fig, png_path, 'Resolution', 300);
    exportgraphics(fig, pdf_path, 'ContentType', 'vector');
    savefig(fig, fig_path);
    close(fig);

    output = struct('project_root', project_root, ...
        'source_csv', source_csv, 'output_dir', output_dir, ...
        'values', plotted_values, 'png', png_path, ...
        'pdf', pdf_path, 'fig', fig_path);
    save(fullfile(output_dir, 'figure_output.mat'), 'output');
    fprintf('Focused cross-term ablation figure: %s\n', output_dir);
end
