function outputs = generate_measured_cabin_summary_table(result_file, band_summary_file)
%GENERATE_MEASURED_CABIN_SUMMARY_TABLE Export the cabin manuscript table.
%   The table uses the Full_100_4000 rows saved by the production pipeline.

if nargin < 1 || isempty(result_file)
    error('CabinTable:MissingResult', ...
        'Provide final_measured_cabin_results.mat from the intended run.');
end
loaded = load(result_file, 'experiment');
experiment = loaded.experiment;
if nargin >= 2 && ~isempty(band_summary_file)
    summary = readtable(band_summary_file, 'TextType', 'string');
else
    summary = experiment.band_summary;
end
table_data = summary(summary.Band == "Full_100_4000", ...
    {'Algorithm', 'AlgorithmField', 'AC_Mean_dB', ...
    'AC_WorstCondition_dB', 'NSRE_Mean_dB', 'AE_Mean_dB', ...
    'BZ_SPL_Std_Mean_dB'});

output_dir = fullfile(experiment.output_dir, 'manuscript_assets');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
csv_file = fullfile(output_dir, 'measured_cabin_table_100_4000.csv');
writetable(table_data, csv_file);

tex_file = fullfile(output_dir, 'measured_cabin_table_100_4000.tex');
file_id = fopen(tex_file, 'w');
if file_id < 0
    error('CabinTable:OpenFailed', 'Could not create %s.', tex_file);
end
cleanup = onCleanup(@() fclose(file_id));

[~, best_ac] = max(table_data.AC_Mean_dB);
[~, best_worst] = max(table_data.AC_WorstCondition_dB);
[~, best_nsre] = min(table_data.NSRE_Mean_dB);
[~, best_ae] = min(table_data.AE_Mean_dB);
[~, best_std] = min(table_data.BZ_SPL_Std_Mean_dB);

fprintf(file_id, '%% Generated from %s\n', result_file);
fprintf(file_id, '\\begin{table}[t]\n');
fprintf(file_id, '\\centering\n');
fprintf(file_id, ['\\caption{Quantitative performance over the 59 ', ...
    'held-out measured-cabin realizations (100--4000~Hz).}\n']);
fprintf(file_id, '\\label{tab:measured-cabin-summary}\n');
fprintf(file_id, '\\scriptsize\n');
fprintf(file_id, '\\setlength{\\tabcolsep}{1.8pt}\n');
fprintf(file_id, '\\renewcommand{\\arraystretch}{1.02}\n');
fprintf(file_id, '\\begin{tabular}{@{}lccccc@{}}\n');
fprintf(file_id, '\\toprule\n');
fprintf(file_id, ['Algorithm & Mean AC (dB) $\\uparrow$ & ', ...
    'Worst AC (dB) $\\uparrow$ & NSRE (dB) $\\downarrow$ & ', ...
    'AE (dB) $\\downarrow$ & BZ std. (dB) $\\downarrow$ \\\\\n']);
fprintf(file_id, '\\midrule\n');
for row = 1:height(table_data)
    name = char(table_data.Algorithm(row));
    values = [table_data.AC_Mean_dB(row), ...
        table_data.AC_WorstCondition_dB(row), ...
        table_data.NSRE_Mean_dB(row), table_data.AE_Mean_dB(row), ...
        table_data.BZ_SPL_Std_Mean_dB(row)];
    best = [row == best_ac, row == best_worst, row == best_nsre, ...
        row == best_ae, row == best_std];
    cells = strings(1, 5);
    for column = 1:5
        cells(column) = sprintf('%.2f', values(column));
        if best(column)
            cells(column) = "\textbf{" + cells(column) + "}";
        end
    end
    if table_data.AlgorithmField(row) == "RACC_PM_Subpro"
        name = ['\textbf{', name, '}'];
    end
    fprintf(file_id, '%s & %s & %s & %s & %s & %s \\\\\n', ...
        name, cells(1), cells(2), cells(3), cells(4), cells(5));
end
fprintf(file_id, '\\bottomrule\n');
fprintf(file_id, '\\end{tabular}\n');
fprintf(file_id, '\\par\\vspace{0.5mm}\\scriptsize\\raggedright\n');
fprintf(file_id, ['NSRE uses a per-frequency complex gain estimated from ', ...
    'the nominal BZ control points and then frozen for all held-out ', ...
    'conditions. AC, AE, and BZ standard deviation are unaffected by ', ...
    'this common gain.\\par\n']);
fprintf(file_id, '\\end{table}\n');

outputs = struct('output_dir', output_dir, 'csv_file', csv_file, ...
    'tex_file', tex_file, 'table', table_data);
end
