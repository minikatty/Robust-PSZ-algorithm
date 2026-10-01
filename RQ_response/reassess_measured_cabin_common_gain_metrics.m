function outputs = reassess_measured_cabin_common_gain_metrics(result_file)
%REASSESS_MEASURED_CABIN_COMMON_GAIN_METRICS Common-gain cabin evaluation.
%   A single complex gain is fitted for every algorithm and frequency on
%   the nominal BZ control points, frozen, and then applied to all 59
%   held-out ATFs. This reproduces the comparison convention used by the
%   earlier manuscript cabin table while preserving native metrics in the
%   production result file.

if nargin < 1 || isempty(result_file)
    error('CabinAligned:MissingResult', ...
        'Provide final_measured_cabin_results.mat.');
end
loaded = load(result_file, 'experiment');
experiment = loaded.experiment;
algorithms = string(experiment.algorithm_fields(:)).';
names = string(experiment.display_names(:)).';
frequencies = double(experiment.frequencies(:)).';
validation = double(experiment.validation_indices(:));
source = double(experiment.config.virtual_source_index);
reference = load(experiment.reference_file, 'ATF_BZ', 'ATF_DZ');

gains = struct();
calibrated_filters = struct();
gain_rows = cell(numel(algorithms), 6);
for algorithm_no = 1:numel(algorithms)
    identifier = char(algorithms(algorithm_no));
    filters = experiment.all_filters.(identifier);
    current_gains = complex(nan(1, numel(frequencies)));
    nominal_nsre = nan(1, numel(frequencies));
    for frequency_no = 1:numel(frequencies)
        HB = reference.ATF_BZ(:, :, frequency_no);
        response = HB * filters(:, frequency_no);
        desired = HB(:, source);
        gain = (response' * desired) / ...
            max(real(response' * response), realmin);
        current_gains(frequency_no) = gain;
        error_signal = gain * response - desired;
        nominal_nsre(frequency_no) = 10 * log10(max( ...
            norm(error_signal)^2 / max(norm(desired)^2, realmin), ...
            realmin));
    end
    gains.(identifier) = current_gains;
    calibrated_filters.(identifier) = filters .* current_gains;
    gain_rows(algorithm_no, :) = {names(algorithm_no), algorithms(algorithm_no), ...
        mean(nominal_nsre), median(abs(current_gains)), ...
        min(abs(current_gains)), max(abs(current_gains))};
end
gain_table = cell2table(gain_rows, 'VariableNames', {'Algorithm', ...
    'AlgorithmField', 'NominalAlignedNSRE_dB', 'MedianGainMagnitude', ...
    'MinimumGainMagnitude', 'MaximumGainMagnitude'});

metric_names = {'AC', 'NSRE', 'AE', 'BZ_SPL_Variance', ...
    'BZ_SPL_Std', 'BZ_SPL_Range'};
results = struct();
for algorithm = algorithms
    identifier = char(algorithm);
    for metric_no = 1:numel(metric_names)
        results.(identifier).(metric_names{metric_no}) = ...
            nan(numel(validation), numel(frequencies));
    end
end

for condition_no = 1:numel(validation)
    current = load(fullfile(experiment.config.data_dir, sprintf( ...
        'ATF_%d.mat', validation(condition_no))), 'ATF_BZ', 'ATF_DZ');
    for algorithm = algorithms
        identifier = char(algorithm);
        filters = calibrated_filters.(identifier);
        for frequency_no = 1:numel(frequencies)
            HB = current.ATF_BZ(:, :, frequency_no);
            HD = current.ATF_DZ(:, :, frequency_no);
            weights = filters(:, frequency_no);
            pressure_B = HB * weights;
            pressure_D = HD * weights;
            desired = HB(:, source);
            results.(identifier).AC(condition_no, frequency_no) = ...
                10 * log10(max(mean(abs(pressure_B).^2), realmin) / ...
                max(mean(abs(pressure_D).^2), realmin));
            results.(identifier).NSRE(condition_no, frequency_no) = ...
                10 * log10(max(norm(pressure_B - desired)^2 / ...
                max(norm(desired)^2, realmin), realmin));
            reference_weights = zeros(size(weights));
            reference_weights(source) = 1;
            RB = HB' * HB;
            generated_energy = max(real(weights' * RB * weights), realmin);
            reference_energy = max(real(reference_weights' * RB * ...
                reference_weights), realmin);
            results.(identifier).AE(condition_no, frequency_no) = ...
                10 * log10(max(real(weights' * weights) * ...
                reference_energy / generated_energy, realmin));
            levels = 20 * log10(max(abs(pressure_B), realmin));
            results.(identifier).BZ_SPL_Variance(condition_no, ...
                frequency_no) = var(levels, 1);
            results.(identifier).BZ_SPL_Std(condition_no, ...
                frequency_no) = std(levels, 1);
            results.(identifier).BZ_SPL_Range(condition_no, ...
                frequency_no) = range(levels);
        end
    end
end

per_frequency = local_per_frequency(results, algorithms, names, frequencies);
band_summary = local_band_summary(results, algorithms, names, frequencies, ...
    experiment.config);
run_tag = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
output_dir = fullfile(experiment.output_dir, 'analysis', ...
    ['common_gain_metrics_' run_tag]);
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
writetable(gain_table, fullfile(output_dir, ...
    'nominal_common_gain_calibration.csv'));
writetable(per_frequency, fullfile(output_dir, ...
    'measured_cabin_per_frequency_summary_common_gain.csv'));
writetable(band_summary, fullfile(output_dir, ...
    'measured_cabin_band_summary_common_gain.csv'));

outputs = struct('source_result_file', result_file, ...
    'output_dir', output_dir, 'gains', gains, ...
    'calibrated_filters', calibrated_filters, 'gain_table', gain_table, ...
    'results', results, 'per_frequency_summary', per_frequency, ...
    'band_summary', band_summary);
save(fullfile(output_dir, 'common_gain_metrics.mat'), 'outputs', '-v7.3');
disp(band_summary(band_summary.Band == "Full_100_4000", :));
fprintf('Common-gain cabin reassessment completed: %s\n', output_dir);
end

function summary = local_per_frequency(results, algorithms, names, frequencies)
row_count = numel(algorithms) * numel(frequencies);
Algorithm = strings(row_count, 1);
AlgorithmField = strings(row_count, 1);
Frequency_Hz = nan(row_count, 1);
AC_Mean = nan(row_count, 1);
NSRE_Mean = nan(row_count, 1);
AE_Mean = nan(row_count, 1);
BZ_SPL_Std_Mean = nan(row_count, 1);
row = 0;
for algorithm_no = 1:numel(algorithms)
    identifier = char(algorithms(algorithm_no));
    for frequency_no = 1:numel(frequencies)
        row = row + 1;
        Algorithm(row) = names(algorithm_no);
        AlgorithmField(row) = algorithms(algorithm_no);
        Frequency_Hz(row) = frequencies(frequency_no);
        AC_Mean(row) = mean(results.(identifier).AC(:, frequency_no));
        NSRE_Mean(row) = mean(results.(identifier).NSRE(:, frequency_no));
        AE_Mean(row) = mean(results.(identifier).AE(:, frequency_no));
        BZ_SPL_Std_Mean(row) = ...
            mean(results.(identifier).BZ_SPL_Std(:, frequency_no));
    end
end
summary = table(Algorithm, AlgorithmField, Frequency_Hz, AC_Mean, ...
    NSRE_Mean, AE_Mean, BZ_SPL_Std_Mean);
end

function summary = local_band_summary(results, algorithms, names, ...
        frequencies, config)
rows = cell(numel(algorithms) * numel(config.band_names), 15);
row = 0;
for band_no = 1:numel(config.band_names)
    limits = config.band_limits_hz(band_no, :);
    selected = frequencies >= limits(1) & frequencies <= limits(2);
    for algorithm_no = 1:numel(algorithms)
        identifier = char(algorithms(algorithm_no));
        ac = mean(results.(identifier).AC(:, selected), 2);
        nsre = mean(results.(identifier).NSRE(:, selected), 2);
        ae = mean(results.(identifier).AE(:, selected), 2);
        spatial_std = mean(results.(identifier).BZ_SPL_Std(:, selected), 2);
        spatial_range = mean(results.(identifier).BZ_SPL_Range(:, selected), 2);
        variance = mean(results.(identifier).BZ_SPL_Variance(:, selected), 2);
        row = row + 1;
        rows(row, :) = {names(algorithm_no), algorithms(algorithm_no), ...
            config.band_names(band_no), limits(1), limits(2), mean(ac), ...
            min(ac), mean(nsre), max(nsre), mean(ae), mean(variance), ...
            mean(spatial_std), max(spatial_std), mean(spatial_range), ...
            max(spatial_range)};
    end
end
summary = cell2table(rows, 'VariableNames', {'Algorithm', ...
    'AlgorithmField', 'Band', 'BandStart_Hz', 'BandEnd_Hz', ...
    'AC_Mean_dB', 'AC_WorstCondition_dB', 'NSRE_Mean_dB', ...
    'NSRE_WorstCondition_dB', 'AE_Mean_dB', ...
    'BZ_SPL_Variance_Mean_dB2', 'BZ_SPL_Std_Mean_dB', ...
    'BZ_SPL_Std_WorstCondition_dB', 'BZ_SPL_Range_Mean_dB', ...
    'BZ_SPL_Range_WorstCondition_dB'});
end
