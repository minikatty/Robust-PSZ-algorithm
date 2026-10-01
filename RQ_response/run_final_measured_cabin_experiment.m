function outputs = run_final_measured_cabin_experiment(user_config)
%RUN_FINAL_MEASURED_CABIN_EXPERIMENT Final parameter-consistent cabin test.
%   Every filter bank is redesigned from ATF_2 only. The unchanged filters
%   are then evaluated on the remaining 59 measured ATF realizations. The
%   60 records are treated as operating conditions, not independent human
%   subjects. No held-out realization is used for parameter selection.

if nargin < 1 || isempty(user_config)
    user_config = struct();
end

response_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(response_dir);
config = local_default_config(project_root, response_dir);
config = local_apply_overrides(config, user_config);

if ~strcmpi(pwd, project_root)
    error('FinalCabin:WrongWorkingDirectory', ...
        'Run this function from the project root: %s.', project_root);
end
addpath(fullfile(project_root, 'RQ_response'));
addpath(genpath(fullfile(project_root, 'src')));
addpath(genpath(fullfile(project_root, 'data')));
addpath(fullfile(project_root, 'src', 'evaluations'), '-begin');

resolved_evaluator = which('evaluate_performance');
expected_evaluator = fullfile(project_root, 'src', 'evaluations', ...
    'evaluate_performance.m');
if ~strcmpi(resolved_evaluator, expected_evaluator)
    error('FinalCabin:EvaluatorResolution', ...
        'Unexpected evaluate_performance resolution: %s', ...
        resolved_evaluator);
end
resolved_racc_pm = which('RACC_PM_Sub');
expected_racc_pm = fullfile(project_root, 'src', 'algorithm', ...
    'RACC_PM_Sub.m');
if ~strcmpi(resolved_racc_pm, expected_racc_pm)
    error('FinalCabin:RACCPMResolution', ...
        'Unexpected RACC_PM_Sub resolution: %s', resolved_racc_pm);
end

if isempty(config.resume_run_dir)
    run_id = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    output_dir = fullfile(config.output_root, ['run_', run_id]);
else
    output_dir = config.resume_run_dir;
    [~, run_name] = fileparts(output_dir);
    run_id = erase(run_name, 'run_');
end
figure_dir = fullfile(output_dir, 'figures');
design_checkpoint_dir = fullfile(output_dir, 'design_checkpoints');
evaluation_checkpoint_dir = fullfile(output_dir, ...
    'evaluation_checkpoints');
full_checkpoint_dir = fullfile(output_dir, ...
    'full_wcracc_frequency_checkpoints');
racc_pm_checkpoint_dir = fullfile(output_dir, ...
    'racc_pm_frequency_checkpoints_full_ball_alpha');
local_ensure_directory(output_dir);
local_ensure_directory(figure_dir);
local_ensure_directory(design_checkpoint_dir);
local_ensure_directory(evaluation_checkpoint_dir);
local_ensure_directory(full_checkpoint_dir);
local_ensure_directory(racc_pm_checkpoint_dir);

diary('off');
diary(fullfile(output_dir, 'matlab_diary.log'));
diary_cleanup = onCleanup(@() diary('off'));
total_wall_start = tic;
total_cpu_start = cputime;
run_started = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
fprintf('Final measured-cabin experiment started: %s\n', ...
    string(run_started));
fprintf('Output directory: %s\n', output_dir);
fprintf('Resolved evaluator: %s\n', resolved_evaluator);
fprintf('Resolved RACC-PM implementation: %s\n', resolved_racc_pm);

start_marker = fullfile(output_dir, 'RUN_STARTED.txt');
if ~exist(start_marker, 'file')
    marker_id = fopen(start_marker, 'w');
    marker_cleanup = onCleanup(@() local_safe_fclose(marker_id));
    fprintf(marker_id, 'Started: %s\n', string(run_started));
    fprintf(marker_id, 'MATLAB PID: %d\n', feature('getpid'));
    fprintf(marker_id, 'Script: %s\n', mfilename('fullpath'));
    clear marker_cleanup;
end

metadata = load(fullfile(config.data_dir, 'Nominal_ATF.mat'), ...
    'freq_params');
freq_params = metadata.freq_params;
frequencies = freq_params.target_freqs(:).';
measurement_indices = local_measurement_indices(config.data_dir);
if ~ismember(config.reference_index, measurement_indices)
    error('FinalCabin:MissingReference', ...
        'ATF_%d.mat was not found.', config.reference_index);
end
validation_indices = setdiff(measurement_indices, ...
    config.reference_index, 'stable');
reference_file = fullfile(config.data_dir, sprintf('ATF_%d.mat', ...
    config.reference_index));
reference = load(reference_file, 'ATF_BZ', 'ATF_DZ');
local_validate_reference(reference, frequencies, config);

ATF_BZ_ctrl = reference.ATF_BZ;
ATF_DZ_ctrl = reference.ATF_DZ;
ATF_desired = squeeze(ATF_BZ_ctrl(:, ...
    config.virtual_source_index, :));
if isvector(ATF_desired)
    ATF_desired = reshape(ATF_desired, size(ATF_BZ_ctrl, 1), []);
end

eta_B = local_relative_radii(ATF_BZ_ctrl, config.nu);
eta_D = local_relative_radii(ATF_DZ_ctrl, config.nu);
freq_params.epsilon = struct('B', eta_B, 'D', eta_D);
freq_params.gamma = struct('B', eta_B.^2, 'D', eta_D.^2);
freq_params.scale = 1;
freq_params.fre_idces = 1:numel(frequencies);
freq_params.virtual_src_idx = config.virtual_source_index;
freq_params.kappa = config.kappa;
freq_params.full_racc = struct( ...
    'relative_eta', config.nu, ...
    'diagonal_loading_ratio', config.diagonal_loading_ratio, ...
    'formulation', 'decomposed', ...
    'solution_method', 'normalized', ...
    'solver', 'mosek', ...
    'mosek_threads', config.mosek_threads_per_worker, ...
    'parallel_workers', config.full_parallel_workers, ...
    'checkpoint_dir', full_checkpoint_dir, ...
    'resume_from_checkpoints', true, ...
    'retry_unsolved_checkpoints', true, ...
    'normalization_multipliers', [0.1, 1, 5, 10, 0.01], ...
    'fallback_to_bisection', true, ...
    'recovery_method', 'paper', ...
    'randomization_trials', config.full_randomization_trials, ...
    'random_seed', 20260820, ...
    'refinement_iterations', 300, ...
    'refinement_seeds', 8);
freq_params.potdc = struct( ...
    'iterations', 20, ...
    'tolerance', 1e-8, ...
    'solver', 'mosek', ...
    'mosek_threads', config.mosek_threads_per_worker, ...
    'parallel_workers', config.full_parallel_workers);
freq_params.racc_pm = struct( ...
    'implementation', 'RACC_PM_Sub', ...
    'parallel_workers', config.racc_pm_parallel_workers, ...
    'mosek_threads', config.racc_pm_mosek_threads_per_worker, ...
    'checkpoint_dir', racc_pm_checkpoint_dir);

noct_reference = NoCT_WCRACC(ATF_BZ_ctrl, ATF_DZ_ctrl, freq_params);
% With four control microphones and 36 loudspeakers, nominal NoCT AC can
% be arbitrarily inflated by exact dark-zone nulls. Use the contrast that
% the same nominal NoCT filter attains under the common full ATF ball.
% This remains a zero-margin, ATF-2-control-only target and gives RACC-PM
% a physically matched feasible reference instead of a nominal deep-null
% artifact.
alpha_AC = noct_reference.full_atf_worst_case_alpha(:);
if any(~isfinite(alpha_AC)) || any(alpha_AC <= 0)
    error('FinalCabin:InvalidControlDerivedAlpha', ...
        'The NoCT full-ATF-ball control contrast is invalid.');
end
alpha_source = ['full ATF-ball robust contrast attained by the ', ...
    'nominal NoCT-WCRACC control filter'];
freq_params.mu = config.mu;
freq_params.rho = config.rho;
freq_params.alpha = alpha_AC;
freq_params.para_gamma = config.mu + config.rho * alpha_AC;
fprintf(['Common radii: eta_Z(f)=%.4g||H_Z(f)||_F; ', ...
    'alpha uses the NoCT control filter full-ball contrast with ', ...
    '0 dB margin.\n'], config.nu);

[macc_selection, selected_macc_alpha, macc_filters, ...
    macc_diagnostics, macc_loading_check, macc_selection_wall_seconds] = ...
    local_select_macc(reference, frequencies, freq_params, config);
writetable(macc_selection, fullfile(output_dir, ...
    'macc_nominal_alpha_selection.csv'));
writetable(macc_loading_check, fullfile(output_dir, ...
    'macc_paper_loading_numerical_check.csv'));

algorithm_fields = string(config.algorithms(:));
display_names = arrayfun(@local_display_name, algorithm_fields);
design_data = struct();
design_data.ATF_BZ.ctrl = ATF_BZ_ctrl;
design_data.ATF_DZ.ctrl = ATF_DZ_ctrl;
design_data.ATF_desired = ATF_desired;

signature = struct();
signature.protocol_version = 1;
signature.reference_index = config.reference_index;
signature.virtual_source_index = config.virtual_source_index;
signature.frequencies = frequencies;
signature.algorithms = algorithm_fields;
signature.nu = config.nu;
signature.rho = config.rho;
signature.mu = config.mu;
signature.kappa = config.kappa;
signature.diagonal_loading_ratio = config.diagonal_loading_ratio;
signature.eta_B = eta_B;
signature.eta_D = eta_D;
signature.alpha_AC = alpha_AC;
signature.alpha_source = alpha_source;
signature.macc_alpha = selected_macc_alpha;
signature.macc_beta = config.macc_beta;
signature.macc_beta_mode = config.macc_beta_mode;
signature.reuse_filter_result_file = string(config.reuse_filter_result_file);
signature.reuse_filter_algorithms = string(config.reuse_filter_algorithms(:));

reuse_bank = local_load_reuse_filter_bank(config, frequencies, ...
    algorithm_fields);

all_filters = struct();
algorithm_details = struct();
algorithm_wall_seconds = struct();
algorithm_cpu_seconds = struct();
for algorithm_no = 1:numel(algorithm_fields)
    algorithm = char(algorithm_fields(algorithm_no));
    checkpoint_file = fullfile(design_checkpoint_dir, sprintf( ...
        '%02d_%s.mat', algorithm_no, algorithm));
    fprintf('====== Designing %s (%d/%d) ======\n', algorithm, ...
        algorithm_no, numel(algorithm_fields));
    if exist(checkpoint_file, 'file')
        saved = load(checkpoint_file, 'checkpoint');
        local_validate_design_checkpoint(saved.checkpoint, algorithm, ...
            signature);
        all_filters.(algorithm) = saved.checkpoint.filter;
        algorithm_details.(algorithm) = saved.checkpoint.detail;
        algorithm_wall_seconds.(algorithm) = ...
            saved.checkpoint.wall_seconds;
        algorithm_cpu_seconds.(algorithm) = ...
            saved.checkpoint.cpu_seconds;
        fprintf('%s restored from its design checkpoint.\n', algorithm);
        continue;
    end

    if isfield(reuse_bank.filters, algorithm)
        filter = reuse_bank.filters.(algorithm);
        detail = reuse_bank.details.(algorithm);
        wall_seconds = reuse_bank.wall_seconds.(algorithm);
        cpu_seconds = reuse_bank.cpu_seconds.(algorithm);
        local_validate_filter(filter, algorithm, reference, frequencies);
        all_filters.(algorithm) = filter;
        algorithm_details.(algorithm) = detail;
        algorithm_wall_seconds.(algorithm) = wall_seconds;
        algorithm_cpu_seconds.(algorithm) = cpu_seconds;
        checkpoint = struct('algorithm', algorithm, 'signature', signature, ...
            'filter', filter, 'detail', detail, ...
            'wall_seconds', wall_seconds, 'cpu_seconds', cpu_seconds, ...
            'completed_at', datetime('now'), ...
            'reused_from', string(config.reuse_filter_result_file));
        save(checkpoint_file, 'checkpoint', '-v7.3');
        fprintf('%s reused from %s.\n', algorithm, ...
            config.reuse_filter_result_file);
        continue;
    end

    wall_start = tic;
    cpu_start = cputime;
    if strcmp(algorithm, 'RACC_PM_Subpro')
        detail = run_racc_pm_sub_parallel_checkpointed( ...
            ATF_BZ_ctrl, ATF_DZ_ctrl, ATF_desired, freq_params, ...
            racc_pm_checkpoint_dir, config.racc_pm_parallel_workers);
        filter = detail.w;
    elseif strcmp(algorithm, 'M_ACC')
        filter = macc_filters;
        detail = macc_diagnostics;
        wall_seconds = detail.design_wall_seconds;
        cpu_seconds = detail.design_cpu_seconds;
    else
        [filter, detail] = design_filters(algorithm, design_data, ...
            freq_params);
    end
    if ~strcmp(algorithm, 'M_ACC')
        wall_seconds = toc(wall_start);
        cpu_seconds = cputime - cpu_start;
    end
    local_validate_filter(filter, algorithm, reference, frequencies);
    local_validate_solver_status(algorithm, detail);

    all_filters.(algorithm) = filter;
    algorithm_details.(algorithm) = detail;
    algorithm_wall_seconds.(algorithm) = wall_seconds;
    algorithm_cpu_seconds.(algorithm) = cpu_seconds;
    checkpoint = struct('algorithm', algorithm, 'signature', signature, ...
        'filter', filter, 'detail', detail, ...
        'wall_seconds', wall_seconds, 'cpu_seconds', cpu_seconds, ...
        'completed_at', datetime('now'));
    save(checkpoint_file, 'checkpoint', '-v7.3');
    fprintf('%s design completed in %.3f s.\n', algorithm, wall_seconds);
end

[calibrated_filters, calibration_table, calibration_gains] = ...
    local_calibrate_filters(all_filters, algorithm_fields, reference, ...
    config.virtual_source_index, frequencies);
writetable(calibration_table, fullfile(output_dir, ...
    'nominal_complex_gain_calibration.csv'));

val_results = local_allocate_results(algorithm_fields, ...
    numel(validation_indices), numel(frequencies));
fprintf('Evaluating %d frozen filter banks on %d held-out ATFs.\n', ...
    numel(algorithm_fields), numel(validation_indices));
for validation_no = 1:numel(validation_indices)
    measurement_index = validation_indices(validation_no);
    checkpoint_file = fullfile(evaluation_checkpoint_dir, sprintf( ...
        'ATF_%03d.mat', measurement_index));
    if exist(checkpoint_file, 'file')
        saved = load(checkpoint_file, 'checkpoint');
        local_validate_evaluation_checkpoint(saved.checkpoint, ...
            measurement_index, signature, calibration_gains);
        condition_results = saved.checkpoint.condition_results;
    else
        current = load(fullfile(config.data_dir, sprintf( ...
            'ATF_%d.mat', measurement_index)), 'ATF_BZ', 'ATF_DZ');
        local_validate_operating_atf(current, reference, measurement_index);
        condition_results = struct();
        for algorithm_no = 1:numel(algorithm_fields)
            algorithm = char(algorithm_fields(algorithm_no));
            condition_results.(algorithm) = local_evaluate_filter_bank( ...
                calibrated_filters.(algorithm), current.ATF_BZ, ...
                current.ATF_DZ, config.virtual_source_index, ...
                local_is_target_referenced(algorithm));
        end
        checkpoint = struct('measurement_index', measurement_index, ...
            'signature', signature, ...
            'calibration_gains', calibration_gains, ...
            'condition_results', condition_results, ...
            'completed_at', datetime('now'));
        save(checkpoint_file, 'checkpoint', '-v7');
    end
    val_results = local_insert_condition(val_results, condition_results, ...
        algorithm_fields, validation_no);
    if mod(validation_no, 10) == 0 || ...
            validation_no == numel(validation_indices)
        fprintf('Held-out ATFs completed: %d/%d\n', validation_no, ...
            numel(validation_indices));
    end
end

per_frequency_summary = local_per_frequency_summary(val_results, ...
    algorithm_fields, display_names, frequencies);
writetable(per_frequency_summary, fullfile(output_dir, ...
    'measured_cabin_per_frequency_summary.csv'));
band_summary = local_band_summary(val_results, algorithm_fields, ...
    display_names, frequencies, config, algorithm_wall_seconds, ...
    algorithm_cpu_seconds);
writetable(band_summary, fullfile(output_dir, ...
    'measured_cabin_band_summary.csv'));
paired_differences = local_paired_differences(val_results, ...
    algorithm_fields, display_names, frequencies, config);
writetable(paired_differences, fullfile(output_dir, ...
    'paired_condition_differences.csv'));
reference_diagnostic = local_reference_diagnostic(config, ...
    measurement_indices, frequencies);
writetable(reference_diagnostic, fullfile(output_dir, ...
    'reference_representativeness.csv'));

figure_paths = struct();
figure_paths.frequency_metrics = local_create_frequency_figure( ...
    per_frequency_summary, algorithm_fields, display_names, config, ...
    figure_dir);
figure_paths.band_metrics = local_create_band_figure( ...
    val_results, algorithm_fields, display_names, frequencies, config, ...
    figure_dir);
figure_paths.macc_selection = local_create_alpha_figure( ...
    macc_selection, figure_dir);

total_wall_seconds = toc(total_wall_start);
total_cpu_seconds = cputime - total_cpu_start;
environment = struct('matlab_version', version, ...
    'computer', computer, 'resolved_evaluator', resolved_evaluator, ...
    'resolved_racc_pm', resolved_racc_pm, ...
    'full_parallel_workers', config.full_parallel_workers, ...
    'mosek_threads_per_worker', config.mosek_threads_per_worker, ...
    'racc_pm_parallel_workers', config.racc_pm_parallel_workers, ...
    'racc_pm_mosek_threads_per_worker', ...
    config.racc_pm_mosek_threads_per_worker);

experiment = struct();
experiment.run_id = run_id;
experiment.output_dir = output_dir;
experiment.config = config;
experiment.signature = signature;
experiment.reference_file = reference_file;
experiment.measurement_indices = measurement_indices;
experiment.validation_indices = validation_indices;
experiment.frequencies = frequencies;
experiment.algorithm_fields = algorithm_fields;
experiment.display_names = display_names;
experiment.freq_params = freq_params;
experiment.all_filters = all_filters;
experiment.calibrated_filters = calibrated_filters;
experiment.calibration_gains = calibration_gains;
experiment.calibration_table = calibration_table;
experiment.algorithm_details = algorithm_details;
experiment.algorithm_wall_seconds = algorithm_wall_seconds;
experiment.algorithm_cpu_seconds = algorithm_cpu_seconds;
experiment.val_results = val_results;
experiment.macc_selection = macc_selection;
experiment.selected_macc_alpha = selected_macc_alpha;
experiment.macc_loading_check = macc_loading_check;
experiment.macc_selection_wall_seconds = macc_selection_wall_seconds;
experiment.per_frequency_summary = per_frequency_summary;
experiment.band_summary = band_summary;
experiment.paired_differences = paired_differences;
experiment.reference_diagnostic = reference_diagnostic;
experiment.figure_paths = figure_paths;
experiment.environment = environment;
experiment.total_wall_seconds = total_wall_seconds;
experiment.total_cpu_seconds = total_cpu_seconds;
experiment.run_started = run_started;
experiment.run_completed = datetime('now', ...
    'Format', 'yyyy-MM-dd HH:mm:ss');

result_file = fullfile(output_dir, ...
    'final_measured_cabin_results.mat');
save(result_file, 'experiment', '-v7.3');
marker_id = fopen(fullfile(output_dir, 'RUN_COMPLETE.txt'), 'w');
marker_cleanup = onCleanup(@() local_safe_fclose(marker_id));
fprintf(marker_id, 'Completed: %s\n', string(experiment.run_completed));
fprintf(marker_id, 'Result: %s\n', result_file);
fprintf(marker_id, 'Reference: ATF_%d.mat\n', config.reference_index);
fprintf(marker_id, 'Held-out operating conditions: %d\n', ...
    numel(validation_indices));
fprintf(marker_id, 'Algorithms: %s\n', strjoin(algorithm_fields, ', '));
fprintf(marker_id, 'nu: %.9g\n', config.nu);
fprintf(marker_id, 'rho: %.9g\n', config.rho);
fprintf(marker_id, 'mu: %.9g\n', config.mu);
fprintf(marker_id, 'RACC-PM: no-suffix RACC_PM_Sub\n');
fprintf(marker_id, 'RACC-PM alpha source: %s\n', alpha_source);
fprintf(marker_id, 'M-ACC alpha: %.9g\n', selected_macc_alpha);
fprintf(marker_id, 'Total wall seconds: %.6f\n', total_wall_seconds);
fprintf(marker_id, 'Total CPU seconds: %.6f\n', total_cpu_seconds);
clear marker_cleanup;

fprintf('\nFinal measured-cabin band summary:\n');
disp(band_summary(band_summary.Band == "Full_100_4000", :));
fprintf('Result: %s\n', result_file);
clear diary_cleanup;

outputs = struct('output_dir', output_dir, 'result_file', result_file, ...
    'band_summary', band_summary, ...
    'per_frequency_summary', per_frequency_summary, ...
    'paired_differences', paired_differences, ...
    'reference_diagnostic', reference_diagnostic, ...
    'figure_paths', figure_paths);
end

function config = local_default_config(project_root, response_dir)
config = struct();
config.data_dir = fullfile(project_root, 'data', 'Cabin_Measurements');
config.output_root = fullfile(response_dir, ...
    'final_measured_cabin_results');
config.resume_run_dir = '';
config.reference_index = 2;
config.virtual_source_index = 20;
config.nu = 0.01;
config.rho = 10;
config.mu = 1;
config.kappa = 0.5;
config.diagonal_loading_ratio = 1e-6;
config.full_parallel_workers = 2;
config.mosek_threads_per_worker = 2;
config.racc_pm_parallel_workers = 1;
config.racc_pm_mosek_threads_per_worker = [];
config.full_randomization_trials = 1000;
config.algorithms = {'ACC', 'ACC_Reg', 'PM', 'ACC_PM', 'wcRACC', ...
    'NoCT_WCRACC', 'Full_WCRACC', 'POTDC_RACC', 'RPM', ...
    'RACC_PM_Subpro', 'M_ACC'};
config.macc_alpha_candidates = [0, 1e-4, 1e-3, 1e-2, 1e-1];
config.macc_beta = 1e-6;
config.macc_beta_mode = 'relative_dark_maxeig';
config.macc_paper_beta = 1e-14;
config.macc_selection_band_hz = [100, 2000];
config.analysis_band_hz = [100, 4000];
config.reuse_filter_result_file = '';
config.reuse_filter_algorithms = {};
config.band_names = ["Paper_100_2000", "Full_100_4000", ...
    "Low_100_500", "Mid_600_1600", "High_1700_4000"];
config.band_limits_hz = [100, 2000; 100, 4000; 100, 500; ...
    600, 1600; 1700, 4000];
config.paired_references = {'RACC_PM_Subpro', 'ACC_Reg'};
config.figure_font = 'Times New Roman';
end

function config = local_apply_overrides(config, overrides)
if ~isstruct(overrides) || ~isscalar(overrides)
    error('FinalCabin:InvalidConfig', ...
        'user_config must be a scalar structure.');
end
fields = fieldnames(overrides);
for field_no = 1:numel(fields)
    if ~isfield(config, fields{field_no})
        error('FinalCabin:UnknownConfig', ...
            'Unknown configuration field: %s.', fields{field_no});
    end
    config.(fields{field_no}) = overrides.(fields{field_no});
end
validateattributes(config.nu, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
validateattributes(config.rho, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
validateattributes(config.mu, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
if ~isempty(config.reuse_filter_result_file) && ...
        ~isfile(config.reuse_filter_result_file)
    error('FinalCabin:MissingReuseResult', ...
        'Reuse result file does not exist: %s', ...
        config.reuse_filter_result_file);
end
if isstring(config.reuse_filter_algorithms)
    config.reuse_filter_algorithms = cellstr( ...
        config.reuse_filter_algorithms(:));
end
if ~iscell(config.reuse_filter_algorithms)
    error('FinalCabin:InvalidReuseAlgorithms', ...
        'reuse_filter_algorithms must be a cell array or string array.');
end
target_referenced = {'PM', 'ACC_PM', 'RPM', 'RACC_PM_Subpro'};
invalid_reuse = intersect(config.reuse_filter_algorithms, ...
    target_referenced, 'stable');
if ~isempty(invalid_reuse)
    error('FinalCabin:TargetFilterReuse', ...
        'Target-referenced filters must be redesigned: %s.', ...
        strjoin(invalid_reuse, ', '));
end
end

function reuse = local_load_reuse_filter_bank(config, frequencies, algorithms)
reuse = struct('filters', struct(), 'details', struct(), ...
    'wall_seconds', struct(), 'cpu_seconds', struct());
if isempty(config.reuse_filter_result_file) || ...
        isempty(config.reuse_filter_algorithms)
    return;
end
loaded = load(config.reuse_filter_result_file, 'experiment');
if ~isfield(loaded, 'experiment')
    error('FinalCabin:InvalidReuseResult', ...
        'Reuse result does not contain experiment.');
end
source = loaded.experiment;
required = {'config', 'frequencies', 'all_filters', ...
    'algorithm_details', 'algorithm_wall_seconds', ...
    'algorithm_cpu_seconds'};
for field_no = 1:numel(required)
    if ~isfield(source, required{field_no})
        error('FinalCabin:InvalidReuseResult', ...
            'Reuse experiment lacks %s.', required{field_no});
    end
end
if source.config.reference_index ~= config.reference_index || ...
        ~isequal(source.frequencies(:), frequencies(:)) || ...
        source.config.nu ~= config.nu || ...
        source.config.diagonal_loading_ratio ~= ...
        config.diagonal_loading_ratio || ...
        source.config.macc_beta ~= config.macc_beta || ...
        ~strcmp(source.config.macc_beta_mode, config.macc_beta_mode)
    error('FinalCabin:IncompatibleReuseResult', ...
        'The reuse result does not match the current design protocol.');
end
requested = string(config.reuse_filter_algorithms(:));
available = string(algorithms(:));
if any(~ismember(requested, available))
    error('FinalCabin:UnknownReuseAlgorithm', ...
        'A requested reuse algorithm is absent from the current run.');
end
for algorithm_no = 1:numel(requested)
    algorithm = char(requested(algorithm_no));
    if ~isfield(source.all_filters, algorithm)
        error('FinalCabin:MissingReuseFilter', ...
            'Reuse result lacks filter %s.', algorithm);
    end
    reuse.filters.(algorithm) = source.all_filters.(algorithm);
    reuse.details.(algorithm) = source.algorithm_details.(algorithm);
    reuse.wall_seconds.(algorithm) = ...
        source.algorithm_wall_seconds.(algorithm);
    reuse.cpu_seconds.(algorithm) = ...
        source.algorithm_cpu_seconds.(algorithm);
end
fprintf('Reusing %d source-independent filter banks from %s.\n', ...
    numel(requested), config.reuse_filter_result_file);
end

function indices = local_measurement_indices(data_dir)
listing = dir(fullfile(data_dir, 'ATF_*.mat'));
indices = nan(1, numel(listing));
keep = false(1, numel(listing));
for file_no = 1:numel(listing)
    token = regexp(listing(file_no).name, '^ATF_(\d+)\.mat$', ...
        'tokens', 'once');
    if ~isempty(token)
        indices(file_no) = str2double(token{1});
        keep(file_no) = true;
    end
end
indices = sort(indices(keep));
if numel(indices) ~= 60 || ~isequal(indices, 1:60)
    error('FinalCabin:UnexpectedMeasurementSet', ...
        'Expected ATF_1.mat through ATF_60.mat; found %s.', ...
        mat2str(indices));
end
end

function local_validate_reference(reference, frequencies, config)
if ~isequal(size(reference.ATF_BZ), size(reference.ATF_DZ)) || ...
        size(reference.ATF_BZ, 1) ~= 4 || ...
        size(reference.ATF_BZ, 2) ~= 36 || ...
        size(reference.ATF_BZ, 3) ~= numel(frequencies)
    error('FinalCabin:UnexpectedReferenceDimensions', ...
        'Expected four points, 36 loudspeakers, and the saved grid.');
end
if config.virtual_source_index > size(reference.ATF_BZ, 2)
    error('FinalCabin:InvalidVirtualSource', ...
        'virtual_source_index exceeds the loudspeaker count.');
end
if any(~isfinite(reference.ATF_BZ(:))) || ...
        any(~isfinite(reference.ATF_DZ(:)))
    error('FinalCabin:NonfiniteReference', ...
        'The reference ATFs contain nonfinite values.');
end
end

function local_validate_operating_atf(current, reference, index)
if ~isequal(size(current.ATF_BZ), size(reference.ATF_BZ)) || ...
        ~isequal(size(current.ATF_DZ), size(reference.ATF_DZ)) || ...
        any(~isfinite(current.ATF_BZ(:))) || ...
        any(~isfinite(current.ATF_DZ(:)))
    error('FinalCabin:InvalidOperatingATF', ...
        'ATF_%d has invalid dimensions or values.', index);
end
end

function radii = local_relative_radii(H, relative_level)
radii = nan(size(H, 3), 1);
for frequency_no = 1:size(H, 3)
    radii(frequency_no) = relative_level * ...
        norm(H(:, :, frequency_no), 'fro');
end
end

function [selection, selected_alpha, selected_filters, ...
        selected_diagnostics, loading_check, total_wall_seconds] = ...
        local_select_macc(reference, frequencies, freq_params, config)
total_start = tic;
candidate_count = numel(config.macc_alpha_candidates);
AC_dB = nan(candidate_count, 1);
BZ_SPL_Variance_dB2 = nan(candidate_count, 1);
BZ_SPL_Std_dB = nan(candidate_count, 1);
BZ_SPL_Range_dB = nan(candidate_count, 1);
CandidateWallSeconds = nan(candidate_count, 1);
MinimumDenominatorRcond = nan(candidate_count, 1);
MaximumStationarityResidual = nan(candidate_count, 1);
selection_band = frequencies >= config.macc_selection_band_hz(1) & ...
    frequencies <= config.macc_selection_band_hz(2);
for candidate_no = 1:candidate_count
    para = freq_params;
    para.macc.alpha = config.macc_alpha_candidates(candidate_no);
    para.macc.beta = config.macc_beta;
    para.macc.beta_mode = config.macc_beta_mode;
    candidate_start = tic;
    [candidate_filter, diagnostics] = M_ACC(reference.ATF_BZ, ...
        reference.ATF_DZ, para);
    CandidateWallSeconds(candidate_no) = toc(candidate_start);
    nominal = local_unscaled_spatial_metrics(candidate_filter, ...
        reference.ATF_BZ, reference.ATF_DZ);
    AC_dB(candidate_no) = mean(nominal.AC(selection_band));
    BZ_SPL_Variance_dB2(candidate_no) = ...
        mean(nominal.BZ_SPL_Variance(selection_band));
    BZ_SPL_Std_dB(candidate_no) = ...
        mean(nominal.BZ_SPL_Std(selection_band));
    BZ_SPL_Range_dB(candidate_no) = ...
        mean(nominal.BZ_SPL_Range(selection_band));
    MinimumDenominatorRcond(candidate_no) = ...
        min(diagnostics.denominator_rcond);
    MaximumStationarityResidual(candidate_no) = ...
        max(diagnostics.stationarity_residual);
end
normalized_ac_loss = local_normalize(max(AC_dB) - AC_dB);
normalized_uniformity = (local_normalize(BZ_SPL_Std_dB) + ...
    local_normalize(BZ_SPL_Range_dB)) / 2;
SelectionScore = hypot(normalized_ac_loss, normalized_uniformity);
[~, selected_index] = min(SelectionScore);
Selected = false(candidate_count, 1);
Selected(selected_index) = true;
Alpha = config.macc_alpha_candidates(:);
selection = table(Alpha, AC_dB, BZ_SPL_Variance_dB2, ...
    BZ_SPL_Std_dB, BZ_SPL_Range_dB, normalized_ac_loss, ...
    normalized_uniformity, SelectionScore, Selected, ...
    CandidateWallSeconds, MinimumDenominatorRcond, ...
    MaximumStationarityResidual);
selected_alpha = Alpha(selected_index);

selected_para = freq_params;
selected_para.macc.alpha = selected_alpha;
selected_para.macc.beta = config.macc_beta;
selected_para.macc.beta_mode = config.macc_beta_mode;
selected_start = tic;
selected_cpu_start = cputime;
[selected_filters, selected_diagnostics] = M_ACC( ...
    reference.ATF_BZ, reference.ATF_DZ, selected_para);
selected_diagnostics.design_wall_seconds = toc(selected_start);
selected_diagnostics.design_cpu_seconds = cputime - selected_cpu_start;
selected_diagnostics.selected_alpha = selected_alpha;
selected_diagnostics.selection_rule = [ ...
    'minimum Euclidean distance to normalized ATF-2 control ideal ', ...
    '(maximum AC, minimum BZ SPL standard deviation and range)'];

paper_para = freq_params;
paper_para.macc.alpha = selected_alpha;
paper_para.macc.beta = config.macc_paper_beta;
paper_para.macc.beta_mode = 'absolute';
paper_status = "Solved";
try
    [~, paper_diagnostics] = M_ACC(reference.ATF_BZ, ...
        reference.ATF_DZ, paper_para);
    paper_min_rcond = min(paper_diagnostics.denominator_rcond);
    paper_max_residual = max(paper_diagnostics.stationarity_residual);
catch exception
    paper_status = "Failed: " + string(exception.identifier);
    paper_min_rcond = NaN;
    paper_max_residual = NaN;
end
loading_check = table(["Paper absolute loading"; ...
    "Selected scale-normalized loading"], ...
    [config.macc_paper_beta; config.macc_beta], ...
    ["absolute"; string(config.macc_beta_mode)], ...
    [paper_min_rcond; min(selected_diagnostics.denominator_rcond)], ...
    [paper_max_residual; ...
    max(selected_diagnostics.stationarity_residual)], ...
    [paper_status; "Solved"], 'VariableNames', {'Setting', ...
    'RequestedBeta', 'BetaMode', 'MinimumDenominatorRcond', ...
    'MaximumStationarityResidual', 'Status'});
total_wall_seconds = toc(total_start);
end

function normalized = local_normalize(values)
span = max(values) - min(values);
if span <= eps(max([1; abs(values(:))]))
    normalized = zeros(size(values));
else
    normalized = (values - min(values)) / span;
end
end

function metrics = local_unscaled_spatial_metrics(filters, HB, HD)
frequency_count = size(filters, 2);
metrics.AC = nan(1, frequency_count);
metrics.BZ_SPL_Variance = nan(1, frequency_count);
metrics.BZ_SPL_Std = nan(1, frequency_count);
metrics.BZ_SPL_Range = nan(1, frequency_count);
for frequency_no = 1:frequency_count
    pressure_B = HB(:, :, frequency_no) * filters(:, frequency_no);
    pressure_D = HD(:, :, frequency_no) * filters(:, frequency_no);
    metrics.AC(frequency_no) = 10 * log10( ...
        max(mean(abs(pressure_B).^2), realmin) / ...
        max(mean(abs(pressure_D).^2), realmin));
    levels = 20 * log10(max(abs(pressure_B), realmin));
    metrics.BZ_SPL_Variance(frequency_no) = var(levels, 1);
    metrics.BZ_SPL_Std(frequency_no) = std(levels, 1);
    metrics.BZ_SPL_Range(frequency_no) = range(levels);
end
end

function local_validate_filter(filter, algorithm, reference, frequencies)
expected = [size(reference.ATF_BZ, 2), numel(frequencies)];
if ~isequal(size(filter), expected) || any(~isfinite(filter(:)))
    error('FinalCabin:InvalidFilter', ...
        '%s returned an invalid filter bank.', algorithm);
end
end

function local_validate_solver_status(algorithm, detail)
if ~ismember(string(algorithm), ...
        ["Full_WCRACC", "POTDC_RACC", "RACC_PM_Subpro"])
    return;
end
if ~isstruct(detail) || ~isfield(detail, 'status')
    error('FinalCabin:MissingSolverStatus', ...
        '%s did not return solver statuses.', algorithm);
end
statuses = string(detail.status);
solved = statuses == "Solved" | statuses == "Inaccurate/Solved";
if any(~solved)
    error('FinalCabin:UnsolvedFrequency', ...
        '%s has unsolved frequency indices: %s.', algorithm, ...
        mat2str(find(~solved)));
end
end

function local_validate_design_checkpoint(checkpoint, algorithm, signature)
required = {'algorithm', 'signature', 'filter', 'detail', ...
    'wall_seconds', 'cpu_seconds'};
for field_no = 1:numel(required)
    if ~isfield(checkpoint, required{field_no})
        error('FinalCabin:IncompleteDesignCheckpoint', ...
            '%s checkpoint lacks %s.', algorithm, required{field_no});
    end
end
saved_signature = checkpoint.signature;
expected_signature = signature;
% Alpha affects only RACC-PM. This permits a scientifically valid resume
% of the already completed alpha-invariant designs after correcting the
% measured-cabin alpha source from nominal deep-null AC to the matched
% full-ATF-ball NoCT control contrast.
if ~strcmp(algorithm, 'RACC_PM_Subpro')
    saved_signature.alpha_AC = expected_signature.alpha_AC;
    saved_signature.alpha_source = expected_signature.alpha_source;
end
if ~strcmp(checkpoint.algorithm, algorithm) || ...
        ~isequaln(saved_signature, expected_signature)
    error('FinalCabin:DesignCheckpointMismatch', ...
        '%s checkpoint belongs to a different protocol.', algorithm);
end
end

function [calibrated, table_out, gains] = local_calibrate_filters( ...
        filters, algorithms, reference, virtual_source_index, frequencies)
calibrated = struct();
gains = struct();
rows = cell(numel(algorithms), 7);
for algorithm_no = 1:numel(algorithms)
    algorithm = char(algorithms(algorithm_no));
    raw = filters.(algorithm);
    current_gains = complex(zeros(1, numel(frequencies)));
    calibrated_bank = complex(zeros(size(raw)));
    nominal_nsre = nan(1, numel(frequencies));
    is_target_referenced = local_is_target_referenced(algorithm);
    if is_target_referenced
        preparation_mode = "Native target solution";
    else
        preparation_mode = "Contrast-only control-point LS anchor";
    end
    for frequency_no = 1:numel(frequencies)
        HB = reference.ATF_BZ(:, :, frequency_no);
        response = HB * raw(:, frequency_no);
        desired = HB(:, virtual_source_index);
        if is_target_referenced
            gain = 1;
        else
            gain = (response' * desired) / ...
                max(real(response' * response), realmin);
        end
        current_gains(frequency_no) = gain;
        calibrated_bank(:, frequency_no) = raw(:, frequency_no) * gain;
        error_signal = HB * calibrated_bank(:, frequency_no) - desired;
        nominal_nsre(frequency_no) = 10 * log10(max( ...
            norm(error_signal)^2 / max(norm(desired)^2, realmin), ...
            realmin));
    end
    calibrated.(algorithm) = calibrated_bank;
    gains.(algorithm) = current_gains;
    rows(algorithm_no, :) = {algorithm, char(preparation_mode), ...
        mean(nominal_nsre), ...
        median(abs(current_gains)), min(abs(current_gains)), ...
        max(abs(current_gains)), ...
        median(abs(rad2deg(angle(current_gains))))};
end
table_out = cell2table(rows, 'VariableNames', {'Algorithm', ...
    'PreparationMode', 'NominalPreparedError_dB', 'MedianGainMagnitude', ...
    'MinimumGainMagnitude', 'MaximumGainMagnitude', ...
    'MedianAbsGainPhase_deg'});
end

function results = local_allocate_results(algorithms, conditions, frequencies)
metric_fields = local_metric_fields();
results = struct();
for algorithm = algorithms(:).'
    for metric_no = 1:numel(metric_fields)
        results.(char(algorithm)).(metric_fields{metric_no}) = ...
            nan(conditions, frequencies);
    end
end
end

function fields = local_metric_fields()
fields = {'AC', 'NSRE', 'ContrastCalibratedError', 'AE', 'BZ_SPL_Variance', ...
    'BZ_SPL_Std', 'BZ_SPL_Range'};
end

function metrics = local_evaluate_filter_bank(filters, HB_all, HD_all, ...
        virtual_source_index, is_target_referenced)
frequency_count = size(filters, 2);
fields = local_metric_fields();
for field_no = 1:numel(fields)
    metrics.(fields{field_no}) = nan(1, frequency_count);
end
for frequency_no = 1:frequency_count
    HB = HB_all(:, :, frequency_no);
    HD = HD_all(:, :, frequency_no);
    weights = filters(:, frequency_no);
    pressure_B = HB * weights;
    pressure_D = HD * weights;
    desired = HB(:, virtual_source_index);
    metrics.AC(frequency_no) = 10 * log10( ...
        max(mean(abs(pressure_B).^2), realmin) / ...
        max(mean(abs(pressure_D).^2), realmin));
    error_db = 10 * log10(max(norm(pressure_B - desired)^2 / ...
        max(norm(desired)^2, realmin), realmin));
    if is_target_referenced
        metrics.NSRE(frequency_no) = error_db;
    else
        metrics.ContrastCalibratedError(frequency_no) = error_db;
    end
    reference_weights = zeros(size(weights));
    reference_weights(virtual_source_index) = 1;
    RB = HB' * HB;
    generated_energy = max(real(weights' * RB * weights), realmin);
    reference_energy = max(real(reference_weights' * RB * ...
        reference_weights), realmin);
    metrics.AE(frequency_no) = 10 * log10(max( ...
        real(weights' * weights) * reference_energy / ...
        generated_energy, realmin));
    levels = 20 * log10(max(abs(pressure_B), realmin));
    metrics.BZ_SPL_Variance(frequency_no) = var(levels, 1);
    metrics.BZ_SPL_Std(frequency_no) = std(levels, 1);
    metrics.BZ_SPL_Range(frequency_no) = range(levels);
end
end

function tf = local_is_target_referenced(algorithm)
tf = ismember(string(algorithm), ["PM", "ACC_PM", "WCRPM", "RPM", ...
    "RACC_PM", "RACC_PM_Sub", "RACC_PM_Subpro", "RACC_PM_GLS"]);
end

function local_validate_evaluation_checkpoint(checkpoint, index, ...
        signature, gains)
if checkpoint.measurement_index ~= index || ...
        ~isequaln(checkpoint.signature, signature) || ...
        ~isequaln(checkpoint.calibration_gains, gains)
    error('FinalCabin:EvaluationCheckpointMismatch', ...
        'ATF_%d evaluation checkpoint belongs to another protocol.', index);
end
end

function results = local_insert_condition(results, condition, ...
        algorithms, row)
fields = local_metric_fields();
for algorithm = algorithms(:).'
    name = char(algorithm);
    for field_no = 1:numel(fields)
        field = fields{field_no};
        results.(name).(field)(row, :) = condition.(name).(field);
    end
end
end

function summary = local_per_frequency_summary(results, algorithms, ...
        display_names, frequencies)
fields = local_metric_fields();
row_count = numel(algorithms) * numel(frequencies);
Algorithm = strings(row_count, 1);
AlgorithmField = strings(row_count, 1);
Frequency_Hz = nan(row_count, 1);
data = struct();
for field_no = 1:numel(fields)
    field = fields{field_no};
    data.([field, '_Mean']) = nan(row_count, 1);
    data.([field, '_P10']) = nan(row_count, 1);
    data.([field, '_P90']) = nan(row_count, 1);
end
row = 0;
for algorithm_no = 1:numel(algorithms)
    algorithm = char(algorithms(algorithm_no));
    for frequency_no = 1:numel(frequencies)
        row = row + 1;
        Algorithm(row) = display_names(algorithm_no);
        AlgorithmField(row) = algorithms(algorithm_no);
        Frequency_Hz(row) = frequencies(frequency_no);
        for field_no = 1:numel(fields)
            field = fields{field_no};
            values = results.(algorithm).(field)(:, frequency_no);
            data.([field, '_Mean'])(row) = mean(values, 'omitnan');
            data.([field, '_P10'])(row) = prctile(values, 10);
            data.([field, '_P90'])(row) = prctile(values, 90);
        end
    end
end
summary = table(Algorithm, AlgorithmField, Frequency_Hz);
for field_no = 1:numel(fields)
    field = fields{field_no};
    summary.([field, '_Mean']) = data.([field, '_Mean']);
    summary.([field, '_P10']) = data.([field, '_P10']);
    summary.([field, '_P90']) = data.([field, '_P90']);
end
end

function summary = local_band_summary(results, algorithms, display_names, ...
        frequencies, config, wall_seconds, cpu_seconds)
rows = cell(numel(algorithms) * numel(config.band_names), 26);
row = 0;
for band_no = 1:numel(config.band_names)
    limits = config.band_limits_hz(band_no, :);
    band = frequencies >= limits(1) & frequencies <= limits(2);
    for algorithm_no = 1:numel(algorithms)
        row = row + 1;
        algorithm = char(algorithms(algorithm_no));
        ac = mean(results.(algorithm).AC(:, band), 2, 'omitnan');
        nsre = mean(results.(algorithm).NSRE(:, band), 2, 'omitnan');
        ae = mean(results.(algorithm).AE(:, band), 2, 'omitnan');
        variance = mean(results.(algorithm).BZ_SPL_Variance(:, band), ...
            2, 'omitnan');
        spatial_std = mean(results.(algorithm).BZ_SPL_Std(:, band), ...
            2, 'omitnan');
        spatial_range = mean(results.(algorithm).BZ_SPL_Range(:, band), ...
            2, 'omitnan');
        rows(row, :) = {display_names(algorithm_no), ...
            algorithms(algorithm_no), config.band_names(band_no), ...
            limits(1), limits(2), mean(ac), prctile(ac, 10), ...
            prctile(ac, 90), min(ac), mean(nsre), prctile(nsre, 10), ...
            prctile(nsre, 90), max(nsre), mean(ae), ...
            prctile(ae, 90), mean(variance), prctile(variance, 90), ...
            max(variance), mean(spatial_std), ...
            prctile(spatial_std, 90), max(spatial_std), ...
            mean(spatial_range), prctile(spatial_range, 90), ...
            max(spatial_range), wall_seconds.(algorithm), ...
            cpu_seconds.(algorithm)};
    end
end
summary = cell2table(rows, 'VariableNames', {'Algorithm', ...
    'AlgorithmField', 'Band', 'BandStart_Hz', 'BandEnd_Hz', ...
    'AC_Mean_dB', 'AC_P10_dB', 'AC_P90_dB', ...
    'AC_WorstCondition_dB', 'NSRE_Mean_dB', 'NSRE_P10_dB', ...
    'NSRE_P90_dB', 'NSRE_WorstCondition_dB', 'AE_Mean_dB', ...
    'AE_P90_dB', 'BZ_SPL_Variance_Mean_dB2', ...
    'BZ_SPL_Variance_P90_dB2', ...
    'BZ_SPL_Variance_WorstCondition_dB2', ...
    'BZ_SPL_Std_Mean_dB', 'BZ_SPL_Std_P90_dB', ...
    'BZ_SPL_Std_WorstCondition_dB', 'BZ_SPL_Range_Mean_dB', ...
    'BZ_SPL_Range_P90_dB', 'BZ_SPL_Range_WorstCondition_dB', ...
    'DesignWallSeconds', 'DesignCPUSeconds'});
end

function table_out = local_paired_differences(results, algorithms, ...
        display_names, frequencies, config)
metrics = local_metric_fields();
row_count = numel(config.paired_references) * numel(config.band_names) * ...
    (numel(algorithms) - 1) * numel(metrics);
rows = cell(row_count, 15);
row = 0;
for reference_no = 1:numel(config.paired_references)
    reference = config.paired_references{reference_no};
    if ~ismember(string(reference), algorithms)
        error('FinalCabin:UnknownPairwiseReference', ...
            'Pairwise reference %s is not in the algorithm set.', reference);
    end
    reference_display = local_display_name(string(reference));
    for band_no = 1:numel(config.band_names)
        limits = config.band_limits_hz(band_no, :);
        band = frequencies >= limits(1) & frequencies <= limits(2);
        for algorithm_no = 1:numel(algorithms)
            algorithm = char(algorithms(algorithm_no));
            if strcmp(algorithm, reference)
                continue;
            end
            for metric_no = 1:numel(metrics)
                metric = metrics{metric_no};
                algorithm_values = mean(results.(algorithm).(metric)(:, band), ...
                    2, 'omitnan');
                reference_values = mean(results.(reference).(metric)(:, band), ...
                    2, 'omitnan');
                difference = algorithm_values - reference_values;
                row = row + 1;
                rows(row, :) = {display_names(algorithm_no), ...
                    algorithms(algorithm_no), reference_display, ...
                    string(reference), config.band_names(band_no), ...
                    limits(1), limits(2), string(metric), ...
                    numel(difference), mean(difference), median(difference), ...
                    prctile(difference, 10), prctile(difference, 90), ...
                    min(difference), max(difference)};
            end
        end
    end
end
rows = rows(1:row, :);
table_out = cell2table(rows, 'VariableNames', {'Algorithm', ...
    'AlgorithmField', 'ReferenceAlgorithm', 'ReferenceAlgorithmField', ...
    'Band', 'BandStart_Hz', 'BandEnd_Hz', 'Metric', ...
    'PairedConditionCount', 'MeanDifference', 'MedianDifference', ...
    'P10Difference', 'P90Difference', 'MinimumDifference', ...
    'MaximumDifference'});
end

function diagnostic = local_reference_diagnostic(config, indices, frequencies)
band = frequencies >= config.analysis_band_hz(1) & ...
    frequencies <= config.analysis_band_hz(2);
all_B = cell(1, numel(indices));
all_D = cell(1, numel(indices));
for index_no = 1:numel(indices)
    loaded = load(fullfile(config.data_dir, sprintf('ATF_%d.mat', ...
        indices(index_no))), 'ATF_BZ', 'ATF_DZ');
    all_B{index_no} = loaded.ATF_BZ;
    all_D{index_no} = loaded.ATF_DZ;
end
n = numel(indices);
values = nan(n, 9);
for reference_no = 1:n
    bright = [];
    dark = [];
    for operating_no = 1:n
        if operating_no == reference_no
            continue;
        end
        for frequency_no = find(band)
            bright(end + 1, 1) = norm( ...
                all_B{operating_no}(:, :, frequency_no) - ...
                all_B{reference_no}(:, :, frequency_no), 'fro') / ...
                max(norm(all_B{reference_no}(:, :, frequency_no), ...
                'fro'), realmin); %#ok<AGROW>
            dark(end + 1, 1) = norm( ...
                all_D{operating_no}(:, :, frequency_no) - ...
                all_D{reference_no}(:, :, frequency_no), 'fro') / ...
                max(norm(all_D{reference_no}(:, :, frequency_no), ...
                'fro'), realmin); %#ok<AGROW>
        end
    end
    combined = [bright; dark];
    values(reference_no, :) = [mean(bright), median(bright), ...
        prctile(bright, 90), max(bright), mean(dark), median(dark), ...
        prctile(dark, 90), max(dark), median(combined)];
end
[~, order] = sort(values(:, 9), 'ascend');
rank = nan(n, 1);
rank(order) = 1:n;
percentile = 100 * (rank - 1) / max(n - 1, 1);
ReferenceIndex = indices(:);
SelectedReference = ReferenceIndex == config.reference_index;
diagnostic = table(ReferenceIndex, SelectedReference, values(:, 1), ...
    values(:, 2), values(:, 3), values(:, 4), values(:, 5), ...
    values(:, 6), values(:, 7), values(:, 8), values(:, 9), rank, ...
    percentile, 'VariableNames', {'ReferenceIndex', ...
    'SelectedReference', 'BZ_MeanRelativeDelta', ...
    'BZ_MedianRelativeDelta', 'BZ_P90RelativeDelta', ...
    'BZ_MaxRelativeDelta', 'DZ_MeanRelativeDelta', ...
    'DZ_MedianRelativeDelta', 'DZ_P90RelativeDelta', ...
    'DZ_MaxRelativeDelta', 'CombinedMedianRelativeDelta', ...
    'CentralityRank_LowerIsMoreCentral', ...
    'CentralityPercentile_LowerIsMoreCentral'});
end

function paths = local_create_frequency_figure(summary, algorithms, ...
        display_names, config, output_dir)
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 26, 19]);
layout = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
metrics = {'AC', 'NSRE', 'AE', 'BZ_SPL_Variance', ...
    'BZ_SPL_Std', 'BZ_SPL_Range'};
ylabels = {'AC (dB)', 'Aligned NSRE (dB)', 'AE (dB)', ...
    'BZ SPL variance (dB^2)', 'BZ SPL std. (dB)', ...
    'BZ SPL range (dB)'};
colors = lines(numel(algorithms));
styles = {'-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--', '-.'};
handles = gobjects(numel(algorithms), 1);
for metric_no = 1:numel(metrics)
    ax = nexttile(layout);
    hold(ax, 'on');
    for algorithm_no = 1:numel(algorithms)
        subset = summary(summary.AlgorithmField == ...
            algorithms(algorithm_no), :);
        handles(algorithm_no) = semilogx(ax, subset.Frequency_Hz, ...
            subset.([metrics{metric_no}, '_Mean']), ...
            'LineStyle', styles{algorithm_no}, 'LineWidth', 1.25, ...
            'Color', colors(algorithm_no, :), ...
            'DisplayName', display_names(algorithm_no));
    end
    xlim(ax, config.analysis_band_hz);
    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, ylabels{metric_no});
    local_style_axes(ax, config.figure_font);
end
legend_handle = legend(handles, display_names, ...
    'Orientation', 'horizontal', 'NumColumns', 4);
legend_handle.Layout.Tile = 'south';
paths = local_save_figure(fig, output_dir, ...
    'final_measured_cabin_frequency_metrics');
close(fig);
end

function paths = local_create_band_figure(results, algorithms, ...
        display_names, frequencies, config, output_dir)
band = frequencies >= config.analysis_band_hz(1) & ...
    frequencies <= config.analysis_band_hz(2);
metrics = {'AC', 'NSRE', 'AE', 'BZ_SPL_Std', 'BZ_SPL_Range'};
ylabels = {'AC (dB)', 'Aligned NSRE (dB)', 'AE (dB)', ...
    'BZ SPL std. (dB)', 'BZ SPL range (dB)'};
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 27, 16]);
layout = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
colors = lines(numel(algorithms));
for metric_no = 1:numel(metrics)
    ax = nexttile(layout);
    hold(ax, 'on');
    for algorithm_no = 1:numel(algorithms)
        algorithm = char(algorithms(algorithm_no));
        values = mean(results.(algorithm).(metrics{metric_no})(:, band), ...
            2, 'omitnan');
        center = mean(values);
        lower = prctile(values, 10);
        upper = prctile(values, 90);
        errorbar(ax, algorithm_no, center, center - lower, ...
            upper - center, 'o', 'Color', colors(algorithm_no, :), ...
            'MarkerFaceColor', colors(algorithm_no, :), ...
            'LineWidth', 1, 'CapSize', 4);
    end
    xlim(ax, [0.4, numel(algorithms) + 0.6]);
    xticks(ax, 1:numel(algorithms));
    xticklabels(ax, display_names);
    xtickangle(ax, 45);
    ylabel(ax, ylabels{metric_no});
    title(ax, '100--4000 Hz: mean and P10--P90');
    local_style_axes(ax, config.figure_font);
end
nexttile(layout, 6);
axis off;
text(0.02, 0.92, { ...
    '59 held-out operating conditions', ...
    '100--4000-Hz average per condition', ...
    'Marker: condition mean', ...
    'Error bar: condition-level P10--P90'}, ...
    'FontName', config.figure_font, 'FontSize', 9, ...
    'VerticalAlignment', 'top', 'Interpreter', 'none');
paths = local_save_figure(fig, output_dir, ...
    'final_measured_cabin_band_metrics');
close(fig);
end

function paths = local_create_alpha_figure(selection, output_dir)
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 18, 11]);
ax = axes(fig);
x = 1:height(selection);
yyaxis(ax, 'left');
plot(ax, x, selection.AC_dB, '-o', 'LineWidth', 1.5, ...
    'MarkerSize', 6, 'DisplayName', 'AC');
ylabel(ax, 'Nominal AC (dB)');
yyaxis(ax, 'right');
plot(ax, x, selection.BZ_SPL_Std_dB, '-s', 'LineWidth', 1.5, ...
    'MarkerSize', 6, 'DisplayName', 'SPL std.');
hold(ax, 'on');
plot(ax, x, selection.BZ_SPL_Range_dB, '-d', 'LineWidth', 1.5, ...
    'MarkerSize', 6, 'DisplayName', 'SPL range');
selected_index = find(selection.Selected);
xline(ax, selected_index, ':k', 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');
ylabel(ax, 'BZ spatial nonuniformity (dB)');
xticks(ax, x);
xticklabels(ax, compose('%.0e', selection.Alpha));
xlabel(ax, '\alpha');
legend(ax, 'Location', 'northeast');
local_style_axes(ax, 'Times New Roman');
paths = local_save_figure(fig, output_dir, ...
    'final_macc_nominal_alpha_selection');
close(fig);
end

function local_style_axes(ax, font)
set(ax, 'FontName', font, 'FontSize', 9, 'LineWidth', 0.8, ...
    'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.2);
ax.Toolbar.Visible = 'off';
end

function paths = local_save_figure(fig, output_dir, basename)
paths = {fullfile(output_dir, [basename, '.png']), ...
    fullfile(output_dir, [basename, '.pdf']), ...
    fullfile(output_dir, [basename, '.fig'])};
exportgraphics(fig, paths{1}, 'Resolution', 300);
exportgraphics(fig, paths{2}, 'ContentType', 'vector');
savefig(fig, paths{3});
end

function name = local_display_name(field)
switch char(field)
    case 'ACC_Reg', name = "ACC-Reg";
    case 'ACC_PM', name = "ACC-PM";
    case 'wcRACC', name = "WCRACC";
    case 'NoCT_WCRACC', name = "NoCT-WCRACC";
    case 'Full_WCRACC', name = "Full-WCRACC";
    case 'POTDC_RACC', name = "POTDC-RACC";
    case 'RACC_PM_Subpro', name = "RACC-PM";
    case 'M_ACC', name = "M-ACC";
    otherwise, name = strrep(string(field), '_', '-');
end
end

function local_ensure_directory(directory)
if ~exist(directory, 'dir')
    mkdir(directory);
end
end

function local_safe_fclose(file_id)
if isnumeric(file_id) && isscalar(file_id) && file_id >= 0
    fclose(file_id);
end
end
