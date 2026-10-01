function outputs = run_planarity_perceptual_experiment(mode, ...
    evaluation_temperature_celsius, filter_file)
%RUN_PLANARITY_PERCEPTUAL_EXPERIMENT Reviewer-oriented planarity experiment.
%   MODE = 'frequency' generates the new RIRs and frequency-domain results.
%   MODE = 'quick' additionally evaluates S01 and M06.
%   MODE = 'full' evaluates all ten speech and ten music excerpts.
%   EVALUATION_TEMPERATURE_CELSIUS selects the evaluation RIR/ATF.
%   FILTER_FILE optionally selects a completed step-2 result. By default the
%   ten-algorithm rho=10, nu=0.01 production result is used. A filter file
%   may additionally provide PERCEPTUAL_CONFIG and CALIBRATION_ATF_BZ_CTRL;
%   these metadata support temperature-informed SICER-VAST snapshots without
%   changing the default fixed-filter protocol.

if nargin < 1
    mode = 'quick';
end
if nargin < 2
    evaluation_temperature_celsius = 22.5;
end
mode = lower(string(mode));
if ~ismember(mode, ["frequency", "quick", "full"])
    error('mode must be ''frequency'', ''quick'', or ''full''.');
end
if ~isscalar(evaluation_temperature_celsius) || ...
        ~isfinite(evaluation_temperature_celsius)
    error('evaluation_temperature_celsius must be a finite scalar.');
end

config = local_configuration(mode, evaluation_temperature_celsius);
[project_root, response_dir] = locate_project();
if nargin < 3 || strlength(string(filter_file)) == 0
    filter_file = fullfile(response_dir, ...
        'racc_pm_parameter_comparison_results', ...
        'rho10_nu0p01_all_run_20260827_002614', ...
        'performance_matrices_20260827_021428.mat');
else
    filter_file = char(string(filter_file));
end
if ~isfile(filter_file)
    error('Filter result file does not exist: %s', filter_file);
end
config.filter_file = filter_file;
filter_data = load_filter_data(filter_file);
config = apply_filter_configuration(config, filter_data);
addpath(genpath(fullfile(project_root, 'src')));
addpath(genpath(fullfile(project_root, 'lib')));
addpath(genpath(fullfile(response_dir, 'PEAQ')));

data_dir = fullfile(response_dir, 'planarity_perceptual_data');
run_started = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
temperature_tag = make_temperature_tag(config.evaluation_temperature_celsius);
run_id = sprintf('%s_%s', ...
    char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')), temperature_tag);
result_dir = fullfile(response_dir, 'planarity_perceptual_results', ...
    char(mode), ['run_', run_id]);
figure_dir = fullfile(result_dir, 'figures');
temporary_dir = fullfile(result_dir, 'temporary_wav');
ensure_directory(data_dir);
ensure_directory(result_dir);
ensure_directory(figure_dir);
ensure_directory(temporary_dir);

diary_file = fullfile(result_dir, 'run_log.txt');
diary('off');
diary(diary_file);
diary_cleanup = onCleanup(@() diary('off'));
fprintf('Planarity/perceptual experiment started: %s\n', ...
    string(run_started));
fprintf('Run ID: %s\n', run_id);
fprintf('Mode: %s\n', mode);
fprintf('Design/calibration temperature: %.2f C\n', ...
    config.calibration_temperature_celsius);
fprintf('Evaluation temperature: %.2f C\n', ...
    config.evaluation_temperature_celsius);
fprintf('Filter source: %s\n', config.filter_file);
fprintf('Primary setting: rho=%.6g, nu=%.6g\n', ...
    config.primary_rho, config.primary_nu);
fprintf('Calibration protocol: %s\n', config.calibration_protocol);
if strlength(config.method_note) > 0
    fprintf('Filter method: %s\n', config.method_note);
end

geometry_file = fullfile(project_root, 'data', 'arrayGeometry', ...
    'array_layout.mat');
nominal_file = fullfile(project_root, 'data', 'SimulateRIR', ...
    'temperature', 'Data_T-22.50.mat');
manifest_file = fullfile(response_dir, 'perceptual_audio', 'manifest.csv');

geometry_data = load(geometry_file, 'roomArray');
nominal_data = load(nominal_file, 'ATF_BZ');
freq_params = filter_data.para.freq_params;
validate_frequency_parameters(freq_params);
validate_filter_data(filter_data, config);

calibration_atf = nominal_data.ATF_BZ.ctrl;
if isfield(filter_data, 'calibration_ATF_BZ_ctrl')
    calibration_atf = filter_data.calibration_ATF_BZ_ctrl;
end

point_definition = define_evaluation_points();
rir_file = fullfile(data_dir, sprintf('ear_pair_rir_%s.mat', ...
    temperature_tag));
evaluation_data = load_or_generate_rirs(rir_file, point_definition, ...
    geometry_data.roomArray, freq_params, config);

reference = build_plane_wave_reference(point_definition, ...
    geometry_data.roomArray, freq_params, config);
[calibrated_weights, calibration_table] = calibrate_filters( ...
    filter_data, calibration_atf, geometry_data.roomArray.BZ_ctrl, ...
    reference, config);

writetable(point_definition.table, fullfile(data_dir, ...
    'ear_pair_coordinates.csv'));
writetable(calibration_table, fullfile(result_dir, ...
    'frequency_calibration.csv'));
calibration_tag = make_temperature_tag( ...
    config.calibration_temperature_celsius);
save(fullfile(result_dir, sprintf('calibrated_filters_%s.mat', ...
    calibration_tag)), ...
    'calibrated_weights', 'calibration_table', 'reference', ...
    'point_definition', 'freq_params', 'config', 'filter_file', '-v7.3');

[frequency_table, frequency_summary] = evaluate_frequency_domain( ...
    calibrated_weights, evaluation_data.ATF, filter_data, reference, ...
    freq_params, config);
writetable(frequency_table, fullfile(result_dir, ...
    'frequency_metrics.csv'));
writetable(frequency_summary, fullfile(result_dir, ...
    'frequency_summary.csv'));

create_geometry_figure(geometry_data.roomArray, point_definition, ...
    figure_dir, config);
create_frequency_figures(frequency_table, figure_dir, config);

outputs = struct();
outputs.run_id = run_id;
outputs.result_dir = result_dir;
outputs.design_temperature_celsius = config.design_temperature_celsius;
outputs.evaluation_temperature_celsius = ...
    config.evaluation_temperature_celsius;
outputs.rir_file = rir_file;
outputs.filter_file = filter_file;
outputs.frequency_metrics = frequency_table;
outputs.frequency_summary = frequency_summary;
outputs.calibration = calibration_table;

if mode ~= "frequency"
    [renderers, reconstruction_table] = build_broadband_renderers( ...
        calibrated_weights, evaluation_data.IR, evaluation_data.ATF, ...
        freq_params, reference);
    writetable(reconstruction_table, fullfile(result_dir, ...
        'broadband_realization_diagnostics.csv'));
    save(fullfile(result_dir, 'broadband_renderers.mat'), ...
        'renderers', 'reconstruction_table', 'reference', ...
        'point_definition', 'freq_params', '-v7.3');

    manifest = readtable(manifest_file, 'TextType', 'string');
    manifest = select_manifest_rows(manifest, mode);
    perceptual = evaluate_audio_set(manifest, response_dir, result_dir, ...
        temporary_dir, renderers, reference, config);
    outputs.perceptual = perceptual;
    outputs.reconstruction = reconstruction_table;
end

remove_temporary_wavs(temporary_dir);
fprintf('\nFrequency-domain summary:\n');
disp(frequency_summary);
fprintf('Experiment completed: %s\n', ...
    string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf('Output: %s\n', result_dir);
write_completion_marker(result_dir, mode, run_id, config);
clear diary_cleanup;
end

function config = local_configuration(mode, evaluation_temperature_celsius)
config = struct();
config.mode = mode;
config.design_temperature_celsius = 22.5;
config.calibration_temperature_celsius = 22.5;
config.evaluation_temperature_celsius = evaluation_temperature_celsius;
config.target_spl_db = 76;
config.reference_pressure_pa = 20e-6;
config.target_pressure_pa = config.reference_pressure_pa * ...
    10^(config.target_spl_db / 20);
config.desired_sound_speed_mps = 343;
config.virtual_source_index = 13;
config.causal_margin_s = 0.010;
config.aliasing_frequency_hz = 1630;
config.summary_upper_frequency_hz = 4000;
config.rir_fade_time_ms = 8;
config.force_regenerate_rir = false;
config.primary_rho = 10;
config.primary_nu = 0.01;
config.method_note = "";
config.calibration_protocol = ...
    "Target-referenced methods retain their native unit-target solution " + ...
    "and receive only the common 76-dB physical scale. Contrast-only " + ...
    "methods use nominal control ATFs to resolve arbitrary narrowband " + ...
    "level and global phase. Evaluation pairs are never used for scaling.";
config.python_executable = getenv('PESQ_PYTHON');
config.listening_item_ids = ["S01", "M06"];
config.algorithm_fields = ["ACC", "ACC_Reg", "PM", "ACC_PM", ...
    "wcRACC", "NoCT_WCRACC", "Full_WCRACC", "POTDC_RACC", ...
    "RPM", "RACC_PM_Subpro"];
config.target_referenced_algorithms = ["PM", "ACC_PM", "RPM", ...
    "RACC_PM_Subpro"];
config.listening_algorithms = config.algorithm_fields;
config.display_algorithms = config.algorithm_fields;
config.display_names = ["ACC", "ACC-Reg", "PM", "ACC-PM", ...
    "WCRACC", "NoCT-WCRACC", "Full-WCRACC", "POTDC-RACC", ...
    "RPM", "RACC-PM"];
config.figure_font = 'Times New Roman';
config.figure_size_cm = [18, 13];
end

function filter_data = load_filter_data(filter_file)
required_names = {'all_filters', 'algorithms_to_test', 'para', ...
    'val_results'};
optional_names = {'perceptual_config', 'calibration_ATF_BZ_ctrl', ...
    'calibration_temperature_celsius', 'filter_temperature_celsius', ...
    'source_result_file'};
contents = whos('-file', filter_file);
available_names = string({contents.name});
missing = string(required_names(~ismember(string(required_names), ...
    available_names)));
if ~isempty(missing)
    error('Missing filter-result variable(s): %s.', strjoin(missing, ', '));
end
present_optional = optional_names(ismember(string(optional_names), ...
    available_names));
filter_data = load(filter_file, required_names{:}, present_optional{:});
end

function config = apply_filter_configuration(config, filter_data)
if ~isfield(filter_data, 'perceptual_config')
    return;
end
override = filter_data.perceptual_config;
required = {'algorithm_fields', 'display_algorithms', ...
    'listening_algorithms', 'display_names'};
for field_no = 1:numel(required)
    if ~isfield(override, required{field_no})
        error('PERCEPTUAL_CONFIG is missing %s.', required{field_no});
    end
end
config.algorithm_fields = string(override.algorithm_fields(:)).';
config.display_algorithms = string(override.display_algorithms(:)).';
config.listening_algorithms = string(override.listening_algorithms(:)).';
config.display_names = string(override.display_names(:)).';
if numel(config.display_algorithms) ~= numel(config.display_names)
    error('PERCEPTUAL_CONFIG display algorithms and names must align.');
end
if isfield(override, 'design_temperature_celsius')
    config.design_temperature_celsius = ...
        double(override.design_temperature_celsius);
end
if isfield(override, 'calibration_temperature_celsius')
    config.calibration_temperature_celsius = ...
        double(override.calibration_temperature_celsius);
elseif isfield(filter_data, 'calibration_temperature_celsius')
    config.calibration_temperature_celsius = ...
        double(filter_data.calibration_temperature_celsius);
end
if isfield(override, 'method_note')
    config.method_note = string(override.method_note);
end
if isfield(override, 'target_referenced_algorithms')
    config.target_referenced_algorithms = ...
        string(override.target_referenced_algorithms(:)).';
else
    config.target_referenced_algorithms = ...
        infer_target_referenced_algorithms(config.algorithm_fields);
end
config.calibration_protocol = sprintf([ ...
    'Target-referenced methods retain their native unit-target solution ', ...
    'and receive only the common 76-dB physical scale. Contrast-only ', ...
    'methods use %.2f-C control ATFs to resolve arbitrary narrowband ', ...
    'level/global phase. Evaluation pairs are not used for scaling.'], ...
    config.calibration_temperature_celsius);
end

function algorithms = infer_target_referenced_algorithms(candidates)
known = ["PM", "ACC_PM", "WCRPM", "RPM", "RACC_PM", ...
    "RACC_PM_Sub", "RACC_PM_Subpro", "RACC_PM_GLS", ...
    "SICER_VAST_VL", "GT_VAST_VL"];
algorithms = candidates(ismember(candidates, known));
end

function [project_root, response_dir] = locate_project()
response_dir = fileparts(mfilename('fullpath'));
if ~strcmpi(get_last_folder(response_dir), 'RQ_response')
    error('PerceptualExperiment:UnexpectedLocation', ...
        'Place this script in the repository RQ_response directory.');
end
project_root = fileparts(response_dir);
end

function validate_frequency_parameters(freq_params)
required = {'fs', 'rir_len', 'nfft', 'target_freqs', ...
    'f_target_indices', 'fade_samples', 'w_tail'};
for field_no = 1:numel(required)
    if ~isfield(freq_params, required{field_no})
        error('Missing frequency parameter: %s.', required{field_no});
    end
end
if freq_params.nfft ~= freq_params.rir_len
    error('This experiment expects nfft and rir_len to be identical.');
end
end

function validate_filter_data(filter_data, config)
required_top_level = {'all_filters', 'algorithms_to_test', 'para', ...
    'val_results'};
for field_no = 1:numel(required_top_level)
    if ~isfield(filter_data, required_top_level{field_no})
        error('Missing filter-result variable: %s.', ...
            required_top_level{field_no});
    end
end

actual_algorithms = string(filter_data.algorithms_to_test(:)).';
expected_algorithms = config.algorithm_fields;
if ~isequal(actual_algorithms, expected_algorithms)
    error(['The selected filter result algorithm order does not match ', ...
        'the active perceptual configuration.']);
end
if any(~ismember(config.target_referenced_algorithms, expected_algorithms))
    error(['A target-referenced algorithm is absent from the selected ', ...
        'filter result.']);
end

frequency_count = numel(filter_data.para.freq_params.target_freqs);
for algorithm_no = 1:numel(expected_algorithms)
    algorithm = char(expected_algorithms(algorithm_no));
    if ~isfield(filter_data.all_filters, algorithm) || ...
            ~isfield(filter_data.val_results, algorithm)
        error('Missing filter or evaluation result for %s.', algorithm);
    end
    weights = filter_data.all_filters.(algorithm);
    if size(weights, 2) ~= frequency_count || ...
            any(~isfinite(real(weights(:)))) || ...
            any(~isfinite(imag(weights(:))))
        error('Invalid filter matrix for %s.', algorithm);
    end
    planarity = filter_data.val_results.(algorithm).Planarity;
    if size(planarity, 2) ~= frequency_count
        error('Invalid planarity matrix for %s.', algorithm);
    end
end
end

function point_definition = define_evaluation_points()
positions = [ ...
    2.51, 2.00, 1.60; ...
    2.69, 2.00, 1.60; ...
    2.60, 1.91, 1.60; ...
    2.60, 2.09, 1.60; ...
    1.31, 2.00, 1.60; ...
    1.49, 2.00, 1.60; ...
    1.40, 1.91, 1.60; ...
    1.40, 2.09, 1.60];
labels = ["BZ_horizontal_left"; "BZ_horizontal_right"; ...
    "BZ_vertical_lower"; "BZ_vertical_upper"; ...
    "DZ_horizontal_left"; "DZ_horizontal_right"; ...
    "DZ_vertical_lower"; "DZ_vertical_upper"];
zones = [repmat("BZ", 4, 1); repmat("DZ", 4, 1)];
orientations = repmat(["horizontal"; "horizontal"; ...
    "vertical"; "vertical"], 2, 1);
ear_number = repmat([1; 2; 1; 2], 2, 1);

point_definition = struct();
point_definition.version = 2;
point_definition.positions = positions;
point_definition.labels = labels;
point_definition.zones = zones;
point_definition.orientations = orientations;
point_definition.ear_number = ear_number;
point_definition.bz_indices = 1:4;
point_definition.dz_indices = 5:8;
point_definition.horizontal_local_indices = [1, 2];
point_definition.vertical_local_indices = [3, 4];
point_definition.table = table((1:8)', labels, zones, orientations, ...
    ear_number, positions(:, 1), positions(:, 2), positions(:, 3), ...
    'VariableNames', {'PointIndex', 'PointLabel', 'Zone', 'Orientation', ...
    'EarNumber', 'X_m', 'Y_m', 'Z_m'});
end

function evaluation_data = load_or_generate_rirs(output_file, ...
    point_definition, room_array, freq_params, config)

if exist(output_file, 'file') && ~config.force_regenerate_rir
    loaded = load(output_file, 'evaluation_data');
    evaluation_data = loaded.evaluation_data;
    same_positions = isequal(size(evaluation_data.positions), ...
        size(point_definition.positions)) && max(abs( ...
        evaluation_data.positions(:) - point_definition.positions(:))) < 1e-12;
    same_parameters = evaluation_data.fs == freq_params.fs && ...
        evaluation_data.rir_len == freq_params.rir_len && ...
        evaluation_data.temperature_celsius == ...
        config.evaluation_temperature_celsius;
    if same_positions && same_parameters
        fprintf('Using cached exact-point RIRs: %s\n', output_file);
        return;
    end
end

positions = point_definition.positions;
sources = room_array.s;
number_of_points = size(positions, 1);
number_of_sources = size(sources, 1);
rir_length = freq_params.rir_len;
fs = freq_params.fs;
sound_speed = temp2speed(config.evaluation_temperature_celsius);
fade_samples = round(config.rir_fade_time_ms / 1000 * fs);
if fade_samples ~= freq_params.fade_samples
    error('Configured fade length does not match the original experiment.');
end
w_tail = hann(2 * fade_samples);
w_tail = w_tail(fade_samples + 1:end);
if norm(w_tail(:) - freq_params.w_tail(:)) > 1e-12
    error('Recreated Hann fade differs from the stored experiment window.');
end

fprintf('Generating %d exact-point RIRs for each of %d loudspeakers.\n', ...
    number_of_points, number_of_sources);
IR = zeros(number_of_points, number_of_sources, rir_length, 'single');
parfor source_no = 1:number_of_sources
    current = rir_generator(sound_speed, fs, positions, ...
        sources(source_no, :), room_array.roomSize, ...
        freq_params.beta, rir_length);
    current(:, end - fade_samples + 1:end) = ...
        current(:, end - fade_samples + 1:end) .* w_tail.';
    IR(:, source_no, :) = reshape(single(current), ...
        number_of_points, 1, rir_length);
end
ATF = compute_atf(IR, freq_params);

evaluation_data = struct();
evaluation_data.version = 2;
evaluation_data.description = ['RIRs at exact 18-cm horizontal and vertical ' ...
    'ear pairs in both zones.'];
evaluation_data.positions = positions;
evaluation_data.labels = point_definition.labels;
evaluation_data.temperature_celsius = config.evaluation_temperature_celsius;
evaluation_data.sound_speed_mps = sound_speed;
evaluation_data.fs = fs;
evaluation_data.rir_len = rir_length;
evaluation_data.fade_time_ms = config.rir_fade_time_ms;
evaluation_data.fade_samples = fade_samples;
evaluation_data.w_tail = w_tail;
evaluation_data.beta = freq_params.beta;
evaluation_data.room_size = room_array.roomSize;
evaluation_data.source_positions = sources;
evaluation_data.IR = IR;
evaluation_data.ATF = ATF;
save(output_file, 'evaluation_data', '-v7.3');
fprintf('Saved exact-point RIR/ATF data: %s\n', output_file);
end

function reference = build_plane_wave_reference(point_definition, ...
    room_array, freq_params, config)

bright_center = [2.6, 2.0, 1.6];
source_position = room_array.s(config.virtual_source_index, :);
direction = bright_center - source_position;
direction = direction ./ norm(direction);

reference = struct();
reference.fs = freq_params.fs;
reference.direction = direction;
reference.sound_speed_mps = config.desired_sound_speed_mps;
reference.target_pressure_pa = config.target_pressure_pa;
reference.positions = point_definition.positions;
reference.target_atf = calculate_plane_wave_atf( ...
    point_definition.positions, freq_params.target_freqs, direction, ...
    config.desired_sound_speed_mps);
reference.calibration_atf = calculate_plane_wave_atf( ...
    room_array.BZ_eval, freq_params.target_freqs, direction, ...
    config.desired_sound_speed_mps);

geometric_delays = point_definition.positions * direction.' .* ...
    freq_params.fs / config.desired_sound_speed_mps;
causal_margin = round(config.causal_margin_s * freq_params.fs);
common_delay = ceil(max(0, -min(geometric_delays))) + causal_margin;
reference.geometric_delay_samples = geometric_delays;
reference.common_delay_samples = common_delay;
reference.total_delay_samples = common_delay + geometric_delays;
reference.point_labels = point_definition.labels;
end

function atf = calculate_plane_wave_atf(positions, frequencies, ...
    direction, sound_speed)
wave_vectors = (2 * pi .* frequencies(:) / sound_speed) .* direction;
atf = exp(-1i .* (positions * wave_vectors.'));
end

function [calibrated_weights, calibration_table] = calibrate_filters( ...
    filter_data, calibration_atf, calibration_positions, reference, config)

algorithms = string(filter_data.algorithms_to_test(:));
frequencies = filter_data.para.freq_params.target_freqs(:);
desired = calculate_plane_wave_atf(calibration_positions, frequencies, ...
    reference.direction, reference.sound_speed_mps);
calibrated_weights = struct();
rows = cell(numel(algorithms), 9);

for algorithm_no = 1:numel(algorithms)
    algorithm = algorithms(algorithm_no);
    raw = double(filter_data.all_filters.(char(algorithm)));
    weights = zeros(size(raw));
    scales = zeros(numel(frequencies), 1);
    rotations = zeros(numel(frequencies), 1);
    raw_spl = zeros(numel(frequencies), 1);
    calibrated_spl = zeros(numel(frequencies), 1);
    is_target_referenced = ismember(algorithm, ...
        config.target_referenced_algorithms);
    if is_target_referenced
        preparation_mode = "Native target solution + common physical scale";
    else
        preparation_mode = "Contrast-only control-point level/phase anchor";
    end
    for frequency_no = 1:numel(frequencies)
        pressure = calibration_atf(:, :, frequency_no) * raw(:, frequency_no);
        pressure_rms = sqrt(mean(abs(pressure).^2));
        if is_target_referenced
            % The optimization already used the unit-magnitude target field.
            % Preserve its attained amplitude and phase; only map the
            % normalized solution to the stated 76-dB physical pressure.
            scales(frequency_no) = config.target_pressure_pa;
            rotations(frequency_no) = 0;
        else
            % A contrast-only eigenvector has arbitrary complex scale. Fix
            % that indeterminacy from nominal control points before forming
            % a broadband renderer. No evaluation point is used here.
            scales(frequency_no) = config.target_pressure_pa / ...
                max(pressure_rms, eps);
            phase_anchor = sum(conj(pressure) .* desired(:, frequency_no));
            rotations(frequency_no) = angle(phase_anchor);
        end
        weights(:, frequency_no) = raw(:, frequency_no) .* ...
            scales(frequency_no) .* exp(1i * rotations(frequency_no));
        raw_spl(frequency_no) = pressure_to_spl(pressure, ...
            config.reference_pressure_pa);
        calibrated_pressure = calibration_atf(:, :, frequency_no) * ...
            weights(:, frequency_no);
        calibrated_spl(frequency_no) = pressure_to_spl( ...
            calibrated_pressure, config.reference_pressure_pa);
    end
    calibrated_weights.(char(algorithm)) = weights;
    rows(algorithm_no, :) = {char(algorithm), char(preparation_mode), ...
        config.target_spl_db, ...
        median(scales), min(scales), max(scales), ...
        median(abs(rad2deg(rotations))), mean(raw_spl), ...
        max(abs(calibrated_spl - config.target_spl_db))};
end

calibration_table = cell2table(rows, 'VariableNames', ...
    {'Algorithm', 'PreparationMode', 'TargetSPL_dB', ...
    'MedianAmplitudeScale', ...
    'MinimumAmplitudeScale', 'MaximumAmplitudeScale', ...
    'MedianAbsPhaseRotation_deg', 'MeanRawBZSPL_dB', ...
    'MaximumPreparedLevelDeviationFromTarget_dB'});
end

function [frequency_table, summary] = evaluate_frequency_domain( ...
    calibrated_weights, evaluation_atf, filter_data, reference, ...
    freq_params, config)

algorithms = string(filter_data.algorithms_to_test(:));
frequencies = freq_params.target_freqs(:);
temperature = filter_data.para.temperature_vector_celsius(:);
[temperature_error, evaluation_index] = min(abs(temperature - ...
    config.evaluation_temperature_celsius));
if temperature_error > 1e-9
    error('Evaluation temperature %.3f C is absent from stored results.', ...
        config.evaluation_temperature_celsius);
end
rows = cell(numel(algorithms) * numel(frequencies), 13);
horizontal_bz_spl = zeros(numel(algorithms), numel(frequencies));
vertical_bz_spl = zeros(numel(algorithms), numel(frequencies));
horizontal_dz_spl = zeros(numel(algorithms), numel(frequencies));
vertical_dz_spl = zeros(numel(algorithms), numel(frequencies));
row_no = 0;

for algorithm_no = 1:numel(algorithms)
    algorithm = algorithms(algorithm_no);
    weights = calibrated_weights.(char(algorithm));
    planarity = filter_data.val_results.(char(algorithm)).Planarity( ...
        evaluation_index, :);
    for frequency_no = 1:numel(frequencies)
        row_no = row_no + 1;
        pressure = evaluation_atf(:, :, frequency_no) * ...
            weights(:, frequency_no);
        desired = config.target_pressure_pa .* ...
            reference.target_atf(:, frequency_no);
        h_bz = pair_metrics(pressure(1:2), desired(1:2), config);
        v_bz = pair_metrics(pressure(3:4), desired(3:4), config);
        h_dz = pair_metrics(pressure(5:6), [], config);
        v_dz = pair_metrics(pressure(7:8), [], config);
        if ismember(algorithm, config.target_referenced_algorithms)
            standard_nsre_h = h_bz.nsre_db;
            standard_nsre_v = v_bz.nsre_db;
            contrast_error_h = NaN;
            contrast_error_v = NaN;
        else
            standard_nsre_h = NaN;
            standard_nsre_v = NaN;
            contrast_error_h = h_bz.nsre_db;
            contrast_error_v = v_bz.nsre_db;
        end
        horizontal_bz_spl(algorithm_no, frequency_no) = h_bz.spl_db;
        vertical_bz_spl(algorithm_no, frequency_no) = v_bz.spl_db;
        horizontal_dz_spl(algorithm_no, frequency_no) = h_dz.spl_db;
        vertical_dz_spl(algorithm_no, frequency_no) = v_dz.spl_db;
        rows(row_no, :) = {char(algorithm), frequencies(frequency_no), ...
            planarity(frequency_no), h_bz.spl_db - h_dz.spl_db, ...
            v_bz.spl_db - v_dz.spl_db, h_bz.level_imbalance_db, ...
            v_bz.level_imbalance_db, ...
            h_bz.relative_phase_error_deg, ...
            v_bz.relative_phase_error_deg, ...
            standard_nsre_h, standard_nsre_v, ...
            contrast_error_h, contrast_error_v};
    end
end

frequency_table = cell2table(rows, 'VariableNames', ...
    {'Algorithm', 'Frequency_Hz', 'Planarity_percent', ...
    'Horizontal_Attenuation_dB', 'Vertical_Attenuation_dB', ...
    'BZ_Horizontal_LevelImbalance_dB', ...
    'BZ_Vertical_LevelImbalance_dB', ...
    'BZ_Horizontal_PhaseError_deg', 'BZ_Vertical_PhaseError_deg', ...
    'StandardNSRE_Horizontal_dB', 'StandardNSRE_Vertical_dB', ...
    'ContrastCalibratedError_Horizontal_dB', ...
    'ContrastCalibratedError_Vertical_dB'});

summary_rows = cell(numel(algorithms), 13);
analysis_band = frequencies <= config.summary_upper_frequency_hz;
aliasing_band = frequencies <= config.aliasing_frequency_hz;
for algorithm_no = 1:numel(algorithms)
    algorithm = algorithms(algorithm_no);
    subset = frequency_table(string(frequency_table.Algorithm) == algorithm, :);
    summary_rows(algorithm_no, :) = {char(algorithm), ...
        mean(subset.Planarity_percent(analysis_band), 'omitnan'), ...
        mean(subset.Planarity_percent(aliasing_band), 'omitnan'), ...
        mean(subset.BZ_Horizontal_LevelImbalance_dB(analysis_band), 'omitnan'), ...
        mean(subset.BZ_Vertical_LevelImbalance_dB(analysis_band), 'omitnan'), ...
        mean(subset.BZ_Horizontal_PhaseError_deg(analysis_band), 'omitnan'), ...
        mean(subset.BZ_Vertical_PhaseError_deg(analysis_band), 'omitnan'), ...
        mean(subset.StandardNSRE_Horizontal_dB(analysis_band), 'omitnan'), ...
        mean(subset.StandardNSRE_Vertical_dB(analysis_band), 'omitnan'), ...
        mean(subset.ContrastCalibratedError_Horizontal_dB(analysis_band), ...
        'omitnan'), ...
        mean(subset.ContrastCalibratedError_Vertical_dB(analysis_band), ...
        'omitnan'), ...
        aggregate_attenuation(horizontal_bz_spl(algorithm_no, analysis_band), ...
        horizontal_dz_spl(algorithm_no, analysis_band)), ...
        aggregate_attenuation(vertical_bz_spl(algorithm_no, analysis_band), ...
        vertical_dz_spl(algorithm_no, analysis_band))};
end

summary = cell2table(summary_rows, 'VariableNames', ...
    {'Algorithm', 'PlanarityMean_100_4000Hz_percent', ...
    'PlanarityMean_BelowAliasing_percent', ...
    'BZ_Horizontal_MeanLevelImbalance_dB', ...
    'BZ_Vertical_MeanLevelImbalance_dB', ...
    'BZ_Horizontal_MeanPhaseError_deg', ...
    'BZ_Vertical_MeanPhaseError_deg', ...
    'StandardNSRE_Horizontal_Mean_dB', ...
    'StandardNSRE_Vertical_Mean_dB', ...
    'ContrastCalibratedError_Horizontal_Mean_dB', ...
    'ContrastCalibratedError_Vertical_Mean_dB', ...
    'Horizontal_AggregateAttenuation_dB', ...
    'Vertical_AggregateAttenuation_dB'});
end

function metrics = pair_metrics(pressure, desired, config)
metrics = struct();
metrics.spl_db = pressure_to_spl(pressure, config.reference_pressure_pa);
metrics.level_imbalance_db = abs(20 * log10( ...
    (abs(pressure(1)) + eps) / (abs(pressure(2)) + eps)));
if isempty(desired)
    metrics.relative_phase_error_deg = NaN;
    metrics.nsre_db = NaN;
else
    actual_relative = pressure(2) * conj(pressure(1));
    desired_relative = desired(2) * conj(desired(1));
    metrics.relative_phase_error_deg = abs(rad2deg(angle( ...
        actual_relative * conj(desired_relative))));
    metrics.nsre_db = 10 * log10(max(norm(pressure - desired)^2 / ...
        max(norm(desired)^2, realmin), realmin));
end
end

function spl_db = pressure_to_spl(pressure, reference_pressure)
mean_square = mean(abs(pressure(:)).^2);
spl_db = 10 * log10(max(mean_square, realmin) / reference_pressure^2);
end

function value = aggregate_attenuation(bright_spl_db, dark_spl_db)
bright_energy = sum(10.^(bright_spl_db / 10), 'omitnan');
dark_energy = sum(10.^(dark_spl_db / 10), 'omitnan');
value = 10 * log10(bright_energy / max(dark_energy, realmin));
end

function [renderers, reconstruction_table] = build_broadband_renderers( ...
    calibrated_weights, selected_rirs, selected_atfs, freq_params, ...
    reference)

algorithms = string(fieldnames(calibrated_weights));
number_of_points = size(selected_rirs, 1);
number_of_sources = size(selected_rirs, 2);
rir_length = size(selected_rirs, 3);
filter_length = freq_params.nfft;
effective_length = filter_length + rir_length - 1;
convolution_length = 2^nextpow2(effective_length);
full_frequencies = (0:(filter_length / 2)) .* ...
    freq_params.fs / filter_length;
target_frequencies = freq_params.target_freqs(:);
target_indices = freq_params.f_target_indices(:);
phase_shift = exp(-1i * 2 * pi .* full_frequencies .* ...
    reference.common_delay_samples / freq_params.fs);

rir_fft = cell(number_of_points, 1);
for point_no = 1:number_of_points
    point_rirs = reshape(double(selected_rirs(point_no, :, :)), ...
        number_of_sources, rir_length);
    rir_fft{point_no} = fft(point_rirs, convolution_length, 2);
end

renderers = struct();
rows = cell(numel(algorithms), 6);
valid_target = target_indices < (filter_length / 2 + 1);
for algorithm_no = 1:numel(algorithms)
    algorithm = algorithms(algorithm_no);
    weights = calibrated_weights.(char(algorithm));
    one_sided = zeros(number_of_sources, filter_length / 2 + 1);
    for source_no = 1:number_of_sources
        one_sided(source_no, :) = interp1([0; target_frequencies], ...
            [0; weights(source_no, :).'], full_frequencies, 'linear', 0);
    end
    one_sided = one_sided .* phase_shift;
    one_sided(:, 1) = real(one_sided(:, 1));
    one_sided(:, end) = real(one_sided(:, end));
    full_spectrum = [one_sided, conj(one_sided(:, end - 1:-1:2))];
    time_filters = real(ifft(full_spectrum, [], 2));

    reconstructed = fft(time_filters, filter_length, 2);
    expected_filters = weights(:, valid_target) .* ...
        phase_shift(target_indices(valid_target));
    filter_error = norm(reconstructed(:, target_indices(valid_target)) - ...
        expected_filters, 'fro') / max(norm(expected_filters, 'fro'), eps);

    filter_fft = fft(time_filters, convolution_length, 2);
    effective_irs = zeros(number_of_points, effective_length);
    for point_no = 1:number_of_points
        spectrum = sum(rir_fft{point_no} .* filter_fft, 1);
        current = real(ifft(spectrum));
        effective_irs(point_no, :) = current(1:effective_length);
    end

    valid_frequency_indices = find(valid_target);
    expected_point_response = zeros(number_of_points, ...
        numel(valid_frequency_indices));
    for column_no = 1:numel(valid_frequency_indices)
        frequency_no = valid_frequency_indices(column_no);
        expected_point_response(:, column_no) = ...
            selected_atfs(:, :, frequency_no) * weights(:, frequency_no) .* ...
            exp(-1i * 2 * pi * target_frequencies(frequency_no) * ...
            reference.common_delay_samples / freq_params.fs);
    end
    realized_point_response = zeros(size(expected_point_response));
    valid_frequencies = target_frequencies(valid_target);
    for point_no = 1:number_of_points
        realized_point_response(point_no, :) = freqz( ...
            effective_irs(point_no, :), 1, valid_frequencies, ...
            freq_params.fs).';
    end
    point_error = norm(realized_point_response - expected_point_response, ...
        'fro') / max(norm(expected_point_response, 'fro'), eps);

    [peak_times_ms, late_energy_db] = effective_ir_diagnostics( ...
        effective_irs, freq_params.fs);
    renderers.(char(algorithm)).time_filters = time_filters;
    renderers.(char(algorithm)).effective_irs = effective_irs;
    rows(algorithm_no, :) = {char(algorithm), filter_error, point_error, ...
        median(peak_times_ms), max(peak_times_ms), median(late_energy_db)};
end

reconstruction_table = cell2table(rows, 'VariableNames', ...
    {'Algorithm', 'FilterTargetBinRelativeError', ...
    'PointPressureRelativeError', 'MedianEffectiveIRPeakTime_ms', ...
    'MaximumEffectiveIRPeakTime_ms', 'MedianLateEnergyAfter20ms_dB'});
end

function [peak_times_ms, late_energy_db] = effective_ir_diagnostics(irs, fs)
number_of_points = size(irs, 1);
peak_times_ms = zeros(number_of_points, 1);
late_energy_db = zeros(number_of_points, 1);
offset = round(0.020 * fs);
for point_no = 1:number_of_points
    current = irs(point_no, :);
    [~, peak_index] = max(abs(current));
    peak_times_ms(point_no) = (peak_index - 1) / fs * 1000;
    late_start = min(numel(current) + 1, peak_index + offset);
    total_energy = sum(current.^2);
    late_energy = sum(current(late_start:end).^2);
    late_energy_db(point_no) = 10 * log10((late_energy + eps) / ...
        (total_energy + eps));
end
end

function manifest = select_manifest_rows(manifest, mode)
if mode == "quick"
    manifest = manifest(ismember(manifest.item_id, ["S01", "M06"]), :);
end
end

function perceptual = evaluate_audio_set(manifest, response_dir, ...
    result_dir, temporary_dir, renderers, reference, config)

audio_root = fullfile(response_dir, 'perceptual_audio');
listening_dir = fullfile(result_dir, 'rendered_listening');
ensure_directory(listening_dir);
algorithms = string(fieldnames(renderers));
source_rows = cell(0, 9);
speech_rows = cell(0, 8);
music_rows = cell(0, 7);
listening_rows = cell(0, 7);

ensure_pesq_python(config.python_executable);
peaq_dir = fullfile(response_dir, 'PEAQ');
original_directory = pwd;

for item_no = 1:height(manifest)
    manifest_row = manifest(item_no, :);
    [source, source_info] = prepare_source(manifest_row, audio_root, ...
        reference, config);
    source_rows(end + 1, :) = source_info_to_row(manifest_row, ...
        source_info); %#ok<AGROW>

    output_length = numel(source) + ...
        size(renderers.(char(algorithms(1))).effective_irs, 2) - 1;
    reference_bz = config.target_pressure_pa .* render_fractional_delays( ...
        source, reference.total_delay_samples(1:4), output_length);
    rendered = cell(numel(algorithms), 1);
    global_peak = max(abs(reference_bz(:)));
    for algorithm_no = 1:numel(algorithms)
        effective_irs = renderers.(char(algorithms(algorithm_no))).effective_irs;
        rendered{algorithm_no} = render_effective_irs(source, effective_irs);
        global_peak = max(global_peak, max(abs(rendered{algorithm_no}(:))));
    end
    common_gain = min(1, 0.95 / max(global_peak, eps));
    reference_digital = reference_bz .* common_gain;

    for algorithm_no = 1:numel(algorithms)
        algorithm = algorithms(algorithm_no);
        pressure_signals = rendered{algorithm_no};
        digital_signals = pressure_signals .* common_gain;

        if manifest_row.category == "speech"
            for point_no = 1:4
                orientation = char("horizontal");
                if point_no > 2
                    orientation = char("vertical");
                end
                reference_point = reference_digital(:, point_no);
                degraded_point = digital_signals(:, point_no);
                pesq_score = calculate_pesq(reference_point, ...
                    degraded_point, reference.fs);
                stoi_score = stoi(reference_point, degraded_point, ...
                    reference.fs);
                speech_rows(end + 1, :) = {char(manifest_row.item_id), ...
                    char(algorithm), reference.point_labels(point_no), ...
                    orientation, point_no, pesq_score, stoi_score, ...
                    common_gain}; %#ok<AGROW>
            end
        else
            cd(peaq_dir);
            for orientation_no = 1:2
                [orientation, bz_columns, ~] = ...
                    orientation_columns(orientation_no);
                [odg_raw, odg_clipped] = calculate_peaq( ...
                    reference_digital(:, bz_columns), ...
                    digital_signals(:, bz_columns), reference.fs, ...
                    temporary_dir, manifest_row.item_id, algorithm, ...
                    orientation);
                music_rows(end + 1, :) = {char(manifest_row.item_id), ...
                    char(algorithm), orientation, orientation_no, ...
                    odg_raw, odg_clipped, common_gain}; %#ok<AGROW>
            end
            cd(original_directory);
        end

        if ismember(manifest_row.item_id, config.listening_item_ids) && ...
                ismember(algorithm, config.listening_algorithms)
            listening_rows = export_algorithm_audio(listening_rows, ...
                listening_dir, manifest_row.item_id, algorithm, ...
                digital_signals, reference.fs, common_gain);
        end
    end

    if ismember(manifest_row.item_id, config.listening_item_ids)
        listening_rows = export_reference_audio(listening_rows, ...
            listening_dir, manifest_row.item_id, reference_digital, ...
            reference.fs, common_gain);
    end
    fprintf('  Audio item %d/%d complete: %s\n', item_no, ...
        height(manifest), manifest_row.item_id);
end
cd(original_directory);

[source_table, speech_table, music_table, listening_table] = ...
    make_audio_tables(source_rows, speech_rows, music_rows, listening_rows);
[speech_summary, music_summary] = summarize_perceptual_results( ...
    algorithms, speech_table, music_table);
writetable(source_table, fullfile(result_dir, 'source_preprocessing.csv'));
writetable(speech_table, fullfile(result_dir, ...
    'speech_pesq_stoi_detailed.csv'));
writetable(music_table, fullfile(result_dir, ...
    'music_peaq_detailed.csv'));
writetable(speech_summary, fullfile(result_dir, 'speech_summary.csv'));
writetable(music_summary, fullfile(result_dir, 'music_summary.csv'));
writetable(listening_table, fullfile(listening_dir, ...
    'listening_manifest.csv'));
save(fullfile(result_dir, 'perceptual_results.mat'), 'source_table', ...
    'speech_table', 'music_table', 'speech_summary', 'music_summary', ...
    'listening_table', '-v7.3');

perceptual = struct();
perceptual.sources = source_table;
perceptual.speech = speech_table;
perceptual.music = music_table;
perceptual.speech_summary = speech_summary;
perceptual.music_summary = music_summary;
perceptual.listening = listening_table;
fprintf('\nSpeech summary:\n');
disp(speech_summary);
fprintf('\nMusic summary:\n');
disp(music_summary);
end

function ensure_pesq_python(python_executable)
environment = pyenv;
if environment.Status == "NotLoaded" && exist(python_executable, 'file')
    pyenv('Version', python_executable, 'ExecutionMode', 'OutOfProcess');
end
try
    py.importlib.import_module('numpy');
    py.importlib.import_module('pesq');
catch exception
    warning(exception.identifier, '%s', exception.message);
end
end

function [source, info] = prepare_source(manifest_row, audio_root, ...
    reference, ~)
audio_file = fullfile(audio_root, manifest_row.category, ...
    manifest_row.copied_filename);
[source, original_fs] = audioread(audio_file);
first_sample = max(1, floor(manifest_row.excerpt_start_s * original_fs) + 1);
last_sample = min(size(source, 1), ...
    round(manifest_row.excerpt_end_s * original_fs));
source = source(first_sample:last_sample, :);
original_channels = size(source, 2);
if original_channels > 1
    source = mean(source, 2);
end
if original_fs ~= reference.fs
    [p, q] = rat(reference.fs / original_fs, 1e-12);
    source = resample(source, p, q);
end
source = double(source(:));
source = source - mean(source);
[b, a] = butter(6, 100 / (reference.fs / 2), 'high');
source = filtfilt(b, a, source);
active_rms_before = calculate_active_rms(source, reference.fs);
source = source ./ max(active_rms_before, eps);

info = struct();
info.original_fs = original_fs;
info.original_channels = original_channels;
info.output_samples = numel(source);
info.output_duration_s = numel(source) / reference.fs;
info.active_rms_before = active_rms_before;
info.active_rms_after = calculate_active_rms(source, reference.fs);
info.peak_after = max(abs(source));
info.was_resampled = original_fs ~= reference.fs;
end

function active_rms = calculate_active_rms(signal, fs)
frame_length = round(0.02 * fs);
hop_length = round(0.01 * fs);
if numel(signal) < frame_length
    active_rms = sqrt(mean(signal.^2));
    return;
end
number_of_frames = 1 + floor((numel(signal) - frame_length) / hop_length);
frame_rms = zeros(number_of_frames, 1);
for frame_no = 1:number_of_frames
    first = (frame_no - 1) * hop_length + 1;
    samples = signal(first:(first + frame_length - 1));
    frame_rms(frame_no) = sqrt(mean(samples.^2));
end
active = frame_rms >= max(frame_rms) * 10^(-40 / 20);
active_rms = sqrt(mean(frame_rms(active).^2));
end

function row = source_info_to_row(manifest_row, info)
row = {char(manifest_row.item_id), char(manifest_row.category), ...
    info.original_fs, info.original_channels, info.output_samples, ...
    info.output_duration_s, info.active_rms_before, ...
    info.active_rms_after, info.peak_after};
end

function output = render_fractional_delays(source, delays, output_length)
source = double(source(:));
fft_length = 2^nextpow2(output_length + 1024);
source_spectrum = fft(source, fft_length);
if mod(fft_length, 2) == 0
    signed_bins = [0:(fft_length / 2), ...
        (-fft_length / 2 + 1):-1].';
else
    half = floor(fft_length / 2);
    signed_bins = [0:half, -half:-1].';
end
output = zeros(output_length, numel(delays));
for delay_no = 1:numel(delays)
    response = exp(-1i * 2 * pi .* signed_bins .* ...
        delays(delay_no) / fft_length);
    if mod(fft_length, 2) == 0
        response(fft_length / 2 + 1) = cos(pi * delays(delay_no));
    end
    current = real(ifft(source_spectrum .* response));
    output(:, delay_no) = current(1:output_length);
end
end

function output = render_effective_irs(source, effective_irs)
output_length = numel(source) + size(effective_irs, 2) - 1;
fft_length = 2^nextpow2(output_length);
source_spectrum = fft(source, fft_length);
response_spectrum = fft(effective_irs, fft_length, 2);
output = real(ifft(response_spectrum .* source_spectrum.', [], 2)).';
output = output(1:output_length, :);
end

function [orientation, bz_columns, dz_columns] = ...
    orientation_columns(orientation_no)
if orientation_no == 1
    orientation = 'horizontal';
    bz_columns = [1, 2];
    dz_columns = [5, 6];
else
    orientation = 'vertical';
    bz_columns = [3, 4];
    dz_columns = [7, 8];
end
end

function score = calculate_pesq(reference, degraded, fs)
try
    reference = double(reference(:));
    degraded = double(degraded(:));
    score = double(py.pesq.pesq(int32(fs), ...
        py.numpy.array(reference.'), py.numpy.array(degraded.'), 'wb'));
catch exception
    warning(exception.identifier, '%s', exception.message);
    score = NaN;
end
end

function [odg_raw, odg_clipped] = calculate_peaq(reference, degraded, fs, ...
    temporary_dir, item_id, algorithm, orientation)
if fs ~= 16000
    error('PEAQ conversion expects a 16 kHz rendered signal.');
end
reference_48k = resample(reference, 3, 1);
degraded_48k = resample(degraded, 3, 1);
minimum_length = min(size(reference_48k, 1), size(degraded_48k, 1));
reference_48k = reference_48k(1:minimum_length, :);
degraded_48k = degraded_48k(1:minimum_length, :);

safe_algorithm = regexprep(char(algorithm), '[^A-Za-z0-9_]', '_');
reference_file = fullfile(temporary_dir, sprintf('%s_%s_%s_ref.wav', ...
    item_id, safe_algorithm, orientation));
degraded_file = fullfile(temporary_dir, sprintf('%s_%s_%s_deg.wav', ...
    item_id, safe_algorithm, orientation));
audiowrite(reference_file, reference_48k, 48000, 'BitsPerSample', 24);
audiowrite(degraded_file, degraded_48k, 48000, 'BitsPerSample', 24);
if ~exist(reference_file, 'file') || ~exist(degraded_file, 'file')
    error('Temporary PEAQ WAV files were not created.');
end
odg_raw = NaN;
for attempt_no = 1:2
    try
        evalc('odg_raw = PQevalAudio(reference_file, degraded_file);');
        break;
    catch exception
        if attempt_no == 2
            warning('PEAQ failed for %s/%s/%s after retry: %s', ...
                item_id, algorithm, orientation, exception.message);
        else
            pause(0.2);
        end
    end
end
if isnan(odg_raw)
    odg_clipped = NaN;
else
    odg_clipped = min(0, max(-4, odg_raw));
end
if exist(reference_file, 'file')
    delete(reference_file);
end
if exist(degraded_file, 'file')
    delete(degraded_file);
end
end

function rows = export_algorithm_audio(rows, output_dir, item_id, ...
    algorithm, signals, fs, common_gain)
for orientation_no = 1:2
    [orientation, bz_columns, dz_columns] = orientation_columns(orientation_no);
    zone_names = {'BZ', 'DZ'};
    columns = {bz_columns, dz_columns};
    for zone_no = 1:2
        filename = sprintf('%s_%s_%s_%s.wav', item_id, ...
            zone_names{zone_no}, orientation, algorithm);
        filename = regexprep(filename, '[^A-Za-z0-9_.-]', '_');
        audiowrite(fullfile(output_dir, filename), ...
            signals(:, columns{zone_no}), fs, 'BitsPerSample', 24);
        rows(end + 1, :) = {char(item_id), char(algorithm), ...
            zone_names{zone_no}, orientation, filename, fs, common_gain}; %#ok<AGROW>
    end
end
end

function rows = export_reference_audio(rows, output_dir, item_id, ...
    reference_signals, fs, common_gain)
for orientation_no = 1:2
    [orientation, bz_columns, ~] = orientation_columns(orientation_no);
    filename = sprintf('%s_BZ_%s_reference.wav', item_id, orientation);
    audiowrite(fullfile(output_dir, filename), ...
        reference_signals(:, bz_columns), fs, 'BitsPerSample', 24);
    rows(end + 1, :) = {char(item_id), 'reference', 'BZ', ...
        orientation, filename, fs, common_gain}; %#ok<AGROW>
end
end

function [source_table, speech_table, music_table, listening_table] = ...
    make_audio_tables(source_rows, speech_rows, music_rows, listening_rows)
source_table = cell2table(source_rows, 'VariableNames', ...
    {'ItemID', 'Category', 'OriginalFs', 'OriginalChannels', ...
    'OutputSamples', 'OutputDuration_s', 'ActiveRMSBefore', ...
    'ActiveRMSAfter', 'PeakAfter'});
speech_table = cell2table(speech_rows, 'VariableNames', ...
    {'ItemID', 'Algorithm', 'PointLabel', 'Orientation', ...
    'PointNumber', 'PESQ', 'STOI', 'CommonDigitalGain'});
music_table = cell2table(music_rows, 'VariableNames', ...
    {'ItemID', 'Algorithm', 'Orientation', 'PairNumber', ...
    'ODG_Raw', 'ODG_Clipped', 'CommonDigitalGain'});
listening_table = cell2table(listening_rows, 'VariableNames', ...
    {'ItemID', 'Algorithm', 'Zone', 'Orientation', 'Filename', ...
    'SampleRate_Hz', 'CommonDigitalGain'});
end

function [speech_summary, music_summary] = ...
    summarize_perceptual_results(algorithms, speech, music)
speech_rows = cell(numel(algorithms), 5);
music_rows = cell(numel(algorithms), 3);
for algorithm_no = 1:numel(algorithms)
    algorithm = algorithms(algorithm_no);
    s = speech(string(speech.Algorithm) == algorithm, :);
    m = music(string(music.Algorithm) == algorithm, :);
    s_h = s(string(s.Orientation) == "horizontal", :);
    s_v = s(string(s.Orientation) == "vertical", :);
    m_h = m(string(m.Orientation) == "horizontal", :);
    m_v = m(string(m.Orientation) == "vertical", :);
    speech_rows(algorithm_no, :) = {char(algorithm), ...
        mean(s_h.PESQ, 'omitnan'), mean(s_v.PESQ, 'omitnan'), ...
        mean(s_h.STOI, 'omitnan'), mean(s_v.STOI, 'omitnan')};
    music_rows(algorithm_no, :) = {char(algorithm), ...
        mean(m_h.ODG_Raw, 'omitnan'), mean(m_v.ODG_Raw, 'omitnan')};
end
speech_summary = cell2table(speech_rows, 'VariableNames', ...
    {'Algorithm', 'PESQ_Horizontal', 'PESQ_Vertical', ...
    'STOI_Horizontal', 'STOI_Vertical'});
music_summary = cell2table(music_rows, 'VariableNames', ...
    {'Algorithm', 'PEAQ_ODG_Horizontal', 'PEAQ_ODG_Vertical'});
end

function create_geometry_figure(room_array, point_definition, ...
    output_dir, config)
figure_handle = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, config.figure_size_cm]);
axes_handle = axes(figure_handle);
hold(axes_handle, 'on');
scatter(axes_handle, room_array.s(:, 1), room_array.s(:, 2), ...
    28, 'k', '^', 'filled', 'DisplayName', 'Loudspeakers');
theta = linspace(0, 2 * pi, 400);
plot(axes_handle, 2.6 + 0.2 * cos(theta), 2.0 + 0.2 * sin(theta), ...
    '-', 'Color', [0.00, 0.45, 0.74], 'LineWidth', 1.2, ...
    'DisplayName', 'Bright zone');
plot(axes_handle, 1.4 + 0.2 * cos(theta), 2.0 + 0.2 * sin(theta), ...
    '-', 'Color', [0.85, 0.33, 0.10], 'LineWidth', 1.2, ...
    'DisplayName', 'Dark zone');
scatter(axes_handle, point_definition.positions(1:4, 1), ...
    point_definition.positions(1:4, 2), 45, [0.00, 0.45, 0.74], ...
    'o', 'filled', 'DisplayName', 'BZ ear points');
scatter(axes_handle, point_definition.positions(5:8, 1), ...
    point_definition.positions(5:8, 2), 45, [0.85, 0.33, 0.10], ...
    's', 'filled', 'DisplayName', 'DZ ear points');
axis(axes_handle, 'equal');
xlim(axes_handle, [0.3, 3.7]);
ylim(axes_handle, [0.3, 3.7]);
xlabel(axes_handle, 'x [m]');
ylabel(axes_handle, 'y [m]');
title(axes_handle, 'Exact evaluation-point geometry', ...
    'FontSize', 11, 'FontWeight', 'normal');
style_axes(axes_handle, config);
legend(axes_handle, 'Location', 'southoutside', ...
    'Orientation', 'horizontal');
save_figure(figure_handle, output_dir, 'evaluation_point_geometry');
close(figure_handle);
end

function create_frequency_figures(table_data, output_dir, config)
algorithms = config.display_algorithms;
colors = lines(numel(algorithms));

figure_handle = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, config.figure_size_cm]);
axes_handle = axes(figure_handle);
hold(axes_handle, 'on');
for algorithm_no = 1:numel(algorithms)
    subset = table_data(string(table_data.Algorithm) == ...
        algorithms(algorithm_no), :);
    plot(axes_handle, subset.Frequency_Hz, subset.Planarity_percent, ...
        'LineWidth', 1.4, 'Color', colors(algorithm_no, :), ...
        'DisplayName', config.display_names(algorithm_no));
end
xline(axes_handle, config.aliasing_frequency_hz, '--k', ...
    'HandleVisibility', 'off');
xlim(axes_handle, [100, config.summary_upper_frequency_hz]);
ylim(axes_handle, [0, 100]);
xlabel(axes_handle, 'Frequency [Hz]');
ylabel(axes_handle, 'Planarity [%]');
style_axes(axes_handle, config);
legend(axes_handle, 'Location', 'southoutside', ...
    'Orientation', 'horizontal', 'NumColumns', 4);
save_figure(figure_handle, output_dir, 'planarity_by_algorithm');
close(figure_handle);

figure_handle = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 19, 18]);
layout = tiledlayout(figure_handle, 3, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
fields = {'BZ_Horizontal_LevelImbalance_dB', ...
    'BZ_Vertical_LevelImbalance_dB', ...
    'BZ_Horizontal_PhaseError_deg', 'BZ_Vertical_PhaseError_deg', ...
    'Horizontal_Attenuation_dB', 'Vertical_Attenuation_dB'};
ylabels = {'BZ horizontal level imbalance [dB]', ...
    'BZ vertical level imbalance [dB]', ...
    'BZ horizontal phase error [deg]', ...
    'BZ vertical phase error [deg]', ...
    'Horizontal BZ-DZ attenuation [dB]', ...
    'Vertical BZ-DZ attenuation [dB]'};
titles = {'BZ horizontal ear pair', 'BZ vertical ear pair', ...
    'BZ horizontal ear pair', 'BZ vertical ear pair', ...
    'Horizontal ear pairs', 'Vertical ear pairs'};
line_handles = gobjects(numel(algorithms), 1);
for panel_no = 1:numel(fields)
    axes_handle = nexttile(layout);
    hold(axes_handle, 'on');
    for algorithm_no = 1:numel(algorithms)
        subset = table_data(string(table_data.Algorithm) == ...
            algorithms(algorithm_no), :);
        line_handles(algorithm_no) = plot(axes_handle, ...
            subset.Frequency_Hz, subset.(fields{panel_no}), ...
            'LineWidth', 1.2, 'Color', colors(algorithm_no, :), ...
            'DisplayName', config.display_names(algorithm_no));
    end
    xline(axes_handle, config.aliasing_frequency_hz, '--k', ...
        'HandleVisibility', 'off');
    xlim(axes_handle, [100, config.summary_upper_frequency_hz]);
    xlabel(axes_handle, 'Frequency [Hz]');
    ylabel(axes_handle, ylabels{panel_no});
    title(axes_handle, titles{panel_no}, ...
        'FontSize', 11, 'FontWeight', 'normal');
    style_axes(axes_handle, config);
end
legend_handle = legend(line_handles, config.display_names, ...
    'Orientation', 'horizontal', 'NumColumns', 4);
legend_handle.Layout.Tile = 'south';
save_figure(figure_handle, output_dir, 'local_pair_frequency_metrics');
close(figure_handle);
end

function style_axes(axes_handle, config)
set(axes_handle, 'FontName', config.figure_font, 'FontSize', 10, ...
    'LineWidth', 0.8, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
    'XMinorGrid', 'on', 'YMinorGrid', 'on', 'GridAlpha', 0.24, ...
    'MinorGridAlpha', 0.12, 'GridLineStyle', '-', ...
    'MinorGridLineStyle', ':');
end

function save_figure(figure_handle, output_dir, basename)
exportgraphics(figure_handle, fullfile(output_dir, [basename, '.pdf']), ...
    'ContentType', 'vector');
exportgraphics(figure_handle, fullfile(output_dir, [basename, '.png']), ...
    'Resolution', 300);
savefig(figure_handle, fullfile(output_dir, [basename, '.fig']));
end

function remove_temporary_wavs(directory)
files = dir(fullfile(directory, '*.wav'));
for file_no = 1:numel(files)
    delete(fullfile(files(file_no).folder, files(file_no).name));
end
end

function ensure_directory(directory)
if ~exist(directory, 'dir')
    mkdir(directory);
end
end

function write_completion_marker(result_dir, mode, run_id, config)
marker_file = fullfile(result_dir, 'RUN_COMPLETE.txt');
file_id = fopen(marker_file, 'w');
if file_id < 0
    error('Unable to create completion marker: %s', marker_file);
end
cleanup = onCleanup(@() fclose(file_id));
fprintf(file_id, 'Run ID: %s\n', run_id);
fprintf(file_id, 'Mode: %s\n', mode);
fprintf(file_id, 'Design/calibration temperature: %.2f C\n', ...
    config.calibration_temperature_celsius);
fprintf(file_id, 'Evaluation temperature: %.2f C\n', ...
    config.evaluation_temperature_celsius);
fprintf(file_id, 'Filter source: %s\n', config.filter_file);
fprintf(file_id, 'Primary rho: %.6g\n', config.primary_rho);
fprintf(file_id, 'Primary nu: %.6g\n', config.primary_nu);
fprintf(file_id, 'Calibration: %s\n', config.calibration_protocol);
if strlength(config.method_note) > 0
    fprintf(file_id, 'Filter method: %s\n', config.method_note);
end
fprintf(file_id, 'Algorithms: %s\n', ...
    strjoin(cellstr(config.algorithm_fields), ', '));
fprintf(file_id, 'Completed: %s\n', ...
    string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
clear cleanup;
end

function tag = make_temperature_tag(temperature_celsius)
tag = sprintf('T%.2fC', temperature_celsius);
tag = strrep(tag, '-', 'm');
tag = strrep(tag, '.', 'p');
end

function folder = get_last_folder(path_value)
[~, folder] = fileparts(path_value);
end
