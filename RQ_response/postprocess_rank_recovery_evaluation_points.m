function outputs = postprocess_rank_recovery_evaluation_points(rank_run_dir)
%POSTPROCESS_RANK_RECOVERY_EVALUATION_POINTS Audit recovery on monitor grid.
%   This is a read-only post-processing step for a completed
%   run_racc_pm_rank_recovery_diagnostic result. It does not solve an SDP.
%   The saved relaxed matrix and deployed rank-one filter are evaluated on
%   the independent nominal BZ/DZ evaluation grid. AC is scale invariant.
%   For protocol-aligned NSRE, the relaxed representation and deployed
%   filter each receive their own energy-matching/global-phase gain determined
%   exclusively from nominal control points. Both gains are then frozen for
%   the common independent evaluation grid.

    response_dir = fileparts(mfilename('fullpath'));
    project_root = fileparts(response_dir);
    if nargin < 1 || isempty(rank_run_dir)
        rank_run_dir = fullfile(response_dir, ...
            'racc_pm_rank_recovery_results', ...
            'fullband_nu0p01_20260907_231815');
    end
    rank_run_dir = char(string(rank_run_dir));
    result_file = fullfile(rank_run_dir, 'rank_recovery_results.mat');
    if ~isfile(result_file)
        error('RankEvalPostprocess:MissingResult', ...
            'Missing completed rank result: %s', result_file);
    end

    addpath(fullfile(project_root, 'src', 'utils'));
    loaded = load(result_file, 'results');
    results = loaded.results;
    frequencies_hz = results.frequencies_hz(:);
    frequency_indices = results.frequency_indices(:);
    num_frequencies = numel(frequencies_hz);

    data_dir = fullfile(project_root, 'data', 'SimulateRIR', 'temperature');
    metadata = load(fullfile(data_dir, 'para.mat'), 'para');
    freq_params = metadata.para.freq_params;
    nominal_file = get_data_filename(data_dir, 'temperature', ...
        results.config.nominal_temperature_c);
    nominal = load(nominal_file, 'ATF_BZ', 'ATF_DZ');
    HB_control = nominal.ATF_BZ.ctrl(:, :, frequency_indices);
    HB_evaluation = nominal.ATF_BZ.eval(:, :, frequency_indices);
    HD_evaluation = nominal.ATF_DZ.eval(:, :, frequency_indices);

    geometry = load(fullfile(project_root, 'data', 'arrayGeometry', ...
        'array_layout.mat'), 'roomArray');
    if isfield(freq_params, 'virtual_src_idx')
        virtual_source_index = freq_params.virtual_src_idx;
    else
        % The simulation protocol uses loudspeaker 13 to define the target
        % propagation direction. The saved control target below is used as
        % an independent numerical check of this metadata fallback.
        virtual_source_index = 13;
    end
    position = struct('src', geometry.roomArray.s, ...
        'pos_center', [2.6, 2.0, 1.6], ...
        'virtual_src_idx', virtual_source_index);
    position.mics_pos = geometry.roomArray.BZ_ctrl;
    desired_control_all = local_plane_wave_generator(position, freq_params);
    desired_control = desired_control_all(:, frequency_indices);
    saved_desired = load(fullfile(project_root, 'data', ...
        'ATF_desired_plane.mat'), 'ATF_desired_plane');
    stored_control = ...
        saved_desired.ATF_desired_plane(:, frequency_indices);
    desired_mismatch = max(abs(desired_control(:) - stored_control(:)));
    if desired_mismatch > 1e-10
        error('RankEvalPostprocess:DesiredMismatch', ...
            'Generated control target differs from saved target by %.3g.', ...
            desired_mismatch);
    end
    position.mics_pos = geometry.roomArray.BZ_eval;
    desired_evaluation_all = local_plane_wave_generator(position, freq_params);
    desired_evaluation = desired_evaluation_all(:, frequency_indices);

    FrequencyHz = frequencies_hz;
    RecoveryMode = strings(num_frequencies, 1);
    RelaxedAlignmentGainMagnitude = nan(num_frequencies, 1);
    RecoveredAlignmentGainMagnitude = nan(num_frequencies, 1);
    RelaxedEvaluationACdB = nan(num_frequencies, 1);
    RecoveredEvaluationACdB = nan(num_frequencies, 1);
    RecoveredMinusRelaxedEvaluationACdB = nan(num_frequencies, 1);
    RawRelaxedEvaluationNSREdB = nan(num_frequencies, 1);
    RawRecoveredEvaluationNSREdB = nan(num_frequencies, 1);
    RawRecoveredMinusRelaxedEvaluationNSREdB = nan(num_frequencies, 1);
    AlignedRelaxedEvaluationNSREdB = nan(num_frequencies, 1);
    AlignedRecoveredEvaluationNSREdB = nan(num_frequencies, 1);
    AlignedRecoveredMinusRelaxedEvaluationNSREdB = nan(num_frequencies, 1);

    for frequency_no = 1:num_frequencies
        solver_output = results.solver_outputs{frequency_no};
        if ~isfield(solver_output, 'W_tilde') || ...
                ~isfield(solver_output, 'recovery')
            error('RankEvalPostprocess:IncompleteFrequency', ...
                'Missing lifted/recovered result at %.2f Hz.', ...
                frequencies_hz(frequency_no));
        end
        W_tilde = (solver_output.W_tilde + solver_output.W_tilde') / 2;
        L = size(HB_evaluation, 2);
        W = W_tilde(1:L, 1:L);
        lifted_w = W_tilde(1:L, L + 1);
        recovered_w = solver_output.recovery.production_filter;
        RecoveryMode(frequency_no) = string( ...
            solver_output.recovery.production_mode);

        HBc = HB_control(:, :, frequency_no);
        HBe = HB_evaluation(:, :, frequency_no);
        HDe = HD_evaluation(:, :, frequency_no);
        dc = desired_control(:, frequency_no);
        de = desired_evaluation(:, frequency_no);

        RelaxedEvaluationACdB(frequency_no) = ...
            local_relaxed_ac_db(W, HBe, HDe);
        RecoveredEvaluationACdB(frequency_no) = ...
            local_filter_ac_db(recovered_w, HBe, HDe);
        RecoveredMinusRelaxedEvaluationACdB(frequency_no) = ...
            RecoveredEvaluationACdB(frequency_no) - ...
            RelaxedEvaluationACdB(frequency_no);

        RawRelaxedEvaluationNSREdB(frequency_no) = ...
            local_relaxed_nsre_db(W, lifted_w, HBe, de);
        RawRecoveredEvaluationNSREdB(frequency_no) = ...
            local_filter_nsre_db(recovered_w, HBe, de);
        RawRecoveredMinusRelaxedEvaluationNSREdB(frequency_no) = ...
            RawRecoveredEvaluationNSREdB(frequency_no) - ...
            RawRelaxedEvaluationNSREdB(frequency_no);

        relaxed_gain = local_relaxed_alignment_gain(W, lifted_w, HBc, dc);
        recovered_gain = local_alignment_gain(HBc * recovered_w, dc);
        RelaxedAlignmentGainMagnitude(frequency_no) = abs(relaxed_gain);
        RecoveredAlignmentGainMagnitude(frequency_no) = ...
            abs(recovered_gain);
        aligned_W = abs(relaxed_gain)^2 * W;
        aligned_lifted_w = relaxed_gain * lifted_w;
        aligned_recovered_w = recovered_gain * recovered_w;
        AlignedRelaxedEvaluationNSREdB(frequency_no) = ...
            local_relaxed_nsre_db(aligned_W, aligned_lifted_w, HBe, de);
        AlignedRecoveredEvaluationNSREdB(frequency_no) = ...
            local_filter_nsre_db(aligned_recovered_w, HBe, de);
        AlignedRecoveredMinusRelaxedEvaluationNSREdB(frequency_no) = ...
            AlignedRecoveredEvaluationNSREdB(frequency_no) - ...
            AlignedRelaxedEvaluationNSREdB(frequency_no);
    end

    summary = table(FrequencyHz, RecoveryMode, ...
        RelaxedAlignmentGainMagnitude, RecoveredAlignmentGainMagnitude, ...
        RelaxedEvaluationACdB, RecoveredEvaluationACdB, ...
        RecoveredMinusRelaxedEvaluationACdB, ...
        RawRelaxedEvaluationNSREdB, RawRecoveredEvaluationNSREdB, ...
        RawRecoveredMinusRelaxedEvaluationNSREdB, ...
        AlignedRelaxedEvaluationNSREdB, ...
        AlignedRecoveredEvaluationNSREdB, ...
        AlignedRecoveredMinusRelaxedEvaluationNSREdB);
    band_summary = local_band_summary(summary);
    gaussian_bins = summary(summary.RecoveryMode == ...
        "Gaussian randomization", :);

    output_dir = fullfile(rank_run_dir, ...
        'evaluation_point_postprocess_independent_alignment');
    if ~isfolder(output_dir), mkdir(output_dir); end
    writetable(summary, fullfile(output_dir, ...
        'evaluation_point_rank_recovery_summary.csv'));
    writetable(band_summary, fullfile(output_dir, ...
        'evaluation_point_rank_recovery_band_summary.csv'));
    writetable(gaussian_bins, fullfile(output_dir, ...
        'evaluation_point_gaussian_bins.csv'));
    outputs = struct('source_result', result_file, ...
        'nominal_file', nominal_file, 'summary', summary, ...
        'band_summary', band_summary, 'gaussian_bins', gaussian_bins, ...
        'desired_control_validation_error', desired_mismatch, ...
        'output_dir', output_dir);
    save(fullfile(output_dir, ...
        'evaluation_point_rank_recovery_postprocess.mat'), ...
        'outputs', '-v7.3');
    disp(band_summary);
    disp(gaussian_bins(:, {'FrequencyHz', 'RecoveryMode', ...
        'RecoveredMinusRelaxedEvaluationACdB', ...
        'AlignedRecoveredMinusRelaxedEvaluationNSREdB'}));
    fprintf('Completed evaluation-point recovery post-processing: %s\n', ...
        output_dir);
end

function gain = local_relaxed_alignment_gain(W, w, HB, desired)
    response_energy = max(real(trace(HB * W * HB')), realmin);
    target_energy = max(norm(desired)^2, realmin);
    cross_term = w' * HB' * desired;
    if abs(cross_term) <= realmin
        phase_factor = 1;
    else
        phase_factor = exp(1i * angle(cross_term));
    end
    gain = sqrt(target_energy / response_energy) * phase_factor;
end

function gain = local_alignment_gain(response, desired)
    if norm(response) <= realmin || norm(desired) <= realmin
        error('RankEvalPostprocess:DegenerateAlignment', ...
            'Cannot align a degenerate response or target.');
    end
    cross_term = response' * desired;
    if abs(cross_term) <= realmin
        phase_factor = 1;
    else
        phase_factor = exp(1i * angle(cross_term));
    end
    gain = (norm(desired) / norm(response)) * phase_factor;
end

function value = local_relaxed_ac_db(W, HB, HD)
    bright_energy = max(real(trace(HB * W * HB')), 0);
    dark_energy = max(real(trace(HD * W * HD')), 0);
    ratio = size(HD, 1) * bright_energy / max( ...
        size(HB, 1) * dark_energy, realmin);
    value = 10 * log10(max(real(ratio), realmin));
end

function value = local_filter_ac_db(w, HB, HD)
    ratio = size(HD, 1) * norm(HB * w)^2 / max( ...
        size(HB, 1) * norm(HD * w)^2, realmin);
    value = 10 * log10(max(real(ratio), realmin));
end

function value = local_relaxed_nsre_db(W, w, HB, desired)
    generated_energy = max(real(trace(HB * W * HB')), 0);
    error_energy = real(generated_energy - ...
        2 * desired' * HB * w + desired' * desired);
    value = 10 * log10(max(max(error_energy, 0) / ...
        max(norm(desired)^2, realmin), realmin));
end

function value = local_filter_nsre_db(w, HB, desired)
    value = 10 * log10(max(norm(HB * w - desired)^2 / ...
        max(norm(desired)^2, realmin), realmin));
end

function desired = local_plane_wave_generator(position, freq_params)
    source_position = position.src(position.virtual_src_idx, :);
    propagation = position.pos_center - source_position;
    propagation = propagation / norm(propagation);
    wave_numbers = (2 * pi .* freq_params.target_freqs(:)) / 343;
    wave_vectors = wave_numbers .* propagation;
    phase_delay = position.mics_pos * wave_vectors.';
    desired = exp(-1i * phase_delay);
end

function band_summary = local_band_summary(summary)
    definitions = { ...
        "Full_100_4000", 100, 4000; ...
        "Primary_100_1600", 100, 1600; ...
        "Stress_1650_4000", 1650, 4000};
    rows = cell(size(definitions, 1), 1);
    for row_no = 1:size(definitions, 1)
        name = definitions{row_no, 1};
        lower_hz = definitions{row_no, 2};
        upper_hz = definitions{row_no, 3};
        mask = summary.FrequencyHz >= lower_hz & ...
            summary.FrequencyHz <= upper_hz;
        band = summary(mask, :);
        ac_gap = abs(band.RecoveredMinusRelaxedEvaluationACdB);
        raw_nsre_gap = abs( ...
            band.RawRecoveredMinusRelaxedEvaluationNSREdB);
        aligned_nsre_gap = abs( ...
            band.AlignedRecoveredMinusRelaxedEvaluationNSREdB);
        gaussian_count = sum(band.RecoveryMode == ...
            "Gaussian randomization");
        rows{row_no} = table(name, lower_hz, upper_hz, height(band), ...
            gaussian_count, median(ac_gap), max(ac_gap), ...
            median(raw_nsre_gap), max(raw_nsre_gap), ...
            median(aligned_nsre_gap), max(aligned_nsre_gap), ...
            'VariableNames', {'Band', 'LowerHz', 'UpperHz', ...
            'FrequencyCount', 'GaussianRecoveryCount', ...
            'MedianAbsoluteEvaluationACGapdB', ...
            'MaximumAbsoluteEvaluationACGapdB', ...
            'MedianAbsoluteRawEvaluationNSREGapdB', ...
            'MaximumAbsoluteRawEvaluationNSREGapdB', ...
            'MedianAbsoluteAlignedEvaluationNSREGapdB', ...
            'MaximumAbsoluteAlignedEvaluationNSREGapdB'});
    end
    band_summary = vertcat(rows{:});
end
