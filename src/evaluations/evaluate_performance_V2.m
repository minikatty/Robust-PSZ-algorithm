% src/evaluation/evaluate_performance.m

function performance = evaluate_performance_V2(filters_w, operating_data, freq_params, eval_params)
% EVALUATE_PERFORMANCE Calculates all performance metrics for a given set of filters.
%
%   This function is designed to be robust and algorithm-aware. It conditionally
%   calculates metrics like NSRE only for relevant algorithms and provides
%   sensible defaults for critical parameters.
%
% Inputs:
%   filters_w       - Control filters [N_loudspeakers x N_target_freqs]. The filters to be evaluated.
%   algorithm_name  - String. The name of the algorithm that generated the filters (e.g., 'ACC', 'PM').
%   operating_data  - Struct. Contains the ATFs for the operating condition (H_B, H_D).
%   freq_params     - Struct. Contains frequency-related parameters (e.g., target_indices).
%   eval_params     - Struct. Contains parameters specific to evaluation (e.g., S_matrix).

% --- 1. SETUP & DATA EXTRACTION ---
% Extract necessary data structures for clarity.
ATF_BZ_eval = operating_data.ATF_BZ.eval; %eval,ctrl
ATF_DZ_eval = operating_data.ATF_DZ.eval; %eval,ctrl
ATF_BZ_ctrl = operating_data.ATF_BZ.ctrl; % Used for defining the reference sound field.

target_indices   = freq_params.target_indices;
num_target_freqs = length(target_indices);
algorithm_name = eval_params.algorithm_name;
% --- 2. VALIDATE PARAMETERS & SET DEFAULTS ---
% This section makes the function robust to missing optional parameters.

% The reference virtual source index is required for AE (all algorithms) and
% NSRE (reproduction algorithms). Default to 1 if not provided.
if isfield(eval_params, 'virtual_source_idx') && ~isempty(eval_params.virtual_source_idx)
    ref_vsrc_idx = eval_params.virtual_source_idx;
else
    ref_vsrc_idx = 1; % Default to the first loudspeaker.
end

% The S matrix for planarity is mandatory. Error out if not provided.
if ~isfield(eval_params, 'S_matrix')
    error('evaluate_performance:MissingParameter', ...
          'S_matrix for Planarity calculation is not provided in eval_params.');
end
S_matrix_all_freqs = eval_params.S_matrix;


% --- 3. ALGORITHM-SPECIFIC LOGIC ---
% Determine if the current algorithm is reproduction-based. This controls
% whether NSRE is calculated. Add new reproduction algorithms to this list.
is_reproduction_based = any(strcmpi(algorithm_name, {'PM', 'RPM','ACC-PM', 'RACC-PM', 'POTDC_RACC'}));

if is_reproduction_based
    % The target sound field p_d is only needed for reproduction algorithms.
    % It's defined by the reference virtual source at the control points.
    p_d_all_freqs = ATF_BZ_ctrl(:, ref_vsrc_idx, :);
end


% --- 4. INITIALIZE RESULT VECTORS ---
% Pre-allocate memory for all potential metrics.
AC        = zeros(num_target_freqs, 1);
AE        = zeros(num_target_freqs, 1);
Planarity = zeros(num_target_freqs, 1);

% Conditionally initialize NSRE. For non-reproduction algorithms, it will be
% populated with NaN (Not a Number), clearly indicating it's not applicable.
if is_reproduction_based
    NSRE = zeros(num_target_freqs, 1);
else
    NSRE = NaN(num_target_freqs, 1);
end


% --- 5. MAIN CALCULATION LOOP ---
% Iterate through each target frequency to compute metrics.
for i = 1:num_target_freqs
    f_idx = target_indices(i);
    % Extract frequency-dependent data for this iteration.
    w_f   = filters_w(:, i);
    H_B_f = ATF_BZ_eval(:, :, f_idx);
    H_D_f = ATF_DZ_eval(:, :, f_idx);
    S_f   = S_matrix_all_freqs(:, :, i);

    % --- Calculate all applicable metrics ---
    AC(i) = calculate_AC(w_f, H_B_f, H_D_f);
    
    % AE is calculated for ALL algorithms for a fair energy comparison.
    % It uses the same reference source for normalization.
    AE(i) = calculate_AE(w_f, H_B_f, ref_vsrc_idx);

    Planarity(i) = calculate_planarity(w_f, H_B_f, S_f);

    % Only calculate NSRE for reproduction-based algorithms.
    if is_reproduction_based
        p_d_f  = p_d_all_freqs(:, 1, f_idx);
        NSRE(i) = calculate_NSRE(w_f, H_B_f, p_d_f);
    end
end

% --- 6. FINALIZE & PACKAGE OUTPUT ---
% Aggregate results into the output structure.
% Averaging with 'omitnan' ensures that NSRE is handled correctly for all cases.
performance.AC        = real(AC);
% figure;plot(target_freqs,real(AC))
performance.NSRE      = NSRE; % Returns NaN for ACC, value for PM.
performance.AE        = real(AE);
performance.Planarity = Planarity;

end