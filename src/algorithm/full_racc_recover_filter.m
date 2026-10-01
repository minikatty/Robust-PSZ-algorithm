function [w_best, recovery] = full_racc_recover_filter( ...
        W, HB, HD, etaB, etaD, options, freq_idx, deltaB, deltaD)
%FULL_RACC_RECOVER_FILTER Recover a realizable rank-one filter from the SDP.
%
% recovery_method='paper' (used by the official Step-2 experiment) follows
% the manuscript exactly: use the principal eigenvector when W is numerically
% rank one; otherwise apply Gaussian randomization and retain the feasible
% candidate with the largest exact ATF-ball worst-case contrast.
%
% recovery_method='enhanced' preserves the earlier engineering diagnostic:
% it additionally tries nominal-ACC/matched-RACC seeds and local complex-
% sphere refinement. This mode must be labelled separately if used because
% that refinement is not part of the stated SDP rank-recovery procedure.

    if nargin < 7 || isempty(freq_idx)
        freq_idx = 1;
    end
    if nargin < 8 || isempty(deltaB)
        deltaB = 0;
    end
    if nargin < 9 || isempty(deltaD)
        deltaD = 0;
    end
    defaults = struct('recovery_method', 'paper', ...
        'randomization_trials', 300, ...
        'random_seed', 20260820, 'refinement_iterations', 300, ...
        'refinement_seeds', 8, 'refinement_tolerance', 1e-8);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(options, names{k})
            options.(names{k}) = defaults.(names{k});
        end
    end
    recovery_method = lower(char(string(options.recovery_method)));
    if ~ismember(recovery_method, {'paper', 'enhanced'})
        error('FullRACC:UnknownRecoveryMethod', ...
            'recovery_method must be ''paper'' or ''enhanced''.');
    end

    W = (W + W') / 2;
    [V, D] = eig(W);
    eigenvalues = real(diag(D));
    [eigenvalues, order] = sort(eigenvalues, 'descend');
    V = V(:, order);
    positive_values = max(eigenvalues, 0);
    positive_sum = sum(positive_values);
    if positive_sum <= 0
        error('FullRACC:InvalidLiftedSolution', ...
            'The lifted solution has no positive eigenvalue.');
    end

    rank_fraction = positive_values(1) / positive_sum;
    if numel(positive_values) >= 2
        rank_ratio = positive_values(1) / max(positive_values(2), eps);
    else
        rank_ratio = Inf;
    end

    seed_vectors = cell(0, 1);
    seed_sources = strings(0, 1);
    seed_alpha = zeros(0, 1);
    [seed_vectors, seed_sources, seed_alpha] = local_add_seed( ...
        seed_vectors, seed_sources, seed_alpha, V(:, 1), ...
        "SDP principal eigenvector", HB, HD, etaB, etaD, deltaB, deltaD);

    if strcmp(recovery_method, 'enhanced')
        [w_acc, w_matched] = local_deterministic_seeds( ...
            HB, HD, etaB, etaD, deltaB, deltaD);
        [seed_vectors, seed_sources, seed_alpha] = local_add_seed( ...
            seed_vectors, seed_sources, seed_alpha, w_acc, ...
            "Nominal ACC seed", HB, HD, etaB, etaD, deltaB, deltaD);
        [seed_vectors, seed_sources, seed_alpha] = local_add_seed( ...
            seed_vectors, seed_sources, seed_alpha, w_matched, ...
            "Matched RACC seed", HB, HD, etaB, etaD, deltaB, deltaD);
    end

    if rank_fraction < 1 - 1e-6 && options.randomization_trials > 0
        covariance_root = V * diag(sqrt(positive_values));
        stream = RandStream('mt19937ar', 'Seed', options.random_seed + freq_idx);
        random_vectors = cell(options.randomization_trials, 1);
        random_alpha = -inf(options.randomization_trials, 1);
        for trial = 1:options.randomization_trials
            z = (randn(stream, size(W, 1), 1) + ...
                1i * randn(stream, size(W, 1), 1)) / sqrt(2);
            candidate = covariance_root * z;
            candidate_norm = norm(candidate);
            if candidate_norm <= 1e-12
                continue;
            end
            candidate = candidate / candidate_norm;
            [~, random_alpha(trial)] = full_racc_worst_case_contrast( ...
                candidate, HB, HD, etaB, etaD, deltaB, deltaD);
            random_vectors{trial} = candidate;
        end
        [~, random_order] = sort(random_alpha, 'descend');
        if strcmp(recovery_method, 'paper')
            keep_count = min(1, options.randomization_trials);
        else
            keep_count = min(options.refinement_seeds, ...
                options.randomization_trials);
        end
        for k = 1:keep_count
            trial = random_order(k);
            if isempty(random_vectors{trial}) || ~isfinite(random_alpha(trial))
                continue;
            end
            [seed_vectors, seed_sources, seed_alpha] = local_add_seed( ...
                seed_vectors, seed_sources, seed_alpha, random_vectors{trial}, ...
                "Gaussian randomization", HB, HD, etaB, etaD, ...
                deltaB, deltaD);
        end
    end

    [best_seed_alpha, best_idx] = max(seed_alpha);
    w_best = seed_vectors{best_idx};
    best_alpha = best_seed_alpha;
    best_source = seed_sources(best_idx);
    best_refinement_iterations = 0;

    if strcmp(recovery_method, 'enhanced')
        [~, seed_order] = sort(seed_alpha, 'descend');
        refine_count = min(numel(seed_order), options.refinement_seeds);
        for k = 1:refine_count
            idx = seed_order(k);
            [candidate, candidate_alpha, iterations] = ...
                local_refine_on_sphere(seed_vectors{idx}, HB, HD, ...
                etaB, etaD, deltaB, deltaD, options);
            if candidate_alpha > best_alpha
                w_best = candidate;
                best_alpha = candidate_alpha;
                best_source = seed_sources(idx) + " + sphere refinement";
                best_refinement_iterations = iterations;
            end
        end
    end

    recovery = struct('rank_fraction', rank_fraction, ...
        'rank_ratio', rank_ratio, 'selected_alpha', best_alpha, ...
        'best_seed_alpha', best_seed_alpha, ...
        'source', char(best_source), ...
        'method', recovery_method, ...
        'refinement_iterations', best_refinement_iterations, ...
        'num_seeds', numel(seed_vectors));
end

function [vectors, sources, scores] = local_add_seed( ...
        vectors, sources, scores, w, source, HB, HD, etaB, etaD, ...
        deltaB, deltaD)
    if isempty(w) || norm(w) <= 1e-12
        return;
    end
    w = w / norm(w);
    [~, alpha] = full_racc_worst_case_contrast( ...
        w, HB, HD, etaB, etaD, deltaB, deltaD);
    vectors{end + 1, 1} = w; %#ok<AGROW>
    sources(end + 1, 1) = source; %#ok<AGROW>
    scores(end + 1, 1) = alpha; %#ok<AGROW>
end

function [w_acc, w_matched] = local_deterministic_seeds( ...
        HB, HD, etaB, etaD, deltaB, deltaD)
    L = size(HB, 2);
    RB = (HB' * HB + (HB' * HB)') / 2;
    RD = (HD' * HD + (HD' * HD)') / 2;
    floor_value = max(1e-12 * max(real(trace(RD)) / L, 1), eps);

    [V, D] = eig(RB + deltaB * eye(L), ...
        RD + (deltaD + floor_value) * eye(L));
    [~, idx] = max(real(diag(D)));
    w_acc = V(:, idx);
    w_acc = w_acc / norm(w_acc);

    A = RB + (deltaB - etaB^2) * eye(L);
    B = RD + (deltaD + etaD^2 + floor_value) * eye(L);
    [V, D] = eig((A + A') / 2, (B + B') / 2);
    [~, idx] = max(real(diag(D)));
    w_matched = V(:, idx);
    w_matched = w_matched / norm(w_matched);
end

function [w, alpha, iterations] = local_refine_on_sphere( ...
        w, HB, HD, etaB, etaD, deltaB, deltaD, options)
    RB = HB' * HB;
    RD = HD' * HD;
    w = w / norm(w);
    [objective, gradient] = local_log_objective_gradient( ...
        w, RB, RD, etaB, etaD, deltaB, deltaD);
    iterations = 0;
    if ~isfinite(objective)
        [~, alpha] = full_racc_worst_case_contrast( ...
            w, HB, HD, etaB, etaD, deltaB, deltaD);
        return;
    end

    for iteration = 1:options.refinement_iterations
        tangent = gradient - w * real(w' * gradient);
        tangent_norm = norm(tangent);
        if tangent_norm <= options.refinement_tolerance
            break;
        end
        direction = tangent / tangent_norm;
        directional_derivative = 2 * real(gradient' * direction);
        step = 1;
        accepted = false;
        for backtrack = 1:30
            candidate = w + step * direction;
            candidate = candidate / norm(candidate);
            candidate_objective = local_log_objective_gradient( ...
                candidate, RB, RD, etaB, etaD, deltaB, deltaD);
            if candidate_objective >= objective + ...
                    1e-4 * step * directional_derivative
                accepted = true;
                break;
            end
            step = step / 2;
        end
        if ~accepted
            break;
        end
        improvement = candidate_objective - objective;
        w = candidate;
        objective = candidate_objective;
        [~, gradient] = local_log_objective_gradient( ...
            w, RB, RD, etaB, etaD, deltaB, deltaD);
        iterations = iteration;
        if improvement <= options.refinement_tolerance
            break;
        end
    end
    [~, alpha] = full_racc_worst_case_contrast( ...
        w, HB, HD, etaB, etaD, deltaB, deltaD);
end

function [objective, gradient] = local_log_objective_gradient( ...
        w, RB, RD, etaB, etaD, deltaB, deltaD)
    w_norm = norm(w);
    bright_norm = sqrt(max(real(w' * RB * w), 0));
    dark_norm = sqrt(max(real(w' * RD * w), 0));
    bright_margin_raw = bright_norm - etaB * w_norm;
    bright_margin = max(bright_margin_raw, 0);
    dark_margin = dark_norm + etaD * w_norm;
    bright_energy = bright_margin^2 + deltaB * w_norm^2;
    dark_energy = dark_margin^2 + deltaD * w_norm^2;
    if bright_energy <= 0 || dark_energy <= 0 || w_norm <= 0
        objective = -Inf;
        gradient = zeros(size(w));
        return;
    end
    objective = log(bright_energy) - log(dark_energy);
    if bright_margin_raw > 0 && bright_norm > 0
        bright_data_gradient = RB * w / bright_norm - etaB * w / w_norm;
    else
        bright_data_gradient = zeros(size(w));
    end
    bright_gradient = ...
        (bright_margin * bright_data_gradient + deltaB * w) / bright_energy;
    if dark_norm > 0
        dark_data_gradient = RD * w / dark_norm + etaD * w / w_norm;
    else
        dark_data_gradient = etaD * w / w_norm;
    end
    dark_gradient = ...
        (dark_margin * dark_data_gradient + deltaD * w) / dark_energy;
    gradient = bright_gradient - dark_gradient;
end
