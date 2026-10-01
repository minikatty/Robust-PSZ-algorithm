function result = run_racc_pm_sub_parallel_checkpointed( ...
        HB, HD, desired, para, checkpoint_dir, parallel_workers)
%RUN_RACC_PM_SUB_PARALLEL_CHECKPOINTED Run maintained RACC_PM_Sub safely.
%   Each frequency is solved by the no-suffix RACC_PM_Sub implementation and
%   saved to a separate checkpoint. Process workers isolate CVX/MOSEK state,
%   while the worker count and MOSEK threads are kept deliberately small to
%   avoid the out-of-memory failures observed for decomposed LMIs.

    arguments
        HB (:,:,:) double
        HD (:,:,:) double
        desired (:,:) double
        para (1,1) struct
        checkpoint_dir (1,:) char
        parallel_workers (1,1) double {mustBeInteger, mustBePositive} = 2
    end

    implementation_path = which('RACC_PM_Sub');
    project_root = fileparts(fileparts(mfilename('fullpath')));
    expected_path = fullfile(project_root, 'src', 'algorithm', ...
        'RACC_PM_Sub.m');
    if ~strcmpi(implementation_path, expected_path)
        error('RACCPMCheckpoint:WrongImplementation', ...
            'Expected %s but MATLAB resolves RACC_PM_Sub to %s.', ...
            expected_path, implementation_path);
    end

    if ~exist(checkpoint_dir, 'dir')
        mkdir(checkpoint_dir);
    end

    [~, loudspeakers, num_frequencies] = size(HB);
    if size(HD, 3) ~= num_frequencies || ...
            size(desired, 2) ~= num_frequencies
        error('RACCPMCheckpoint:FrequencyDimensionMismatch', ...
            'HB, HD, and desired data must contain the same frequencies.');
    end

    filters = zeros(loudspeakers, num_frequencies);
    status = strings(1, num_frequencies);
    normalization_multiplier = nan(1, num_frequencies);
    normalization_attempt_status = cell(1, num_frequencies);
    pending = zeros(1, num_frequencies);
    pending_count = 0;

    for frequency_index = 1:num_frequencies
        checkpoint_file = local_checkpoint_file( ...
            checkpoint_dir, frequency_index);
        if exist(checkpoint_file, 'file')
            loaded = load(checkpoint_file, 'checkpoint');
            local_validate_checkpoint(loaded.checkpoint, para, ...
                frequency_index);
            [filters, status, normalization_multiplier, ...
                normalization_attempt_status] = local_insert_checkpoint( ...
                filters, status, normalization_multiplier, ...
                normalization_attempt_status, loaded.checkpoint);
            fprintf('RACC-PM frequency %d/%d restored from checkpoint.\n', ...
                frequency_index, num_frequencies);
        else
            pending_count = pending_count + 1;
            pending(pending_count) = frequency_index;
        end
    end
    pending = pending(1:pending_count);

    if ~isempty(pending)
        completed = cell(1, numel(pending));
        if parallel_workers == 1
            fprintf(['Solving %d pending RACC-PM frequencies ', ...
                'sequentially; each frequency is checkpointed.\n'], ...
                numel(pending));
            for task_index = 1:numel(pending)
                completed{task_index} = local_run_one( ...
                    HB, HD, desired, para, checkpoint_dir, ...
                    pending(task_index));
            end
        else
            pool = gcp('nocreate');
            if isempty(pool)
                pool = parpool('Processes', parallel_workers);
            end
            workers_to_use = min(parallel_workers, pool.NumWorkers);
            fprintf(['Solving %d pending RACC-PM frequencies with %d ', ...
                'process workers; each frequency is checkpointed.\n'], ...
                numel(pending), workers_to_use);
            parfor (task_index = 1:numel(pending), workers_to_use)
                completed{task_index} = local_run_one( ...
                    HB, HD, desired, para, checkpoint_dir, ...
                    pending(task_index));
            end
        end

        for task_index = 1:numel(completed)
            [filters, status, normalization_multiplier, ...
                normalization_attempt_status] = local_insert_checkpoint( ...
                filters, status, normalization_multiplier, ...
                normalization_attempt_status, completed{task_index});
        end
    end

    solved_mask = arrayfun(@(value) local_status_is_solved(char(value)), ...
        status);
    if any(~solved_mask)
        failed = find(~solved_mask);
        error('RACCPMCheckpoint:UnsolvedFrequencies', ...
            'RACC_PM_Sub failed at frequency indices: %s.', mat2str(failed));
    end

    result = struct();
    result.w = filters;
    result.scale = para.scale;
    result.status = status;
    result.normalization_multiplier = normalization_multiplier;
    result.normalization_attempt_status = normalization_attempt_status;
    result.checkpoint_dir = checkpoint_dir;
    result.implementation_path = implementation_path;
    result.parallel_workers = parallel_workers;
end

function checkpoint = local_run_one( ...
        HB, HD, desired, para, checkpoint_dir, frequency_index)
    rng(20260826 + frequency_index, 'twister');
    solved = RACC_PM_Sub(HB, HD, desired, para, frequency_index);
    checkpoint = local_extract_checkpoint(solved, para, frequency_index);
    local_save_checkpoint(local_checkpoint_file( ...
        checkpoint_dir, frequency_index), checkpoint);
end

function checkpoint = local_extract_checkpoint(solved, para, frequency_index)
    checkpoint = struct();
    checkpoint.frequency_index = frequency_index;
    checkpoint.w = solved.w(:, frequency_index);
    checkpoint.status = solved.status(frequency_index);
    checkpoint.normalization_multiplier = ...
        solved.normalization_multiplier(frequency_index);
    checkpoint.normalization_attempt_status = ...
        solved.normalization_attempt_status{frequency_index};
    checkpoint.epsilon_B = local_pick(para.epsilon.B, frequency_index);
    checkpoint.epsilon_D = local_pick(para.epsilon.D, frequency_index);
    checkpoint.alpha = local_pick(para.alpha, frequency_index);
    checkpoint.rho = para.rho;
    checkpoint.mu = para.mu;
    checkpoint.scale = para.scale;
    checkpoint.completed_at = datetime('now');
end

function [filters, status, multiplier, attempts] = ...
        local_insert_checkpoint(filters, status, multiplier, attempts, ...
        checkpoint)
    index = checkpoint.frequency_index;
    filters(:, index) = checkpoint.w;
    status(index) = string(checkpoint.status);
    multiplier(index) = checkpoint.normalization_multiplier;
    attempts{index} = checkpoint.normalization_attempt_status;
end

function local_validate_checkpoint(checkpoint, para, frequency_index)
    if checkpoint.frequency_index ~= frequency_index || ...
            ~local_equal(checkpoint.epsilon_B, ...
                local_pick(para.epsilon.B, frequency_index)) || ...
            ~local_equal(checkpoint.epsilon_D, ...
                local_pick(para.epsilon.D, frequency_index)) || ...
            ~local_equal(checkpoint.alpha, ...
                local_pick(para.alpha, frequency_index)) || ...
            ~local_equal(checkpoint.rho, para.rho) || ...
            ~local_equal(checkpoint.mu, para.mu) || ...
            ~local_equal(checkpoint.scale, para.scale)
        error('RACCPMCheckpoint:ConfigurationMismatch', ...
            ['Checkpoint %d does not match the current uncertainty, ', ...
             'alpha, rho, mu, or scale configuration.'], frequency_index);
    end
end

function value = local_pick(values, index)
    if isscalar(values)
        value = values;
    else
        value = values(index);
    end
end

function equal = local_equal(a, b)
    tolerance = 1e-12 * max([1, abs(double(a)), abs(double(b))]);
    equal = abs(double(a) - double(b)) <= tolerance;
end

function file = local_checkpoint_file(directory, frequency_index)
    file = fullfile(directory, sprintf( ...
        'racc_pm_frequency_%03d.mat', frequency_index));
end

function local_save_checkpoint(file, checkpoint)
    save(file, 'checkpoint', '-v7');
end

function solved = local_status_is_solved(status)
    solved = strcmp(status, 'Solved') || ...
        strcmp(status, 'Inaccurate/Solved');
end
