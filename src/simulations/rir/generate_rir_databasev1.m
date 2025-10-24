function generate_rir_databasev1(mode, varargin)
% GENERATE_RIR_DATABASE - A flexible function to generate RIR datasets for different experimental scenarios.
%
% Syntax:
%   generate_rir_database('snr', 'snr_range', [10 30], 'num_steps', 21, 'output_dir', 'results/SNR_Experiment');
%   generate_rir_database('temperature', 'temp_range', [-10 40], 'num_steps', 7, 'output_dir', 'results/Temperature_Experiment');
%   generate_rir_database('position', 'num_runs', 50, 'max_offset', 0.05, 'output_dir', 'results/Position_Experiment');
%
% Inputs:
%   mode - string: The experimental scenario. Must be one of:
%          'snr':         Simulates RIRs with varying levels of measurement noise.
%          'temperature': Simulates RIRs under varying room temperatures.
%          'position':    Simulates RIRs with random microphone position perturbations.
%
% Optional Name-Value Pair Arguments:
%   'array_file'    - string: Path to the array geometry .mat file. Default: 'data/array.mat'.
%   'output_dir'    - string: Directory to save the generated data.
%   --- SNR Mode Specific ---
%   'snr_range'     - [1x2]: Min and max SNR in dB. Default: [10 30].
%   'num_steps'     - scalar: Number of steps for SNR/Temperature. Default: 21 for SNR, 7 for Temp.
%   --- Temperature Mode Specific ---
%   'temp_range'    - [1x2]: Min and max temperature in Celsius. Default: [-10 40].
%   --- Position Mode Specific ---
%   'num_runs'      - scalar: Number of random perturbation runs. Default: 50.
%   'max_offset'    - scalar: Max perturbation radius in meters. Default: 0.05.
%   'base_temp'     - scalar: Base temperature for position experiment. Default: 20.

    %% 1. Default Parameters & Input Parsing
    p = inputParser;
    addRequired(p, 'mode', @(x) ismember(x, {'snr', 'temperature', 'position'}));
    addParameter(p, 'array_file', 'data/arrayGeometry/array.mat', @ischar);
    addParameter(p, 'output_dir', '', @ischar);

    % Scenario-specific parameters
    addParameter(p, 'snr_range', [10 30], @isnumeric);
    addParameter(p, 'temp_range', [-10 40], @isnumeric);
    addParameter(p, 'num_steps', 0, @isnumeric); 
    addParameter(p, 'num_runs', 50, @isnumeric);
    addParameter(p, 'max_offset', 0.05, @isnumeric);
    addParameter(p, 'base_temp', 20, @isnumeric);
    parse(p, mode, varargin{:});

    % Assign parsed results to variables
    array_file = p.Results.array_file;
    output_dir = p.Results.output_dir;
    
    % Set default num_steps if not provided
    if p.Results.num_steps == 0
        switch mode
            case 'snr', num_steps = 50;
            case 'temperature', num_steps = 50;
            otherwise, num_steps = 50;
        end
    else
        num_steps = p.Results.num_steps;
    end

    fprintf('====================================================\n');
    fprintf('    RIR Database Generation Tool\n');
    fprintf('====================================================\n');
    fprintf('Selected Mode: %s\n', mode);

    if isempty(output_dir)
        base_dir = fullfile('data', 'SimulateRIR');
        output_dir = fullfile(base_dir, mode);
        % if ~exist(output_dir, 'dir'), mkdir(output_dir); end
        fprintf('Output directory not specified. Using default path:\n  %s\n', output_dir);
    end

    %% 2. General Setup
    freq_params = configure_freq_parameters();
    
    % 使用结构体中的值
    fs = freq_params.fs_native;
    p_resample = freq_params.p_resample;
    q_resample = freq_params.q_resample;

    array = load(array_file);
    beta = 0.3; %
    % This section defines all parameters related to time-frequency transformation,
    % including sampling rates, RIR length, FFT settings, and target frequency bins.

    truncatad_time = 128; % unit: ms, include early reverberation
    rir_len = floor(truncatad_time/1e3 * fs);
    
    s = array.array.s;
    bCtrPtsPositions = array.array.bCtrPtsPositions;
    dCtrPtsPositions = array.array.dCtrPtsPositions;
    bEvalPtsPositions = array.array.bPerPtsPositions; 
    dEvalPtsPositions = array.array.dPerPtsPositions; 

    nCtr = size(bCtrPtsPositions, 1);
    % nSrc = size(s, 1);
   
    %% 3. Mode-Specific Logic & Main Loop
    parpool('local', 8);
    switch mode
        case 'snr'
            snr_vec = linspace(p.Results.snr_range(1), p.Results.snr_range(2), num_steps);
            fprintf('Generating %d RIR sets for SNR range [%.1f, %.1f] dB...\n', ...
                num_steps, snr_vec(1), snr_vec(end));
            
            c = temps2speed(p.Results.base_temp); % Use base temperature

            para.mode = 'snr';
            para.snr_vector = snr_vec;
            para.freq_params = freq_params;
            
            for i = 1:num_steps
                snr_val = snr_vec(i);
                fprintf('  Processing SNR = %.2f dB (Step %d/%d)...\n', snr_val, i, num_steps);
                
                % Generate one clean RIR set first
                [IR_BZ, IR_DZ] = generate_single_rir_set(c, fs, rir_len, s, ...
                    bCtrPtsPositions, dCtrPtsPositions, bEvalPtsPositions,...
                    dEvalPtsPositions, roomSize, beta);
                
                % Add noise
                IR_BZ.ctrl = awgn(IR_BZ.ctrl, snr_val, 'measured');
                IR_DZ.ctrl = awgn(IR_DZ.ctrl, snr_val, 'measured');
                IR_BZ.eval = awgn(IR_BZ.eval, snr_val, 'measured'); 
                IR_DZ.eval = awgn(IR_DZ.eval, snr_val, 'measured'); 
                
                % Resample after adding noise
                IR_BZ.ctrl = resample(IR_BZ.ctrl, p_resample, q_resample, 'Dimension', 3);
                IR_DZ.ctrl = resample(IR_DZ.ctrl, p_resample, q_resample, 'Dimension', 3);
                IR_BZ.eval = resample(IR_BZ.eval, p_resample, q_resample, 'Dimension', 3);
                IR_DZ.eval = resample(IR_DZ.eval, p_resample, q_resample, 'Dimension', 3);

                % compute all ATFs
                ATF_BZ.ctrl = compute_atf(IR_BZ.ctrl, freq_params);
                ATF_BZ.eval = compute_atf(IR_BZ.eval, freq_params);
                ATF_DZ.ctrl = compute_atf(IR_DZ.ctrl, freq_params);
                ATF_DZ.eval = compute_atf(IR_DZ.eval, freq_params);

                % convert to single datatype
                IR_BZ = structfun(@single, IR_BZ, 'UniformOutput', false);
                IR_DZ = structfun(@single, IR_DZ, 'UniformOutput', false);
                ATF_BZ = structfun(@single, ATF_BZ, 'UniformOutput', false);
                ATF_DZ = structfun(@single, ATF_DZ, 'UniformOutput', false);
                
                % Save data
                output_filename = fullfile(output_dir, sprintf('Data_SNR-%.2f.mat', snr_val));
                save(output_filename, 'IR_BZ', 'IR_DZ', 'ATF_BZ', 'ATF_DZ','-v7.3');
            end

        case 'temperature'
            temp_vec = linspace(p.Results.temp_range(1), p.Results.temp_range(2), num_steps);
            fprintf('Generating %d RIR sets for Temperature range [%.1f, %.1f] C...\n', ...
                num_steps, temp_vec(1), temp_vec(end));
            
            para.mode = 'temperature';
            para.temperature_vector_celsius = temp_vec;
            para.freq_params = freq_params;

            for i = 1:num_steps
                temp_val = temp_vec(i);
                c = temps2speed(temp_val);
                fprintf(' Processing Temp = %.2f C (Sound Speed = %.2f m/s)...\n', temp_val, c);
                
                [IR_BZ, IR_DZ] = generate_single_rir_set(c, fs, rir_len, ...
                    s, bCtrPtsPositions, dCtrPtsPositions, bEvalPtsPositions, ...
                    dEvalPtsPositions, roomSize, beta);
                
                % Resample
                IR_BZ.ctrl = resample(IR_BZ.ctrl, p_resample, q_resample, 'Dimension', 3);
                IR_DZ.ctrl = resample(IR_DZ.ctrl, p_resample, q_resample, 'Dimension', 3);
                IR_BZ.eval = resample(IR_BZ.eval, p_resample, q_resample, 'Dimension', 3);
                IR_DZ.eval = resample(IR_DZ.eval, p_resample, q_resample, 'Dimension', 3);

                % compute all ATFs
                ATF_BZ.ctrl = compute_atf(IR_BZ.ctrl, freq_params);
                ATF_BZ.eval = compute_atf(IR_BZ.eval, freq_params);
                ATF_DZ.ctrl = compute_atf(IR_DZ.ctrl, freq_params);
                ATF_DZ.eval = compute_atf(IR_DZ.eval, freq_params);

                % convert to single datatype
                IR_BZ = structfun(@single, IR_BZ, 'UniformOutput', false);
                IR_DZ = structfun(@single, IR_DZ, 'UniformOutput', false);
                ATF_BZ = structfun(@single, ATF_BZ, 'UniformOutput', false);
                ATF_DZ = structfun(@single, ATF_DZ, 'UniformOutput', false);
                
                % Save data
                output_filename = fullfile(output_dir, sprintf('Data_T-%.2f.mat', temp_val));
                para.temperature = temp_val;
                para.sound_speed = c;
                para.mode = 'temperature';
                save(output_filename, 'IR_BZ', 'IR_DZ', 'ATF_BZ', 'ATF_DZ','-v7.3');
            end

        case 'position'
            num_steps = p.Results.num_runs; % Re-interpret num_runs as number of levels
            max_offset = p.Results.max_offset;

            % Create a vector of perturbation radii from 0 to max_offset
            perturb_radii = linspace(0, max_offset, num_steps);

            % Pre-allocate metadata storage
            offsets_metadata = zeros(num_steps, nCtr, 3);

            fprintf('Generating %d RIR/ATF sets for stratified position perturbations...\n', num_steps);
            fprintf('Radius levels from %.3f m to %.3f m.\n', perturb_radii(1), perturb_radii(end));
           
            c = temps2speed(p.Results.base_temp); % Use base temperature

            % --- Prepare common metadata for para.mat ---
            para.mode = 'position';
            para.perturbation_radii_m = perturb_radii;
            para.freq_params = freq_params;
            
            for i = 1:num_steps
                current_radius = perturb_radii(i);
                fprintf('  Processing Level %d/%d: Radius = %.4f m...\n', i, num_steps, current_radius);
                
                % Generate a random offset with a FIXED radius and RANDOM direction
                [bCtrPts_perturbed, dCtrPts_perturbed, bEvalPtsPositions,dEvalPtsPositions, run_offsets] = ...
                    perturb_all_positions_fixed_radius(bCtrPtsPositions, dCtrPtsPositions, ...
                                                       bEvalPtsPositions, dEvalPtsPositions, current_radius);

                % Store the generated offset vector for metadata
                offsets_metadata(i, :, :) = run_offsets;
                
                % Generate RIR with perturbed positions
                [IR_BZ, IR_DZ] = generate_single_rir_set(c, fs, rir_len, s, ...
                    bCtrPts_perturbed, dCtrPts_perturbed, bEvalPtsPositions, ...
                    dEvalPtsPositions, roomSize, beta);
                
                % Resample
                IR_BZ.ctrl = resample(IR_BZ.ctrl, p_resample, q_resample, 'Dimension', 3);
                IR_DZ.ctrl = resample(IR_DZ.ctrl, p_resample, q_resample, 'Dimension', 3);
                IR_BZ.eval = resample(IR_BZ.eval, p_resample, q_resample, 'Dimension', 3);
                IR_DZ.eval = resample(IR_DZ.eval, p_resample, q_resample, 'Dimension', 3);

                % compute all ATFs
                ATF_BZ.ctrl = compute_atf(IR_BZ.ctrl, freq_params);
                ATF_BZ.eval = compute_atf(IR_BZ.eval, freq_params);
                ATF_DZ.ctrl = compute_atf(IR_DZ.ctrl, freq_params);
                ATF_DZ.eval = compute_atf(IR_DZ.eval, freq_params);

                % convert to single datatype
                IR_BZ = structfun(@single, IR_BZ, 'UniformOutput', false);
                IR_DZ = structfun(@single, IR_DZ, 'UniformOutput', false);
                ATF_BZ = structfun(@single, ATF_BZ, 'UniformOutput', false);
                ATF_DZ = structfun(@single, ATF_DZ, 'UniformOutput', false);
                
                % Save data
                output_filename = fullfile(output_dir, sprintf('Data_Pos-radius-%.4f.mat', current_radius));                
                save(output_filename, 'IR_BZ', 'IR_DZ', 'ATF_BZ', 'ATF_DZ','-v7.3');
            end
            % Save the metadata file with all offsets
            offset_meta_file = fullfile(output_dir, 'Position_Offsets_Metadata.mat');
            save(offset_meta_file, 'offsets_metadata', 'max_offset');
            fprintf('Detailed offset vectors saved to: %s\n', offset_meta_file);
    end

    para_filename = fullfile(output_dir, 'para.mat');
    save(para_filename, 'para');
    fprintf('\nCommon metadata saved to: %s\n', para_filename);

    delete(gcp('nocreate'));
    fprintf('\nGeneration complete for mode: %s\n', mode);
    fprintf('====================================================\n');
end


%% Other Functions
function [IR_BZ, IR_DZ] = generate_single_rir_set(c, fs, rir_len, s, ...
    b_ctrl_mics, d_ctrl_mics, b_eval_mics, d_eval_mics, room, beta)
    % Generates a full set of RIRs for both control and evaluation points
    nCtr = size(b_ctrl_mics, 1);
    nSrc = size(s, 1);
    
    IR_BZ_ctrl = zeros(nCtr, nSrc, rir_len);
    IR_DZ_ctrl = zeros(nCtr, nSrc, rir_len);
    IR_BZ_eval = zeros(nCtr, nSrc, rir_len);
    IR_DZ_eval = zeros(nCtr, nSrc, rir_len);
    
    for i = 1:nSrc
        si = s(i,:);
        parfor j = 1:nCtr
            % Control Points
            IR_BZ_ctrl(j,i,:) = rir_generator(c, fs, b_ctrl_mics(j,:), si, room, beta, rir_len);
            IR_DZ_ctrl(j,i,:) = rir_generator(c, fs, d_ctrl_mics(j,:), si, room, beta, rir_len);
            % Evaluation Points
            IR_BZ_eval(j,i,:) = rir_generator(c, fs, b_eval_mics(j,:), si, room, beta, rir_len);
            IR_DZ_eval(j,i,:) = rir_generator(c, fs, d_eval_mics(j,:), si, room, beta, rir_len);
        end
    end
    
    % Organize into structures
    IR_BZ.ctrl = IR_BZ_ctrl;
    IR_BZ.eval = IR_BZ_eval;
    IR_DZ.ctrl = IR_DZ_ctrl;
    IR_DZ.eval = IR_DZ_eval;
end


function [b_ctrl_p, d_ctrl_p, b_eval_p, d_eval_p, offsets] = ...
    perturb_all_positions_fixed_radius(b_ctrl, d_ctrl, b_eval, d_eval, target_r)

% Perturbs all microphone positions with a random direction but a FIXED radius.
    nCtr = size(b_ctrl, 1);
    offsets = zeros(nCtr, 3);
    
    for i = 1:nCtr
        % Radius is now the fixed target radius
        r = target_r;
        
        % Direction generation remains the same (random on a sphere)
        theta = 2 * pi * rand();
        phi = asin(2*rand() - 1);
        
        % Handle the special case of r=0 to avoid generating NaN from phi
        if r == 0
            dx = 0; dy = 0; dz = 0;
        else
            dx = r * cos(phi) * cos(theta);
            dy = r * cos(phi) * sin(theta);
            dz = r * sin(phi);
        end
        
        offsets(i, :) = [dx, dy, dz];
    end
    
    % Apply the same offset to control and evaluation points with the same index
    b_ctrl_p = b_ctrl + offsets;
    d_ctrl_p = d_ctrl + offsets;
    b_eval_p = b_eval + offsets;
    d_eval_p = d_eval + offsets;
end