function generate_rir_database(mode, varargin)
% GENERATE_RIR_DATABASE - A flexible function to generate RIR and ATF datasets. (Version 4.0 Refactored)
%
% Syntax:
%   generate_rir_database('snr', 'snr_range', [10 30], 'num_steps', 21);
%   generate_rir_database('temperature', 'temp_range', [-10 40], 'num_steps', 50);
%   generate_rir_database('position', 'num_runs', 51, 'max_offset', 0.05);

    %% 1. Default Parameters & Input Parsing
    p = inputParser;
    addRequired(p, 'mode', @(x) ismember(x, {'snr', 'temperature', 'position'}));
    addParameter(p, 'array_file', 'data/arrayGeometry/array.mat', @ischar);
    addParameter(p, 'output_dir', '', @ischar);
    addParameter(p, 'snr_range', [10 30], @isnumeric);
    addParameter(p, 'temp_range', [-10 40], @isnumeric);
    addParameter(p, 'num_steps', 51, @isnumeric); % Default for snr/temp
    addParameter(p, 'num_runs', 51, @isnumeric);  % Default for position (0 to N-1 steps)
    addParameter(p, 'max_offset', 0.05, @isnumeric);
    addParameter(p, 'base_temp', 20, @isnumeric);
    parse(p, mode, varargin{:});

    % --- Assign parsed results to variables ---
    array_file = p.Results.array_file;
    output_dir = p.Results.output_dir;
    
    fprintf('====================================================\n');
    fprintf('   RIR & ATF Database Generation Tool Is Working!! \n');
    fprintf('====================================================\n');
    fprintf('Selected Mode: %s\n', mode);

    % --- Dynamically set output directory if not specified ---
    if isempty(output_dir)
        base_dir = fullfile('data', 'SimulateRIR');
        output_dir = fullfile(base_dir, mode);
        fprintf('Output directory not specified. Using default path:\n  %s\n', output_dir);
    end
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    
    %% 2. General Setup
    % --- Load Frequency and Geometry Data ---
    freq_params = configure_freq_parameters();
    array_data_wrapper = load(array_file);
    array_geom = array_data_wrapper.array; % Unpack the nested struct

    % --- Define Simulation Constants ---
    beta = 0.3;
    fs = freq_params.fs;
    rir_len = freq_params.rir_len_native;
    p_resample = freq_params.p_resample;
    q_resample = freq_params.q_resample;
    
    % --- Extract Geometry ---
    s = array_geom.s;
    bCtrPtsPositions = array_geom.bCtrPtsPositions;
    dCtrPtsPositions = array_geom.dCtrPtsPositions;
    bEvalPtsPositions = array_geom.bPerPtsPositions; 
    dEvalPtsPositions = array_geom.dPerPtsPositions; 
    roomSize = array_geom.roomSize;
    nCtr = size(bCtrPtsPositions, 1);
   
    %% 3. Mode-Specific Logic & Main Loop
    parpool('local', 8);
    
    switch mode
        case 'snr'
            num_steps = p.Results.num_steps;
            snr_vec = linspace(p.Results.snr_range(1), p.Results.snr_range(2), num_steps);
            fprintf('Generating %d RIR/ATF sets for SNR range [%.1f, %.1f] dB...\n', num_steps, snr_vec(1), snr_vec(end));
            
            c = temp2speed(p.Results.base_temp);
            para.mode = 'snr';
            para.snr_vector = snr_vec;
            para.freq_params = freq_params;

            [IR_BZ, IR_DZ] = generate_full_rir_set(c, fs, rir_len, s, ...
                    bCtrPtsPositions, dCtrPtsPositions, bEvalPtsPositions, dEvalPtsPositions, roomSize, beta);
            
            for i = 1:num_steps
                snr_val = snr_vec(i);
                fprintf('  Processing SNR = %.2f dB (Step %d/%d)...\n', snr_val, i, num_steps);
                     
                % Add noise to both control and evaluation points
                IR_BZ.ctrl = awgn(IR_BZ.ctrl, snr_val, 'measured');
                IR_BZ.eval = awgn(IR_BZ.eval, snr_val, 'measured'); 
                IR_DZ.ctrl = awgn(IR_DZ.ctrl, snr_val, 'measured');
                IR_DZ.eval = awgn(IR_DZ.eval, snr_val, 'measured'); 
                
                output_filename = fullfile(output_dir, sprintf('Data_SNR-%.2f.mat', snr_val));
                process_and_save_data(IR_BZ, IR_DZ, freq_params, p_resample, q_resample, output_filename);
            end

        case 'temperature'
            num_steps = p.Results.num_steps;
            temp_vec = linspace(p.Results.temp_range(1), p.Results.temp_range(2), num_steps);
            fprintf('Generating %d RIR/ATF sets for Temperature range [%.1f, %.1f] C...\n', num_steps, temp_vec(1), temp_vec(end));
            
            para.mode = 'temperature';
            para.temperature_vector_celsius = temp_vec;
            para.freq_params = freq_params;

            for i = 1:num_steps
                temp_val = temp_vec(i);
                c = temp2speed(temp_val);
                fprintf('  Processing Temp = %.2f C (Sound Speed = %.2f m/s)...\n', temp_val, c);
                
                [IR_BZ, IR_DZ] = generate_full_rir_set(c, fs, rir_len, s, ...
                    bCtrPtsPositions, dCtrPtsPositions, bEvalPtsPositions, dEvalPtsPositions, roomSize, beta);
                
                output_filename = fullfile(output_dir, sprintf('Data_T-%.2f.mat', temp_val));
                process_and_save_data(IR_BZ, IR_DZ, freq_params, p_resample, q_resample, output_filename);
            end

        case 'position'
            num_levels = p.Results.num_runs;
            max_offset = p.Results.max_offset;
            perturb_radii = linspace(0, max_offset, num_levels);
            offsets_metadata = zeros(num_levels, nCtr, 3);

            fprintf('Generating %d RIR/ATF sets for stratified position perturbations...\n', num_levels);
            fprintf('Radius levels from %.3f m to %.3f m.\n', perturb_radii(1), perturb_radii(end));
           
            c = temp2speed(p.Results.base_temp);
            para.mode = 'position';
            para.perturbation_radii_m = perturb_radii;
            para.freq_params = freq_params;
            
            for i = 1:num_levels
                current_radius = perturb_radii(i);
                fprintf('  Processing Level %d/%d: Radius = %.4f m...\n', i, num_levels, current_radius);
                
                [b_ctrl_p, d_ctrl_p, b_eval_p, d_eval_p, run_offsets] = ...
                    perturb_all_positions_fixed_radius(bCtrPtsPositions, dCtrPtsPositions, ...
                                                       bEvalPtsPositions, dEvalPtsPositions, current_radius);
                offsets_metadata(i, :, :) = run_offsets;
                
                [IR_BZ, IR_DZ] = generate_full_rir_set(c, fs, rir_len, s, ...
                    b_ctrl_p, d_ctrl_p, b_eval_p, d_eval_p, roomSize, beta);
                
                output_filename = fullfile(output_dir, sprintf('Data_Pos-%.4f.mat', current_radius));
                process_and_save_data(IR_BZ, IR_DZ, freq_params, p_resample, q_resample, output_filename);
            end
            
            offset_meta_file = fullfile(output_dir, 'Position_Offsets_Metadata.mat');
            save(offset_meta_file, 'offsets_metadata');
            fprintf('Detailed offset vectors saved to: %s\n', offset_meta_file);
    end

    para_filename = fullfile(output_dir, 'para.mat');
    save(para_filename, 'para');
    fprintf('\nCommon metadata saved to: %s\n', para_filename);

    delete(gcp('nocreate'));
    fprintf('\nGeneration complete for mode: %s\n', mode);
    fprintf('====================================================\n');
end


%% Helper Functions
function [IR_BZ, IR_DZ] = generate_full_rir_set(c, fs, rir_len, s, ...
    b_ctrl_mics, d_ctrl_mics, b_eval_mics, d_eval_mics, room, beta)

    nCtr = size(b_ctrl_mics, 1);
    nSrc = size(s, 1);
    
    IR_BZ_ctrl = zeros(nCtr, nSrc, rir_len);
    IR_DZ_ctrl = zeros(nCtr, nSrc, rir_len);
    IR_BZ_eval = zeros(nCtr, nSrc, rir_len);
    IR_DZ_eval = zeros(nCtr, nSrc, rir_len);
    
    for i = 1:nSrc
        si = s(i,:);
        parfor j = 1:nCtr
            IR_BZ_ctrl(j,i,:) = rir_generator(c, fs, b_ctrl_mics(j,:), si, room, beta, rir_len);
            IR_DZ_ctrl(j,i,:) = rir_generator(c, fs, d_ctrl_mics(j,:), si, room, beta, rir_len);
            IR_BZ_eval(j,i,:) = rir_generator(c, fs, b_eval_mics(j,:), si, room, beta, rir_len);
            IR_DZ_eval(j,i,:) = rir_generator(c, fs, d_eval_mics(j,:), si, room, beta, rir_len);
        end
    end
    
    IR_BZ.ctrl = IR_BZ_ctrl;
    IR_BZ.eval = IR_BZ_eval;
    IR_DZ.ctrl = IR_DZ_ctrl;
    IR_DZ.eval = IR_DZ_eval;
end

function [b_ctrl_p, d_ctrl_p, b_eval_p, d_eval_p, offsets] = perturb_all_positions_fixed_radius(b_ctrl, d_ctrl, b_eval, d_eval, target_r)
    nCtr = size(b_ctrl, 1);
    offsets = zeros(nCtr, 3);
    
    for i = 1:nCtr
        r = target_r;
        theta = 2 * pi * rand();
        phi = asin(2*rand() - 1);
        
        if r == 0
            dx = 0; dy = 0; dz = 0;
        else
            dx = r * cos(phi) * cos(theta);
            dy = r * cos(phi) * sin(theta);
            dz = r * sin(phi);
        end
        
        offsets(i, :) = [dx, dy, dz];
    end
    
    b_ctrl_p = b_ctrl + offsets;
    d_ctrl_p = d_ctrl + offsets;
    b_eval_p = b_eval + offsets;
    d_eval_p = d_eval + offsets;
end

function process_and_save_data(IR_BZ, IR_DZ, freq_params, p_resample, q_resample, output_filename)
    % This function handles resampling, ATF computation, data type conversion, and saving.
    
    % 1. Resample
    IR_BZ.ctrl = resample(IR_BZ.ctrl, p_resample, q_resample, 'Dimension', 3);
    IR_BZ.eval = resample(IR_BZ.eval, p_resample, q_resample, 'Dimension', 3);
    IR_DZ.ctrl = resample(IR_DZ.ctrl, p_resample, q_resample, 'Dimension', 3);
    IR_DZ.eval = resample(IR_DZ.eval, p_resample, q_resample, 'Dimension', 3);

    % 2. Compute all ATFs
    ATF_BZ.ctrl = compute_atf(IR_BZ.ctrl, freq_params);
    ATF_BZ.eval = compute_atf(IR_BZ.eval, freq_params);
    ATF_DZ.ctrl = compute_atf(IR_DZ.ctrl, freq_params);
    ATF_DZ.eval = compute_atf(IR_DZ.eval, freq_params);

    % 3. Convert to single datatype
    IR_BZ = structfun(@single, IR_BZ, 'UniformOutput', false);
    IR_DZ = structfun(@single, IR_DZ, 'UniformOutput', false);
    ATF_BZ = structfun(@single, ATF_BZ, 'UniformOutput', false);
    ATF_DZ = structfun(@single, ATF_DZ, 'UniformOutput', false);
    
    % 4. Save data
    save(output_filename, 'IR_BZ', 'IR_DZ', 'ATF_BZ', 'ATF_DZ','-v7.3');
end