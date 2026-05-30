function generate_rir_ATF_database(mode, params, varargin)
% GENERATE_RIR_DATABASE - A flexible function to generate RIR and ATF datasets. (Version 4.0 Refactored)
%
% Syntax (updated):
%   generate_rir_database('snr', 'snr_range', [10 30], 'num_steps', 21);
%   generate_rir_database('temperature', 'temp_range', [-10 40], 'num_steps', 50);
%   generate_rir_database('position', 'num_runs', 51, 'max_offset', 0.05);
%
% mode             : 'snr' | 'temperature' | 'position'

%% 1. Default Parameters & Input Parsing
p = inputParser;
addRequired(p, 'mode', @(x) ismember(x, {'snr', 'temperature', 'position'}));
addParameter(p, 'array_file', 'data/arrayGeometry/array_layout.mat', @ischar);
addParameter(p, 'output_dir', '', @ischar);
addParameter(p, 'snr_range', [10 30], @isnumeric);
addParameter(p, 'temp_range', [20 25], @isnumeric);
addParameter(p, 'num_steps', 51, @isnumeric); % Default for snr/temp
addParameter(p, 'num_runs', 51, @isnumeric);  % Default for position (0 to N-1 steps)
addParameter(p, 'max_offset', 0.2, @isnumeric); % the perturbation range:0~0.2m
addParameter(p, 'base_temp', 20, @isnumeric);
% IMPORTANT:
parse(p, mode, varargin{:});
% --- Assign parsed results to variables ---
array_file = p.Results.array_file;
output_dir = p.Results.output_dir;

% --- Dynamically set output directory if not specified ---
if isempty(output_dir)
    base_dir = fullfile('data', 'SimulateRIR');
    output_dir = fullfile(base_dir, mode);
end
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% 2. General Setup
% --- Geometry Data ---
array_data_wrapper = load(array_file);
array_geom = array_data_wrapper.roomArray; % Unpack the nested struct

% --- Define Simulation Constants ---
beta = params.beta;
fs = params.fs;
rir_len = params.rir_len;


% --- Extract Geometry ---
s = array_geom.s;
BZ.BZ_ctrl = array_geom.BZ_ctrl;
DZ.DZ_ctrl = array_geom.DZ_ctrl;
BZ.BZ_eval = array_geom.BZ_eval;
DZ.DZ_eval = array_geom.DZ_eval;
roomSize = array_geom.roomSize;
nCtr = size(BZ.BZ_ctrl, 1);

%% 3. Mode-Specific Logic & Main Loop
switch mode
    case 'snr'
        num_steps = p.Results.num_steps;
        snr_vec = linspace(p.Results.snr_range(1), p.Results.snr_range(2), num_steps);
        c = temp2speed(p.Results.base_temp);
        para.mode = 'snr';
        para.snr_vector = snr_vec;
        para.freq_params = params;
        [base_IR_BZ, base_IR_DZ] = generate_full_rir_set(c, fs, rir_len, s, BZ, DZ, roomSize, beta, params);

        for i = 1:num_steps
            snr_val = snr_vec(i);
            % Add noise to both control and evaluation points
            IR_BZ.ctrl = awgn(base_IR_BZ.ctrl, snr_val, 'measured');
            IR_BZ.eval = awgn(base_IR_BZ.eval, snr_val, 'measured');
            IR_DZ.ctrl = awgn(base_IR_DZ.ctrl, snr_val, 'measured');
            IR_DZ.eval = awgn(base_IR_DZ.eval, snr_val, 'measured');
            output_filename = fullfile(output_dir, sprintf('Data_SNR-%.2f.mat', snr_val));
            process_and_save_data(IR_BZ, IR_DZ, params, output_filename);
        end

    case 'temperature'
        num_steps = p.Results.num_steps;
        temp_vec = linspace(p.Results.temp_range(1), p.Results.temp_range(2), num_steps);
        para.mode = 'temperature';
        para.temperature_vector_celsius = temp_vec;
        para.freq_params = params;
        for i = 1:num_steps
            temp_val = temp_vec(i);
            c = temp2speed(temp_val);
            [IR_BZ, IR_DZ] = generate_full_rir_set(c, fs, rir_len, s, BZ, DZ, roomSize, beta, params);
            output_filename = fullfile(output_dir, sprintf('Data_T-%.2f.mat', temp_val));
            process_and_save_data(IR_BZ, IR_DZ, params, output_filename);
        end

    case 'position'
        num_levels = p.Results.num_runs;
        max_offset = p.Results.max_offset;
        perturb_radii = linspace(0, max_offset, num_levels);
        offsets_metadata = zeros(num_levels, nCtr, 3);
        c = temp2speed(p.Results.base_temp);
        para.mode = 'position';
        para.perturbation_radii_m = perturb_radii;
        para.freq_params = params;
        for i = 1:num_levels
            current_radius = perturb_radii(i);
            [BZ_p, DZ_p, run_offsets] = perturb_positions_fixed_radius(BZ, DZ,current_radius);
            offsets_metadata(i, :, :) = run_offsets;
            [IR_BZ, IR_DZ] = generate_full_rir_set(c, fs, rir_len, s, BZ_p, DZ_p, roomSize, beta, params);
            output_filename = fullfile(output_dir, sprintf('Data_Pos-%.4f.mat', current_radius));
            process_and_save_data(IR_BZ, IR_DZ, params,output_filename);
        end
        offset_meta_file = fullfile(output_dir, 'Position_Offsets_Metadata.mat');
        save(offset_meta_file, 'offsets_metadata');
end
% save shared metadata
para_filename = fullfile(output_dir, 'para.mat');
save(para_filename, 'para');
end


%% Functions
function [IR_BZ, IR_DZ] = generate_full_rir_set(c, fs, rir_len, s, BZ, DZ, room, beta, params)
    BZ_ctrl = BZ.BZ_ctrl; DZ_ctrl = DZ.DZ_ctrl;
    BZ_eval = BZ.BZ_eval; DZ_eval = DZ.DZ_eval;
        
    nCtr = size(BZ_ctrl, 1);
    nSrc = size(s, 1);
    w_tail = params.w_tail;
    fade_samples = params.fade_samples;

    IR_BZ_ctrl = zeros(nCtr, nSrc, rir_len);
    IR_DZ_ctrl = zeros(nCtr, nSrc, rir_len);
    IR_BZ_eval = zeros(nCtr, nSrc, rir_len);
    IR_DZ_eval = zeros(nCtr, nSrc, rir_len);

    parfor i = 1:nSrc
        si = s(i,:);
        RIR_tmp = rir_generator(c, fs, BZ_ctrl, si, room, beta, rir_len);
        RIR_tmp(:,end-fade_samples+1:end) = RIR_tmp(:, end-fade_samples+1:end) .* w_tail';
        RIR_tmp = single(RIR_tmp);
        RIR_tmp = reshape(RIR_tmp, nCtr, 1, rir_len);
        IR_BZ_ctrl(:,i,:) = RIR_tmp;

        RIR_tmp = rir_generator(c, fs, BZ_eval, si, room, beta, rir_len);
        RIR_tmp(:,end-fade_samples+1:end) = RIR_tmp(:, end-fade_samples+1:end) .* w_tail';
        RIR_tmp = single(RIR_tmp);
        RIR_tmp = reshape(RIR_tmp, nCtr, 1, rir_len);
        IR_BZ_eval(:,i,:) = RIR_tmp;

        RIR_tmp = rir_generator(c, fs, DZ_ctrl, si, room, beta, rir_len);
        RIR_tmp(:,end-fade_samples+1:end) = RIR_tmp(:, end-fade_samples+1:end) .* w_tail';
        RIR_tmp = single(RIR_tmp);
        RIR_tmp = reshape(RIR_tmp, nCtr, 1, rir_len);
        IR_DZ_ctrl(:,i,:) = RIR_tmp;

        RIR_tmp = rir_generator(c, fs, DZ_eval, si, room, beta, rir_len);
        RIR_tmp(:,end-fade_samples+1:end) = RIR_tmp(:, end-fade_samples+1:end) .* w_tail';
        RIR_tmp = single(RIR_tmp);
        RIR_tmp = reshape(RIR_tmp, nCtr, 1, rir_len);
        IR_DZ_eval(:,i,:) = RIR_tmp;
    end

    IR_BZ.ctrl = IR_BZ_ctrl;
    IR_BZ.eval = IR_BZ_eval;
    IR_DZ.ctrl = IR_DZ_ctrl;
    IR_DZ.eval = IR_DZ_eval;
end


function [BZ_p, DZ_p, offsets] = perturb_positions_fixed_radius(BZ, DZ, target_r)

    BZ_ctrl = BZ.BZ_ctrl; DZ_ctrl = DZ.DZ_ctrl;
    BZ_eval = BZ.BZ_eval; DZ_eval = DZ.DZ_eval;
        
    nCtr = size(BZ_ctrl, 1);
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

    BZ_p.BZ_ctrl = BZ_ctrl + offsets;
    DZ_p.DZ_ctrl = DZ_ctrl + offsets;
    BZ_p.BZ_eval = BZ_eval + offsets;
    DZ_p.DZ_eval = DZ_eval + offsets;
end


function process_and_save_data(IR_BZ, IR_DZ, freq_params, output_filename)
% This function handles ATF computation and saving.
    ATF_BZ.ctrl = compute_atf(IR_BZ.ctrl, freq_params);
    ATF_BZ.eval = compute_atf(IR_BZ.eval, freq_params);
    ATF_DZ.ctrl = compute_atf(IR_DZ.ctrl, freq_params);
    ATF_DZ.eval = compute_atf(IR_DZ.eval, freq_params);
    % Save data
    save(output_filename, 'IR_BZ', 'IR_DZ', 'ATF_BZ', 'ATF_DZ','-v7.3');
end
