%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2025, Lei Zhou
% All rights reserved.
%
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree.
%
% @author: Lei Zhou (zhouleicqupt2016@outlook.com)
% @date: 2025/11/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MAIN SCRIPT - Robust PSZ Algorithm Simulation
% This script serves as the main entry point for the entire project.
% It sets up the environment, runs simulations,executes the core algorithm,
% and visualizes the results.
clear;
clc;
close all;

% Add source code folders to the path
% This ensures MATLAB can find your generate_rir_database function
addpath(genpath('src/'));

rng(2025); % For reproducibility

%% MAIN SCRIPT - Quick Debug & Verification (Refactored)
% This script systematically generate rir databases of all modes 

fprintf('--- Starting Generator for All Generation Modes ---\n\n');

% --- Define all modes in a cell array ---
modes_to_test = {'snr', 'temperature', 'position'};

% --- Loop through each mode and run the generator ---
for i = 1:length(modes_to_test)% <------这里可以更改生成模式
    
    current_mode = modes_to_test{i};
    
    fprintf('--- Testing Mode: %s ---\n', upper(current_mode));
    
    try
        % Define specific small parameters for each mode for a quick test
        switch current_mode
            case 'snr'
                generate_rir_database(current_mode);
            case 'temperature'
                generate_rir_database(current_mode);
            case 'position'
                generate_rir_database(current_mode);
            otherwise
                % Default call if no specific parameters are needed
                generate_rir_database(current_mode);
        end
        
        % --- Success Message ---
        fprintf('SUCCESS: "%s" mode generation is completed.\n', current_mode);
        fprintf('Please check the "data/SimulateRIR/%s/" folder.\n\n', current_mode);
        
    catch ME
        % --- Error Message ---
        fprintf('ERROR in "%s" mode test: %s\n', current_mode, ME.message);
        % Optional: rethrow(ME); % Uncomment if you want the script to stop on error
    end
    
end

fprintf('--- Dataset Generation is Finished ---\n');

% % =========================================================================
% % 1. SETUP ENVIRONMENT
% % =========================================================================
% clear;             % Clear workspace
% clc;               % Clear command window
% close all;         % Close all figures
% 
% fprintf('Setting up environment...\n');
% % Add all source code folders and subfolders to the MATLAB path
% addpath(genpath('src/'));
% % Add the RIR generator library to the path
% addpath('lib/rir_generator/');
% fprintf('Environment setup complete.\n\n');
% 
% % =========================================================================
% % 2. DEFINE PARAMETERS & FILE PATHS
% % =========================================================================
% 
% fprintf('Defining parameters...\n');
% % --- Input Data ---
% paths.array_data = 'data/arrayGeometry/array.mat'; % 假设 array.mat 在这里
% paths.para_data  = 'data/para.mat';
% 
% % --- Simulation & Results ---
% paths.results_dir = 'results/SoundFiled2D/'; % 统一存放结果的目录
% paths.rir_data    = fullfile(paths.results_dir, 'gridRIR_data.mat');
% paths.figure_dir  = fullfile(paths.results_dir, 'figures/');
% 
% % --- Simulation Parameters ---
% sim_params.fs = 16000;
% sim_params.f_target = 1000; % Target frequency for analysis
% 
% % --- Visualization Parameters ---
% vis_params.dynamic_range = 40;
% vis_params.colormap = 'viridis';
% 
% % Create results directories if they don't exist
% if ~exist(paths.results_dir, 'dir'), mkdir(paths.results_dir); end
% if ~exist(paths.figure_dir, 'dir'), mkdir(paths.figure_dir); end
% 
% fprintf('Parameters defined.\n\n');


% =========================================================================
% 3. DATA GENERATION (Simulate RIRs)
% =========================================================================
% Check if RIR data already exists to save time
% if ~exist(paths.rir_data, 'file')
%     fprintf('RIR data not found. Generating now...\n');
% 
%     % 调用您的RIR生成函数 (假设它在 src/simulations/ 中)
%     % 注意：您可能需要根据 generateAndSaveListenRIR 的参数来调整这里的调用
%     generateAndSaveListenRIR(paths.array_data, 'save_filename', paths.rir_data);
% 
%     fprintf('RIR data generation complete.\n\n');
% else
%     fprintf('RIR data found at "%s". Skipping generation.\n\n', paths.rir_data);
% end


% =========================================================================
% 4. ALGORITHM EXECUTION (Calculate Filters)
% =========================================================================
% fprintf('Executing core algorithm to calculate filters...\n');

% --- 在这里，您需要调用您的核心算法 ---
% 这是一个示例，您需要替换为您自己的算法函数
% 假设我们暂时不使用滤波器，先看自然声场
% N_speakers = 64; % 您需要从 array.mat 中获取真实的扬声器数量
% filter_length = 128;
% filters_w = zeros(N_speakers, filter_length);
% filters_w(:, 1) = 1; % 单位脉冲，代表“无滤波器”

% load(paths.rir_data);
% filters_w = your_psz_algorithm(RIR, ...); % <--- 替换为您的真实算法

% fprintf('Filter calculation complete.\n\n');


% =========================================================================
% 5. VISUALIZE RESULTS
% =========================================================================
% fprintf('Visualizing the sound field...\n');

% 调用您的绘图函数 (假设它在 src/visualization/ 中)
% plot_soundfield_pressure(paths.rir_data, filters_w, sim_params.f_target, ...
%     'fs', sim_params.fs, ...
%     'DynamicRange', vis_params.dynamic_range, ...
%     'Colormap', vis_params.colormap, ...
%     'SpeakerMasking', true, ...
%     'SpeakerDataFile', paths.array_data);
% 
% fprintf('Visualization complete. Check the generated figure.\n');
% fprintf('\n--- Script Finished ---\n');


% % grid_RIR_file = ;
% logfile = "D:\CodeManage\Robust-PSZ-algorithm\log\run_2025-10-27_23-41.log";
% 
% send_graphmail("zhouleicqupt2016@outlook.com", ...
%                "Run finished. Log attached.", ...
%                "The job ended. See attached log file.", ...
%                'Attachments', logfile);
