%% MAIN SCRIPT - Cross-Validation Experiment Runner
% This script systematically evaluates the robustness of sound field control
% algorithms under various environmental perturbations. It loads a
% pre-generated RIR/ATF database, runs a full N x N cross-validation,
% and saves the resulting performance matrices for later analysis.
%
% Author: ZhouLei
% Date: 2025/10/24

clear; clc; close all;

%% ========================================================================
%  1. EXPERIMENT CONFIGURATION
%  ========================================================================
%  Users should modify this section to define the experiment to be run.

% --- CHOOSE THE EXPERIMENT TO RUN ---
% This name must match the subfolder in 'data/SimulateRIR/'.
% Options: 'snr', 'temperature', 'position'
experiment_mode = 'snr';

% --- DEFINE ALGORITHMS TO TEST ---
% These names must have a corresponding case in the 'design_filters.m' function.
algorithms_to_test = {'MyRobustPSZ', 'ACC', 'PM', 'WCRACC', 'POTDC_RACC'};

% --- DEFINE METRICS TO EVALUATE ---
% These names must match the fields returned by 'evaluate_performance.m'.
metrics_to_evaluate = {'AC', 'NSRE', 'AE','Planarity ';


%% ========================================================================
%  2. SETUP & DATA LOADING
%  ========================================================================
fprintf('--- Starting Cross-Validation for "%s" experiment ---\n', experiment_mode);

% --- Add source code folders to the MATLAB path ---
addpath(genpath('src/'));

% --- Load the common metadata file ('para.mat') ---
data_dir = fullfile('data', 'SimulateRIR', experiment_mode);
para_file = fullfile(data_dir, 'para.mat');
if ~exist(para_file, 'file')
    error('Metadata file "para.mat" not found in %s. Please generate the data first.', data_dir);
end
load(para_file); % Loads the 'para' struct
fprintf('Loaded metadata from: %s\n', para_file);

% --- Get Experiment-Specific Parameters ---
switch experiment_mode
    case 'snr', param_vector = para.snr_vector;
    case 'temperature', param_vector = para.temperature_vector_celsius;
    case 'position', param_vector = para.perturbation_radii_m;
    otherwise, error('Invalid experiment mode specified.');
end
num_levels = length(param_vector);
freq_params = para.freq_params;

% --- Pre-allocate a structure to store all results ---
results = struct();
for k = 1:length(algorithms_to_test)
    for m = 1:length(metrics_to_evaluate)
        % Initialize each metric matrix with NaNs for better debugging
        results.(algorithms_to_test{k}).(metrics_to_evaluate{m}) = nan(num_levels, num_levels);
    end
end
fprintf('Initialized %dx%d results matrices for %d algorithms.\n\n', num_levels, num_levels, length(algorithms_to_test));


%% ========================================================================
%  3. MAIN CROSS-VALIDATION LOOP
%  ========================================================================
fprintf('Starting main evaluation loop...\n');
tic;

% --- Outer loop: Iterate through each "Design Condition" ---
for i = 1:num_levels
    design_param = param_vector(i);
    fprintf('\n================== Design Level %d/%d (Param: %.2f) ==================\n', i, num_levels, design_param);
    
    % --- Load the RIR/ATF data for the current design condition ---
    design_filename = get_data_filename(data_dir, experiment_mode, design_param);
    design_data = load(design_filename);
    
    % --- Design filters ONCE for this design condition using all algorithms ---
    all_filters = struct();
    for k = 1:length(algorithms_to_test)
        algo_name = algorithms_to_test{k};
        fprintf('  Designing filters with "%s"...\n', algo_name);
        all_filters.(algo_name) = design_filters(algo_name, design_data, freq_params);
    end
    
    % --- Inner loop: Iterate through each "Operating Condition" for evaluation ---
    for j = 1:num_levels
        operating_param = param_vector(j);
        fprintf('    Evaluating on Operating Level %d/%d (Param: %.2f)...\n', j, num_levels, operating_param);
        
        % --- Load the RIR/ATF data for the current operating condition ---
        % Avoid reloading if the operating condition is the same as the design condition
        if i == j
            operating_data = design_data;
        else
            operating_filename = get_data_filename(data_dir, experiment_mode, operating_param);
            operating_data = load(operating_filename);
        end
        
        % --- Evaluate all pre-designed filters on the current operating data ---
        for k = 1:length(algorithms_to_test)
            algo_name = algorithms_to_test{k};
            
            % 1. Get the filters designed for condition 'i'
            filters_w = all_filters.(algo_name);
            
            % 2. Evaluate their performance in condition 'j'
            performance = evaluate_performance(filters_w, operating_data, freq_params);
            
            % 3. Store all resulting metrics into the results structure
            for m = 1:length(metrics_to_evaluate)
                metric_name = metrics_to_evaluate{m};
                results.(algo_name).(metric_name)(i, j) = performance.(metric_name);
            end
        end
    end
end

elapsed_time = toc;
fprintf('\nMain evaluation loop finished. Total time: %.2f minutes.\n', elapsed_time / 60);

%% ========================================================================
%  4. SAVE RESULTS
%  ========================================================================
% Create a dedicated directory for the analysis results of this experiment
results_dir = fullfile('results', [experiment_mode, '_Analysis']);
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

% Save all relevant variables into a single .mat file
results_filename = fullfile(results_dir, 'performance_matrices.mat');
save(results_filename, 'results', 'para', 'algorithms_to_test', 'metrics_to_evaluate', 'experiment_mode');

fprintf('\n--- Cross-Validation Complete ---\n');
fprintf('All performance matrices have been saved to:\n  %s\n', results_filename);