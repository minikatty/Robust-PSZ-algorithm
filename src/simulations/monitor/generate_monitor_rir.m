function generate_monitor_rir(varargin)
% GENERATE_MONITOR_RIR - Generates a dense grid of RIRs for sound field visualization.
% This is a self-contained function that intelligently manages its own parallel pool.
% It generates a large .mat file containing RIRs from all speakers to all 
% points on a 2D grid, suitable for creating contour maps.
%
% Optional Name-Value Pair Arguments:
%   'array_file'    - string: Path to the array geometry file.
%   'output_file'   - string: Path to save the output .mat file.
%   'log_fid'       - file ID: File identifier for logging. Defaults to console (1).
%   'echo_console'  - logical: Whether to echo log messages to the console.
%   'grid_spacing'  - scalar: Spacing for the monitor grid (m).
%   'temperature'   - scalar: Room temperature in Celsius.
%   'beta'          - scalar: Wall reflection coefficient.
%   'margin'        - scalar: Margin from the walls (m).
%   'z_height'      - scalar: Height of the monitor grid (m).

    %% 1. Default Parameters & Input Parsing
    p = inputParser;
    addParameter(p, 'array_file', 'data/arrayGeometry/array.mat', @ischar);
    addParameter(p, 'output_file', 'data/MonitorGrid/gridRIR_data.mat', @ischar);
    addParameter(p, 'log_fid', nan, @isnumeric);
    addParameter(p, 'echo_console', true, @islogical);

    addParameter(p, 'grid_spacing', 0.05, @isnumeric);%目前是0.02
    addParameter(p, 'temperature', 20, @isnumeric);
    addParameter(p, 'beta', 0.3, @isnumeric);
    addParameter(p, 'fs', 16000, @isnumeric);
    addParameter(p, 'z_height', 1.6, @isnumeric);
    parse(p, varargin{:});

    % --- Assign parsed results to local variables ---
    array_file = p.Results.array_file;
    output_file = p.Results.output_file;
    log_fid = p.Results.log_fid;
    echo_console = p.Results.echo_console;
    fs = p.Results.fs;
    
    grid_spacing = p.Results.grid_spacing;
    temperature = p.Results.temperature;
    beta = p.Results.beta;
    % margin = p.Results.margin;
    z_height = p.Results.z_height;

      %% 2. Log File Management
 
    if isnan(log_fid)
        % --- Standalone Mode: Create and manage our own log file ---
        log_dir = 'logs';
        if ~exist(log_dir, 'dir'), mkdir(log_dir); end
        
        run_timestamp = datetime("now","Format","yyyyMMdd_HHmm");
        log_filename = fullfile(log_dir, sprintf('%s_monitor_rir_run.log', char(run_timestamp)));
        
        log_fid = fopen(log_filename, 'a');
        if log_fid == -1
            error('StandaloneRun:LogOpenFailed', '无法打开日志文件: %s', log_filename);
        end
        
        % Create a cleanup object to auto-close the file handle on exit
        cleanupObj = onCleanup(@() fclose(log_fid));
        
        fprintf('Running in standalone mode. Logging to: %s\n', log_filename);
    end
    
    log_message(log_fid, '========== Starting Monitor Grid RIR Generation ==========');

    %% 2. Load Data and Define Parameters
    % --- Load array geometry ---
    try
        array_data_wrapper = load(array_file);
        array_geom = array_data_wrapper.array;
    catch ME
        error('Failed to load array file: %s. Error: %s', array_file, ME.message);
    end
    speakers = array_geom.s;
    room_size = array_geom.roomSize;
    log_message(log_fid, sprintf('Loaded array file: %s (%d speakers)', array_file, size(speakers, 1)), 'INFO', echo_console);

    % --- Signal Processing Parameters from centralized config ---
    truncated_time = 128; % unit: ms, include early reverberation
    rir_len = floor(truncated_time/1e3 * fs); % length of RIR
    
    % --- Physical Parameters ---
    c = temp2speed(temperature);
    log_message(log_fid, sprintf('Temperature: %.1f C, Sound Speed: %.1f m/s', temperature, c), 'INFO', echo_console);
    
    %% 3. Generate Grid Points
    x_vec = 0 : grid_spacing : room_size(1);
    y_vec = 0 : grid_spacing : room_size(2);
    [X, Y] = meshgrid(x_vec, y_vec);
    grid_points = [X(:), Y(:), repmat(z_height, numel(X), 1)];
    
    grid_info.x_vec = x_vec;
    grid_info.y_vec = y_vec;
    grid_info.nx = length(x_vec);
    grid_info.ny = length(y_vec);
    
    N_spk = size(speakers, 1);
    N_grid = size(grid_points, 1);
    log_message(log_fid, sprintf('Generated grid: %d x %d = %d points (Spacing: %.1f cm)', grid_info.nx, grid_info.ny, N_grid, grid_spacing*100), 'INFO', echo_console);

    %% 4. RIR Computation (Parallelized with Self-Managed Pool)
    
    total_rirs = N_spk * N_grid;
    log_message(log_fid, sprintf('Starting computation of %d RIRs...', total_rirs), 'INFO', echo_console);
    
     % --- 简单的并行池管理 ---
    log_message(log_fid, 'Starting parallel pool for this task...', 'INFO', echo_console);
    parpool('local',8); % 直接启动一个新的池

    RIR = zeros(N_spk, N_grid, rir_len);    
    tic;
       % --- Main Parallel Loop ---
    parfor i = 1:N_spk
        % Create a temporary slice for the current speaker to improve parfor efficiency
        temp_rir_slice = zeros(N_grid, rir_len);
        spk_pos_i = speakers(i,:); 
        
        for j = 1:N_grid
            temp_rir_slice(j, :) = rir_generator(c, fs, grid_points(j,:), spk_pos_i, room_size, beta, rir_len);
        end
        
        RIR(i, :, :) = temp_rir_slice;
        
        % Progress indicator (less frequent in parfor, prints to main console)
        if mod(i, 300) == 0 && echo_console
            fprintf('  Speaker %d/%d processed.\n', i, N_spk);
        end
    end
    elapsed = toc;
    log_message(log_fid, sprintf('Computation finished! Time elapsed: %.1f minutes', elapsed/60), 'INFO', echo_console);
        % --- 正常结束时，关闭池 ---
    log_message(log_fid, 'Shutting down parallel pool.', 'INFO', echo_console);
    delete(gcp('nocreate'));
    
    %% 5. Save Results
    log_message(log_fid, sprintf('Saving results to: %s', output_file), 'INFO', echo_console);
    
    % Convert to single precision before saving to save space
    RIR = single(RIR);
    
    % Ensure the output directory exists
    output_dir = fileparts(output_file);
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    
    save(output_file, 'RIR', 'grid_points', 'grid_info', '-v7.3');
    
    file_info = dir(output_file);
    log_message(log_fid, sprintf('File saved successfully (Size: %.1f MB)', file_info.bytes/1024^2), 'INFO', echo_console);
    log_message(log_fid, '==========================================================', 'INFO', echo_console);
end
