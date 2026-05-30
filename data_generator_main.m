%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2025, Lei Zhou
% All rights reserved.
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MAIN SCRIPT - Automated RIR Database Generation for Server
% Sets up environment, generates all RIR datasets, with logging and email notifications.

clear; clc; close all;

%% 1. CONFIGURATION
% -------------------------------------------------------------------------
addpath(genpath('src/'));
rng(2025); % For global reproducibility

% --- Logging & Notification ---
log_dir = 'logs';
if ~exist(log_dir, 'dir')
    mkdir(log_dir);
end

% 生成运行时间戳，如 "20251027_1912"
run_timestamp = datetime("now","Format","yyyyMMdd_HHmm");

% 日志完整路径，比如 logs/20251027_1912_generation_run.log
log_filename = fullfile(log_dir, sprintf('%s_generation_run.log', char(run_timestamp)));

% 打开日志文件，追加模式
log_fid = fopen(log_filename, 'a');
if log_fid == -1
    error('main_1:LogOpenFailed', '无法打开日志文件: %s', log_filename);
end

% 确保无论函数怎么退出，都会关掉这个fid
cleanupObj = onCleanup(@() safe_close(log_fid));

% 服务器后台跑 (nohup) 建议设 false，就不会刷到stdout/nohup.out
% 本地调试可以设 true，看实时输出
ECHO_TO_CONSOLE = false;

recipient_email = 'zhouleicqupt2016@outlook.com'; % <-- Set your email
enable_notifications = true;

% --- Modes to Generate ---
modes_to_generate = {'position'}; % 可调节接口：'snr', 'temperature', 

%% 2. MAIN EXECUTION
% -------------------------------------------------------------------------
try
    log_message(log_fid, '--- RIR Generation Script Started ---', 'INFO', ECHO_TO_CONSOLE);

    total_timer = tic;

    % --- Loop through each generation mode ---
    for i = 1:length(modes_to_generate)
        current_mode = modes_to_generate{i};
        log_message(log_fid, ...
            sprintf('--- Starting Mode: %s ---', upper(current_mode)), ...
            'INFO', ECHO_TO_CONSOLE);

        try
            % 生成数据的函数在这里
            generate_rir_database(current_mode,log_fid, ECHO_TO_CONSOLE);

            log_message(log_fid, ...
                sprintf('SUCCESS: "%s" mode completed.', current_mode), ...
                'INFO', ECHO_TO_CONSOLE);

        catch ME_mode
            % 这个分支：单个mode失败，但不会直接杀整个脚本
            error_msg_short = sprintf('ERROR in "%s" mode: %s', ...
                                      current_mode, ME_mode.message);
            log_message(log_fid, error_msg_short, 'ERROR', ECHO_TO_CONSOLE);

            % 记录栈（basic）
            log_message(log_fid, getReport(ME_mode, 'basic'), ...
                        'DEBUG', ECHO_TO_CONSOLE);
        end
    end

    % --- Final Success Notification ---
    elapsed_time_total_hours = toc(total_timer) / 3600;
    final_message = sprintf('All generation tasks finished in %.2f hours.', ...
                            elapsed_time_total_hours);

    log_message(log_fid, final_message, 'INFO', ECHO_TO_CONSOLE);

    if enable_notifications
        % 附件里把完整日志带走
        send_graphmail(recipient_email, ...
            '[MATLAB Job Finished] Success', ...
            final_message, ...
            'Attachments', string(log_filename));
    end

    % 如果你在本地调试时希望一个总结行，下面这个printf是安全的；
    % 如果nohup后台，stdout会进nohup.out
    if ECHO_TO_CONSOLE
        fprintf('\n--- All Tasks Finished Successfully ---\n');
    end

catch ME_main
    % --- Global Error Handling & Notification ---
    final_error_msg = sprintf([ ...
        'CRITICAL FAILURE: The script has stopped.\n\n' ...
        'Error: %s'], ME_main.message);

    log_message(log_fid, final_error_msg, 'FATAL', ECHO_TO_CONSOLE);

    % 更详细的报错（包含堆栈），去掉超链接，方便纯文本环境/邮件
    log_message(log_fid, ...
        getReport(ME_main, 'extended', 'hyperlinks', 'off'), ...
        'DEBUG', ECHO_TO_CONSOLE);

    if enable_notifications
        send_graphmail(recipient_email, ...
            '[MATLAB Job FAILED]', ...
            final_error_msg, ...
            {log_filename});
    end

    if ECHO_TO_CONSOLE
        fprintf('\n--- Script Terminated Due to a Critical Error ---\n');
    end

    rethrow(ME_main);
end

%% 3. LOCAL FUNCTIONS
% -------------------------------------------------------------------------
function safe_close(fid)
    % 安全关闭文件句柄（给 onCleanup 用）
    if ~(isscalar(fid) && isnumeric(fid) && fid == -1)
        try
            fclose(fid);
        catch
            % 清理阶段静默处理
        end
    end
end


