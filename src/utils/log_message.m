function log_message(log_fid, message, level, echo_to_console)
% LOG_MESSAGE - 写一行带时间戳的日志到已打开的日志文件，
%               可选地回显到标准输出（命令行 / nohup.out）
%
% 用法:
%   log_message(log_fid, 'started');
%   log_message(log_fid, 'bad thing happened', 'ERROR');
%   log_message(log_fid, 'silent write', 'INFO', false);

    if nargin < 3 || isempty(level)
        level = 'INFO';
    end
    if nargin < 4 || isempty(echo_to_console)
        echo_to_console = true;
    end

    % 1. 时间戳（到秒）
    timestamp = datetime("now", "Format", "yyyy-MM-dd HH:mm:ss");

    % 2. 生成一整行日志
    %    e.g. [2025-10-27 19:40:12.387] [INFO] Program started
    log_entry = sprintf('[%s] [%s] %s\n', ...
                        char(timestamp), upper(level), char(message));

    % 3. 写入日志文件
    try
        % 检查 log_fid 是否有效（正整数句柄）
        if ~(isscalar(log_fid) && isnumeric(log_fid) && log_fid > 0)
            warning('log_message:InvalidFID', ...
                    'Failed to write to log file: invalid file handle.');
        else
            fprintf(log_fid, '%s', log_entry);

            % 旧版 MATLAB 没有 fflush(fid)，
            % 所以我们不强制flush，等 fclose 的时候自然落盘。
            % 如果你以后升级到带 fflush 的版本，可手动加回：
            %    if exist('fflush','builtin') || exist('fflush','file')
            %        fflush(log_fid);
            %    end
        end
    catch ME
        warning('log_message:WriteFailed', ...
                'Failed to write to log file: %s', ME.message);
    end
    
    % 4. 打到标准输出
    %    （服务器 nohup 模式你把 echo_to_console 设成 false，
    %     这样就不会把所有日志复制到 nohup.out 里）
    if echo_to_console
        fprintf('%s', log_entry);
    end
end
