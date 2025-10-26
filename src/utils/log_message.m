function log_message(log_filepath, message, level)
% LOG_MESSAGE - Writes a formatted message with a timestamp to a log file.
%
% Usage:
%   log_message('my_log.txt', 'Starting process...');
%   log_message('my_log.txt', 'An error occurred!', 'ERROR');

    if nargin < 3
        level = 'INFO'; % Default log level
    end

    % 1. Create the timestamp string
    timestamp = datetime("now", "Format", "yyyy-MM-dd HH:mm:ss");
    
    % 2. Format the full log entry
    log_entry = sprintf('[%s] [%s] %s\n', timestamp, upper(level), message);
    
    % 3. Open the file in 'append' mode ('a') and write the entry
    try
        fid = fopen(log_filepath, 'a');
        if fid == -1
            error('Could not open log file: %s', log_filepath);
        end
        fprintf(fid, log_entry);
        fclose(fid);
    catch ME
        warning(ME.identifier, 'Failed to write to log file: %s', ME.message);
    end
    % 4. Also display the message in the command window for real-time feedback
    fprintf(log_entry);
end