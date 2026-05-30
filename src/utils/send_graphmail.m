function send_graphmail(to_addr, subject_txt, body_txt, varargin)
%SEND_GRAPHMAIL  通过 Microsoft Graph 发送邮件 (基于 send_notification.py)
%
% 基本用法 (最常见):
%   send_graphmail("you@example.com", "Job Done", "Task finished");
%
% 进阶用法 (带抄送/附件/HTML):
%   send_graphmail("you@example.com", "Alert", "<b>Overvoltage!</b>", ...
%                  'CC', ["boss@example.com","log@example.com"], ...
%                  'Attachments', ["D:\logs\daq.txt", "D:\plots\scope.png"], ...
%                  'HTML', true);
%
% 参数:
%   to_addr        收件人 (字符串)
%   subject_txt    邮件主题 (字符串)
%   body_txt       邮件正文（默认按纯文本发；如果 'HTML',true 则按 HTML 发送）
%
% name-value 可选项:
%   'CC'           string 或 string数组 / cellstr
%   'Attachments'  string 或 string数组 / cellstr
%   'HTML'         logical true/false (默认 false)
%
% 依赖:
%   1. 你机器上的 Python 环境中已成功运行过 send_notification.py 并授权
%   2. send_notification.py 使用 consumers authority，且 saveToSentItems=False
%   3. msal_token_cache_outlook.json 在同目录可用(或将自动创建/刷新)
%
% 注意:
%   - 如果 token 缓存被你删掉，下一次第一次调用可能会在 Python 输出里提示 device code 登录
%   - 本函数本身不会弹交互框，但你需要在那时手动去浏览器完成一次授权

    % ----------- 固定环境配置: 你的 Python 解释器和脚本路径 -----------
    pythonExe     = 'C:\Users\ZhouLei\anaconda3\envs\python3.9\python.exe';
    python_script = 'D:\CodeManage\Robust-PSZ-algorithm\src\utils\send_notification.py';

    % ----------- 解析可选参数 -----------
    p = inputParser;
    addParameter(p, 'CC', [], @(x) isstring(x) || ischar(x) || iscellstr(x));
    addParameter(p, 'Attachments', [], @(x) isstring(x) || ischar(x) || iscellstr(x));
    addParameter(p, 'HTML', false, @(x) islogical(x) || isnumeric(x));
    parse(p, varargin{:});

    cc_list        = p.Results.CC;
    attach_list    = p.Results.Attachments;
    html_flag      = logical(p.Results.HTML);

    % ----------- 把 CC/Attachments 变成命令行友好的字符串 -----------
    % Python 这边支持:
    %   --cc "a@x.com,b@y.com"
    %   --attach "C:\log.txt;C:\img.png"
    %   --html  (仅作为flag)
    %
    % 我们只在有值的时候才加这些参数
    extra_args = "";

    % CC
    if ~isempty(cc_list)
        cc_list = string(cc_list); % 统一成 string array
        cc_joined = strjoin(cc_list, ",");  % 用逗号拼接
        extra_args = extra_args + sprintf(' --cc "%s"', cc_joined);
    end

    % Attachments
    if ~isempty(attach_list)
        attach_list = string(attach_list);
        % 我们用分号分隔，以配合 Python 端 split_list_arg 里的正则 ([;,])
        attach_joined = strjoin(attach_list, ";");
        extra_args = extra_args + sprintf(' --attach "%s"', attach_joined);
    end

    % HTML flag
    if html_flag
        extra_args = extra_args + " --html";
    end

    % ----------- 拼最终命令行 -----------
    % 注意引号，确保路径/主题/正文里如果有空格也能正确传进 Python
    % cmd = sprintf('"%s" "%s" "%s" "%s" "%s"%s', ...
    %     pythonExe, ...
    %     python_script, ...
    %     to_addr, ...
    %     subject_txt, ...
    %     body_txt, ...
    %     extra_args);

    cmd = sprintf('"%s" "%s" "%s" "%s" "%s" %s', ...
    pythonExe, ...
    python_script, ...
    escape_shell(to_addr), ...
    escape_shell(subject_txt), ...
    escape_shell(body_txt), ...
    extra_args);


    % ----------- 执行并收集输出 -----------
    [status, out] = system(cmd);

    % 打印 Python 那边输出，方便调试
    % fprintf('--- send_notification.py output ---\n%s\n-------------------------------\n', out);

    % 若 Python 返回非 0，就认为失败（可能是 401 / 403 / 网络问题等）
    if status ~= 0
        warning('send_graphmail:failed', ...
            'Email send failed or Python returned error code %d.\n%s', status, out);
    end
end

function out = escape_shell(str)
    % 转义双引号，确保字符串安全地传入 system(cmd)
    out = strrep(str, '"', '\"');
end
