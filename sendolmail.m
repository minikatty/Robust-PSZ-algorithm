function sendolmail(to, subject, body, attachments)
    %Sends email using MS Outlook. Handles single or multiple recipients.
    
    % Create object and set parameters.
    h = actxserver('outlook.Application');
    mail = h.CreateItem('olMailItem');
    mail.Subject = subject;
    
    % --- 【最终修正】稳健地处理收件人 'to' ---
    if iscell(to)
        % If 'to' is a cell array, join its elements with a semicolon.
        % This correctly handles both {'a@a.com'} and {'a@a.com', 'b@b.com'}.
        to_string = strjoin(to, ';');
    elseif ischar(to) || isstring(to)
        % If 'to' is already a single string or char vector, use it directly.
        to_string = to;
    else
        % Handle incorrect input type.
        error('Recipient "to" must be a character vector, a string, or a cell array of character vectors.');
    end
    
    mail.To = to_string;
    
    mail.BodyFormat = 'olFormatHTML';
    mail.HTMLBody = body;
    
    % Add attachments, if specified.
    if nargin == 4 && ~isempty(attachments)
        for i = 1:length(attachments)
            mail.attachments.Add(attachments{i});
        end
    end
    
    % Send message and release object.
    mail.Send;
    h.release;
end