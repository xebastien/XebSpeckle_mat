function progmeter(x, message)
% PROGMETER displays the progress of completion of a task in the MATLAB 
% command window. The meter is displayed as text that can be updated with
% subsequent calls. Its usage is similar to <a href="matlab: help
% waitbar">WAITBAR</a>.
% 
% PROGMETER(x, message) will create a new progress meter by displaying the string variable
% message and the progress represented as a percentage based on x. x can
% take values between 0 and 1 where 1 implies 100% completion. If a
% progress meter exists, this will update the existing progress meter with
% the new message and progress value calculated from x.
%
% PROGMETER(x) will update the progress meter with the new value in x without
% changing the message of an existing progress meter. If a new progress
% meter is created it will not display any message, just the progress.
%
% PROGMETER done will update the progress meter to "Done" and reset the internal state 
% PROGMETER clear will erase the meter and message and reset the internal state 
% PROGMETER reset will not update the meter but reset the internal state
%
% When the internal state is reset, subsequent calls to progmeter will
% display on a new line. Note that any other text displayed to the MATLAB
% command window between calls to progmeter can result in undesired
% behavior.
%
% EXAMPLES:
%
% SINGLE LOOP:
% clc
% progmeter(0, 'Running through pauses')
% for i = 1:20
%     pause(.2);
%     progmeter(i/20);
% end
% progmeter done
%
% MULTIPLE LOOPS:
% % This loop is 67% of the processing
% for i = 1:20
%     progmeter((i-1)/20*.67, 'Processing Task (loop 1)');
%     % Do some operation and pretend it takes .2 seconds
%     pause(.2);
% end
% % This loop is 33% of the processing
% for i = 1:10
%     progmeter((i-1)/10*.33 + .67, 'Processing Task (loop 2)');
%     % Do something else and pretend it takes .2 seconds
%     pause(.2);
% end
% progmeter clear

persistent clearlen messlen 
if isempty(messlen)
    messlen = 0;
    clearlen = 0;
end

% Remove previous meter text
if clearlen ~= 0 && ~(ischar(x) && strncmpi(x,'r',1)) % dont erase if a 'reset' is issued
    fprintf(repmat('\b', 1, clearlen));
end

if nargin > 1 % A new message has been input
    if messlen > 0 
        % erase previous message
        fprintf(repmat('\b', 1, messlen));
    end
    fprintf('%s: ',message);
    messlen = length(message) + 2;
end

if ischar(x) && ~isempty(x) % reset, clear or done
    clear clearlen
    clear messlen
    switch lower(x(1))
        case 'c'
            fprintf(repmat('\b', 1, messlen));
        case 'd'
            fprintf('Done\n');
        case 'r'
            fprintf('\n');
        otherwise
            error('The first input must be a numeric scalar between 0 and 1 or the strings ''done'', ''clear'' or ''reset''');
    end
elseif isnumeric(x) && ~isempty(x)
    progress = int2str(round(x*100));
    clearlen = length(progress) + 1;
    fprintf('%s%%', progress);
else
    error('The first input must be a numeric scalar between 0 and 1 or the strings ''done'', ''clear'' or ''reset''');
end