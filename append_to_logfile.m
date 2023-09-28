% Copyright © 2023 Martin Schonger
% This software is licensed under the GPLv3.


function append_to_logfile(logfile_fid, files_to_log_cellarray)

for i = 1:length(files_to_log_cellarray)
    content_to_append = fileread(files_to_log_cellarray{i});
    fprintf(logfile_fid, '\n\n\n[BEGIN] %s\n\n%s\n[END] %s', files_to_log_cellarray{i}, content_to_append, files_to_log_cellarray{i});
end

end