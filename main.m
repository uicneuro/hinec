function main(imgpath, nimpath)
arguments
    % Path to NiFTI image file
    imgpath string
    
    % Path to save processed .mat file (must end in `.mat`)
    nimpath string
end

% include folders to path
addpath('nim_utils');
addpath('utils');
addpath('bfgs');


start_time = datetime('now', 'Format', 'yyyy-MM-dd hh:mm:ss');
fprintf("HINEC START: %s\n", string(start_time));

nim = nim_read(imgpath);
nim = nim_dt_spd(nim);
%nim = nim_eig(nim);
nim = nim_fa(nim);
nim_save(nim, nimpath);

end_time = datetime('now', 'Format', 'yyyy-MM-dd hh:mm:ss');
fprintf("HINEC END: %s\n", string(end_time));
end
