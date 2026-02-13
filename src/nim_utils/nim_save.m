function nim_save(nim, nimpath)
    arguments
        % nim struct
        nim

        % Path to save processed .mat file. Must end in `.mat`
        nimpath
    end

    % Convert to char array for save() function compatibility
    nimpath = char(nimpath);

    disp("Saving '" + nimpath + "'...");

    % Estimate size and use v7.3 format for large files (>2GB)
    % v7.3 format supports variables larger than 2GB and uses HDF5
    info = whos('nim');
    size_mb = info.bytes / 1024 / 1024;

    if size_mb > 1900  % Use 1900MB threshold to be safe
        fprintf('Large file detected (%.1f MB), using MAT-file v7.3 format...\n', size_mb);
        save(nimpath, 'nim', '-v7.3');
    else
        save(nimpath, 'nim');
    end

    disp("Done");
end
