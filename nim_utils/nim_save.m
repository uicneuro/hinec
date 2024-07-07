function nim_save(nim, nimpath)
    arguments
        % nim struct
        nim

        % Path to save processed .mat file. Must end in `.mat`
        nimpath string
    end
    disp("Saving '" + nimpath + "'...");
    save(nimpath, 'nim');
    disp("Done");
end
