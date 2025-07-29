function nim = nim_parcellation(nim, parcellation_file)
% nim_parcellation: Simple wrapper for nim_parcellation_fixed
% This maintains compatibility with existing code that calls nim_parcellation
arguments
    nim
    parcellation_file {mustBeFile} % Path to the parcellation NIfTI file
end

% Delegate to the more robust fixed version
nim = nim_parcellation_fixed(nim, parcellation_file);

end