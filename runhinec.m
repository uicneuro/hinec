main nifti_sample/sample sample.mat
load sample.mat

% interpolation order
p = 3;
nim_interp(nim);