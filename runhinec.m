main nifti_sample/parcellation_sample/sample sample_parcellated.mat
load sample_parcellated.mat
% main nifti_sample/original_sample/sample sample.mat
% load sample.mat

% plot all
% nim_plotall(nim);
% nim_plotall_opt(nim);
nim_plotparcelall(nim);

% interpolation order
% p = 3;
% nim_interp(nim,p);