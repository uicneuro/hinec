% main nifti_sample/parcellation_sample/sample sample_parcellated.mat
load sample_parcellated.mat
% main nifti_sample/original_sample/sample sample.mat
% load sample.mat

% plot all
% nim_plot(nim, 'mode', 'grid');             % Plot all regions in grid
nim_plot(nim, 'mode', 'parcels');            % Plot all parcels in separate figures

% interpolation order
% p = 3;
% nim_interp(nim,p);