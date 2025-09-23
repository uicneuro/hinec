load sample.mat
% nim_plot(nim)                                    % Plot all voxels
indx = 50:60;
indy = 50:60;
indz = 20:30;
nim_plot(nim, 'indx', indx, 'indy', indy, 'indz', indz, 'figindex', 100, 'show_figure', true)