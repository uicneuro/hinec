function nim_plotall_opt(nim, downsample_factor)
arguments
  % nim struct
  nim
  downsample_factor = 2; % Default downsample factor
end

Nvox_x = nim.hdr.ImageSize(1);
Nvox_y = nim.hdr.ImageSize(2);
Nvox_z = nim.hdr.ImageSize(3);

figindex = 1;

hx = 16;
hy = 16;
hz = 7;
for i = 1:Nvox_x/hx
  for j = 1:Nvox_y/hy
    for k = 1:Nvox_z/hz
      
      indx = ((i-1)*hx) + (1:hx);
      indy = ((j-1)*hy) + (1:hy);
      indz = ((k-1)*hz) + (1:hz);
      
      nim_plot_opt(nim, indx, indy, indz, figindex, downsample_factor);
      figindex = figindex + 1;
      fprintf('Plotting %d of %d\n', figindex, (Nvox_x/hx)*(Nvox_y/hy)*(Nvox_z/hz));
    end
  end
end
end
