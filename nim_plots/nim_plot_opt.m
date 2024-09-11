function nim_plot_opt(nim, indx, indy, indz, figindex, downsample_factor)
arguments
  % nim struct
  nim
  
  indx = 1:nim.hdr.ImageSize(1);
  indy = 1:nim.hdr.ImageSize(2);
  indz = 1:nim.hdr.ImageSize(3);
  figindex = 1;
  downsample_factor = 2; % Default downsample factor
end

Nvox_x = nim.hdr.ImageSize(1);
Nvox_y = nim.hdr.ImageSize(2);
Nvox_z = nim.hdr.ImageSize(3);

% Downsample indices
indx = indx(1:downsample_factor:end);
indy = indy(1:downsample_factor:end);
indz = indz(1:downsample_factor:end);

% Vertex points
[X,~] = zwuni(Nvox_x);     % [-1 1]
[Y,~] = zwuni(Nvox_y);
[Z,~] = zwuni(Nvox_z);

X = X .* floor(Nvox_x/2);  % Scale. [-Nvox/2 Nvox/2]
Y = Y .* floor(Nvox_y/2);
Z = Z .* floor(Nvox_z/2);

% Center points
Xc = 0.5.*(X(2:end) + X(1:end-1));
Yc = 0.5.*(Y(2:end) + Y(1:end-1));
Zc = 0.5.*(Z(2:end) + Z(1:end-1));

% Vectors for each voxel
Nvox = nim.hdr.ImageSize(1:3);
Vx = reshape(nim.evec(:, :, :, 1, 1), Nvox);
Vy = reshape(nim.evec(:, :, :, 1, 2), Nvox);
Vz = reshape(nim.evec(:, :, :, 1, 3), Nvox);

% Voxel Vertices
[XXc,YYc,ZZc] = meshgrid(Xc,Yc,Zc);

XXcind = XXc(indx,indy,indz);
YYcind = YYc(indx,indy,indz);
ZZcind = ZZc(indx,indy,indz);
Vxind = Vx(indx,indy,indz);
Vyind = Vy(indx,indy,indz);
Vzind = Vz(indx,indy,indz);

% Get colors based on direction
colors = vector_to_color(Vxind, Vyind, Vzind);

figure(figindex)
hold on
quiver3(XXcind, YYcind, ZZcind, Vxind, Vyind, Vzind, 'AutoScale', 'off');
for i = 1:numel(Vxind)
  quiver3(XXcind(i), YYcind(i), ZZcind(i), Vxind(i), Vyind(i), Vzind(i), 'Color', colors(i, :), 'LineWidth', 2);
end
hold off
axis([-Nvox_x/2, Nvox_x/2, -Nvox_y/2, Nvox_y/2, -Nvox_z/2, Nvox_z/2])
drawnow
pause
end
