function nim_plotparcel(nim, parcel_id, figindex)
arguments
    nim
    parcel_id
    figindex = 1;
end

Nvox_x = nim.hdr.ImageSize(1);
Nvox_y = nim.hdr.ImageSize(2);
Nvox_z = nim.hdr.ImageSize(3);

% Extract the parcellation mask and the brain mask
mask = nim.parcellation_mask;
brain_mask = nim.mask > 0;

% Find the indices of the voxels belonging to the specified parcel within the brain mask
[indx, indy, indz] = ind2sub(size(mask), find(mask == parcel_id & brain_mask));

% Vertex points
[X,~] = zwuni(Nvox_x);     % [-1 1]
[Y,~] = zwuni(Nvox_y);
[Z,~] = zwuni(Nvox_z);

X = X .* floor(Nvox_x / 2);  % Scale. [-Nvox/2 Nvox/2]
Y = Y .* floor(Nvox_y / 2);
Z = Z .* floor(Nvox_z / 2);

% Center points
Xc = 0.5 .* (X(2:end) + X(1:end-1));
Yc = 0.5 .* (Y(2:end) + Y(1:end-1));
Zc = 0.5 .* (Z(2:end) + Z(1:end-1));

% Vectors for each voxel
Nvox = nim.hdr.ImageSize(1:3);
Vx = reshape(nim.evec(:, :, :, 1, 1), Nvox);
Vy = reshape(nim.evec(:, :, :, 1, 2), Nvox);
Vz = reshape(nim.evec(:, :, :, 1, 3), Nvox);

% Voxel Vertices
[XXc, YYc, ZZc] = meshgrid(Xc, Yc, Zc);

% Convert subscripts to linear indices for consistent indexing
linear_indices = sub2ind(size(mask), indx, indy, indz);

XXcind = XXc(linear_indices);
YYcind = YYc(linear_indices);
ZZcind = ZZc(linear_indices);
Vxind = Vx(linear_indices);
Vyind = Vy(linear_indices);
Vzind = Vz(linear_indices);

% Get parcel name if available
parcel_name = '';
if isfield(nim, 'atlas_labels') && isfield(nim.atlas_labels, 'map')
    try
        if nim.atlas_labels.map.isKey(parcel_id)
            parcel_name = nim.atlas_labels.map(parcel_id);
        end
    catch
        % If any error occurs, continue without the name
    end
end

% Create figure title
if ~isempty(parcel_name)
    title_str = ['Parcel ' num2str(parcel_id) ': ' parcel_name];
else
    title_str = ['Parcel ' num2str(parcel_id)];
end

figure(figindex);
hold on;
quiver3(XXcind, YYcind, ZZcind, Vxind, Vyind, Vzind, 'AutoScale', 'off');
% hold off;

axis([-Nvox_x / 2, Nvox_x / 2, -Nvox_y / 2, Nvox_y / 2, -Nvox_z / 2, Nvox_z / 2]);
xlabel('X');
ylabel('Y');
zlabel('Z');
title(title_str);
drawnow;
pause;
end
