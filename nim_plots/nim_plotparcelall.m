function nim_plotparcelall(nim)
arguments
  nim
end

% Get the number of parcels
num_parcels = max(nim.parcellation_mask(:));

% Plot each parcel (starting from 1, skip background = 0)
figindex = 1;
for parcel_id = 1:num_parcels
  % Check if this parcel actually has voxels
  if sum(nim.parcellation_mask(:) == parcel_id) > 0
    nim_plotparcel(nim, parcel_id, figindex);
    figindex = figindex + 1;
  end
end
end
