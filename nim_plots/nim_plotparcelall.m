function nim_plotparcelall(nim)
arguments
  nim
end

% Get the number of parcels
num_parcels = max(nim.parcellation_mask(:));

% Plot each parcel
figindex = 1;
for parcel_id = 1:num_parcels
  if parcel_id > 0 % Skip non-brain areas marked as -1
    nim_plotparcel(nim, parcel_id, figindex);
    figindex = figindex + 1;
  end
end
end
