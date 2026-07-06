function nmi = compute_normalized_mutual_information(img1, img2)
% Compute normalized mutual information between two image vectors

% Number of bins for histogram
nbins = 64;

% Normalize intensities to [0, nbins-1]
img1_norm = round((img1 - min(img1)) / (max(img1) - min(img1)) * (nbins - 1)) + 1;
img2_norm = round((img2 - min(img2)) / (max(img2) - min(img2)) * (nbins - 1)) + 1;

% Ensure valid range
img1_norm = max(1, min(nbins, img1_norm));
img2_norm = max(1, min(nbins, img2_norm));

% Compute joint histogram
joint_hist = accumarray([img1_norm(:), img2_norm(:)], 1, [nbins, nbins]);
joint_hist = joint_hist / sum(joint_hist(:));

% Compute marginal histograms
hist1 = sum(joint_hist, 2);
hist2 = sum(joint_hist, 1);

% Compute entropies
H1 = -sum(hist1(hist1 > 0) .* log2(hist1(hist1 > 0)));
H2 = -sum(hist2(hist2 > 0) .* log2(hist2(hist2 > 0)));
H12 = -sum(joint_hist(joint_hist > 0) .* log2(joint_hist(joint_hist > 0)));

% Compute normalized mutual information
nmi = (H1 + H2) / H12;

end
