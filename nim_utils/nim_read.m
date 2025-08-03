function nim = nim_read(datapath, opts)
arguments
    % Path to NIfTI-1 file
    datapath string

    % --- Optional file paths ---
    % Override for b-values file path
    opts.BvalPath string = ""
    % Override for b-vectors file path
    opts.BvecPath string = ""
    % Override for mask file path
    opts.MaskPath string = ""

    % --- Toggles for reading files ---
    % Should read binary mask data
    opts.Mask {mustBeMember(opts.Mask, ["on", "off"])} = "on"
    % Should read b-values from ".bval" files
    opts.Bval {mustBeMember(opts.Bval, ["on", "off"])} = "on"
    % Should read b-vectors from ".bvec" files
    opts.Bvec {mustBeMember(opts.Bvec, ["on", "off"])} = "on"

    % --- DTI parameters ---
    % Specify b-value threshold (inclusive) for determining b0 images
    opts.B0Threshold = 5
end

% Ensure datapath is a character vector for robust file operations
datapath = char(datapath);

% Read NIfTI-1 header
nim.hdr = niftiinfo(datapath);
nim.xdim = nim.hdr.ImageSize(1);
nim.ydim = nim.hdr.ImageSize(2);
nim.zdim = nim.hdr.ImageSize(3);
nim.size3 = prod(nim.hdr.ImageSize(1:3));
nim.size_b0 = 1;
nim.size_bi = 1;

% Read image data
nim.img = niftiread(datapath);
disp("Loaded " + nim.hdr.Version + " file of dimensions: " + nim.xdim + "*" + nim.ydim + "*" + nim.zdim);

% Strip the .nii.gz extension to get the base file path for fallbacks
if endsWith(datapath, '.nii.gz')
    base_path = datapath(1:end-7);
else
    base_path = datapath;
end

% Read the b-matrix
if opts.Bval == "on"
    bval_filepath = "";
    if opts.BvalPath ~= ""
        bval_filepath = opts.BvalPath;
    else
        bval_filepath = base_path + ".bval";
    end
    nim.bval = transpose(str2num(string(fileread(bval_filepath))));
  
  % Number of b0 images
  nim.size_b0 = sum(nim.bval == 0);
  nim.size_bi = nim.hdr.ImageSize(4)-nim.size_b0;
  
  % Separate b0 and bi images
  nim.img_b0 = mean(nim.img(:, :, :, nim.bval < opts.B0Threshold), 4);
  nim.img_bi = double(nim.img(:, :, :, nim.bval >= opts.B0Threshold));
  nim.thrsh_b0 = opts.B0Threshold;
  
  disp("Found b-values");
else
  nim.bval = zeros(1, 0);
  nim.size_b0 = 0;
  nim.size_bi = 0;
  nim.img_b0 = uint8(1, 0);
  nim.img_bi = uint8(1, 0);
end

if opts.Bvec == "on"
    bvec_filepath = "";
    if opts.BvecPath ~= ""
        bvec_filepath = opts.BvecPath;
    else
        bvec_filepath = base_path + ".bvec";
    end
  % Each line contains x, y, z elements
  lines = splitlines(fileread(bvec_filepath));
  gx = str2num(string(lines(1))); %#ok<*ST2NM>
  gy = str2num(string(lines(2)));
  gz = str2num(string(lines(3)));
  nim.bvec = transpose([ gx; gy; gz; ]);
  
  disp("Found b-vectors");
else
  nim.bvec = zeros(1, 0);
end

% Read mask
if opts.Mask == "on"
    mask_filepath = "";
    if opts.MaskPath ~= ""
        mask_filepath = opts.MaskPath;
    else
        mask_filepath = base_path + "_M.nii.gz";
    end
    nim.mask = niftiread(mask_filepath);
    disp("Found brain mask");
else
  nim.mask = zeros(1, 0);
end

end
