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

  % Read bval file and detect format
  bval_content = fileread(bval_filepath);
  bval_lines = splitlines(bval_content);
  bval_lines = bval_lines(~cellfun('isempty', bval_lines)); % Remove empty lines

  if length(bval_lines) == 1
      % Standard format: Single line with space-separated values
      % Example: "0 1000 1000 1000 ..."
      nim.bval = transpose(str2num(string(bval_lines(1))));
      disp("Found b-values (standard format: single line)");
  elseif length(bval_lines) > 1
      % IronTract format: One value per line
      % Example: "0" on line 1, "1601.25" on line 2, etc.
      nim.bval = zeros(length(bval_lines), 1);
      for i = 1:length(bval_lines)
          nim.bval(i) = str2double(string(bval_lines(i)));
      end
      disp("Found b-values (IronTract format: one value per line)");
  else
      error('Unable to read .bval file: %s', bval_filepath);
  end

  % Number of b0 images
  nim.size_b0 = sum(nim.bval < opts.B0Threshold);
  nim.size_bi = nim.hdr.ImageSize(4) - nim.size_b0;

  % Separate b0 and bi images
  nim.img_b0 = mean(nim.img(:, :, :, nim.bval < opts.B0Threshold), 4);
  nim.img_bi = double(nim.img(:, :, :, nim.bval >= opts.B0Threshold));
  nim.thrsh_b0 = opts.B0Threshold;
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

  % Read bvec file and detect format
  lines = splitlines(fileread(bvec_filepath));
  lines = lines(~cellfun('isempty', lines)); % Remove empty lines

  % Detect format by checking first line
  first_line_vals = str2num(string(lines(1))); %#ok<*ST2NM>

  if length(first_line_vals) == 3 && length(lines) > 3
      % IronTract format: One 3D vector per line (x y z)
      % Example: "0.2 0 0" on each line
      bvec_all = zeros(length(lines), 3);
      for i = 1:length(lines)
          bvec_all(i, :) = str2num(string(lines(i)));
      end
      disp("Found b-vectors (IronTract format: one vector per line)");

  elseif length(first_line_vals) > 3 || (length(first_line_vals) >= 1 && length(lines) == 3)
      % Standard format: 3 lines (X, Y, Z), each with space-separated values
      % Example: Line 1 = "0 0.2 0.4 ...", Line 2 = "0 0 0.2 ...", Line 3 = "0 0 0 ..."
      gx = str2num(string(lines(1)));
      gy = str2num(string(lines(2)));
      gz = str2num(string(lines(3)));
      bvec_all = transpose([ gx; gy; gz ]);
      disp("Found b-vectors (standard format: 3 lines of components)");

  else
      error('Unable to determine .bvec file format. Expected either:\n  - 3 lines (X, Y, Z components)\n  - N lines with 3 values each (x y z per line)');
  end

  % IMPORTANT: Filter bvec to match img_bi (non-b0 volumes only)
  % This ensures bvec dimensions match bval filtering in downstream code
  if opts.Bval == "on"
      nim.bvec = bvec_all(nim.bval >= opts.B0Threshold, :);
  else
      nim.bvec = bvec_all;
  end
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
