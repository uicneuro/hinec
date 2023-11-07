function nim = nim_read(datapath, opts)
  arguments
    % Path to NIfTI-1 file
    datapath string

    % Should read binary mask data from "_M" files
    opts.Mask {mustBeMember(opts.Mask, ["on", "off"])} = "on"

    % Should read b-values from ".bval" files
    opts.Bval {mustBeMember(opts.Bval, ["on", "off"])} = "on"

    % Should read b-vectors from ".bvec" files
    opts.Bvec {mustBeMember(opts.Bvec, ["on", "off"])} = "on"

    % Specify b-value threshold (inclusive) for determining b0 images
    % Default is not zero, to account for b0 images taken with low, nonzero b values
    % (e.g. Phillips scanners)
    opts.B0Threshold = 5
  end

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

  % Read the b-matrix
  if opts.Bval == "on"
    nim.bval = transpose(str2num(string(fileread(datapath + ".bval"))));

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
    % Each line contains x, y, z elements
    lines = splitlines(fileread(datapath + ".bvec"));
    gx = str2num(string(lines(1)));
    gy = str2num(string(lines(2)));
    gz = str2num(string(lines(3)));
    nim.bvec = transpose([ gx; gy; gz; ]);

    disp("Found b-vectors");
  else
    nim.bvec = zeros(1, 0);
  end

  % Read mask
  if opts.Mask == "on"
    nim.mask = niftiread(datapath + "_M");
    disp("Found brain mask");
  else
    nim.mask = zeros(1, 0);
  end

end

