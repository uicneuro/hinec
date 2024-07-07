function nim_new = nim_dt(nim, opts)
  arguments
    nim
    opts.Mask {mustBeMember(opts.Mask, ["on", "off"])} = "on"
  end
  enable_mask = opts.Mask == "on";
  if ~enable_mask
    disp("Brain mask is disabled. Analyzing whole brain.");
  end

  disp("Calculating diffusion tensors...");

  % Apparent diffusion coefficient(ADC)
  dim4 = [nim.hdr.ImageSize(1:3) nim.size_bi ];
  img_b_0 = mean(nim.img(:, :, :, nim.bval == 0), 4);
  img_b_i = double(nim.img(:, :, :, 2:end));
  Y = double(zeros(dim4));
  i = 1;
  for t=2:nim.hdr.ImageSize(4)
    if nim.bval(t) ~= 0
      Y(:, :, :, i) = log(img_b_0 ./ img_b_i(:, :, :, t-1)) ./ nim.bval(t);
      i = i + 1;
    end
  end

  % b-matrix (magnetic gradient)
  g_x = nim.bvec(:, 1);
  g_y = nim.bvec(:, 2);
  g_z = nim.bvec(:, 3);
  H = [ g_x .* g_x, g_y .* g_y, g_z .* g_z, ...
    g_x .* g_y .* 2, g_y .* g_z * 2, g_z .* g_x * 2 ];
  H = H(nim.bval > 0, :);

  % Fit diffusion tensor
  warning("off", "MATLAB:rankDeficientMatrix");
  nim.DT = double(zeros([nim.hdr.ImageSize(1:3) 6]));
  i = 0;
  cnt_neg = 0;
  for x=1:nim.xdim
    for y=1:nim.ydim
      for z=1:nim.zdim
        if enable_mask && nim.mask(x, y, z) == 0
          nim.DT(x, y, z, :) = zeros(1, 6);
        else
          Y_i = reshape(Y(x, y, z, :), [nim.size_bi 1]);

          % Solve for Y = Hd
          nim.DT(x, y, z, :) = H \ Y_i;
          i = i + 1;

        end
      end
    end
  end
  warning("on", "MATLAB:rankDeficientMatrix");
  disp("Fit diffusion tensors for " + i + " points.");
  disp(cnt_neg + "/" + i + " (" + (cnt_neg / i * 100) + "%) tensors were not PSD.");
  nim_new = nim;
end

