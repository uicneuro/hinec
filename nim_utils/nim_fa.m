function nim_out = nim_fa(nim, opts)
  arguments
    nim
    opts.Mask {mustBeMember(opts.Mask, ["on", "off"])} = "on"
  end

  enable_mask = opts.Mask == "on";
  if ~enable_mask
    disp("Brain mask is disabled. Analyzing whole brain.");
  end

  disp("Calcuating fractional anisotropy (FA)...");
  nim.FA = zeros([nim.hdr.ImageSize(1:3) 1]);

  for x=1:nim.xdim
    for y=1:nim.ydim
      for z=1:nim.zdim
        if ~enable_mask || nim.mask(x, y, z) == 1
          l = nim.eval(x, y, z, :);
          denom = sqrt(l(1)^2 + l(2)^2 + l(3)^2);
          if denom > 0
            % Prevent divide by zero
            nim.FA(x, y, z) = sqrt(1/2) * sqrt( (l(1) - l(2))^2 + (l(2) - l(3))^2 + ...
              (l(3) - l(1))^2 ) / denom;
          end
        end
      end
    end
  end
  disp("Done.");
  nim_out = nim;
end
