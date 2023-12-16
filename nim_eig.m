function nim_out = nim_eig(nim, opts)
  arguments
    nim
    opts.Mask {mustBeMember(opts.Mask, ["on", "off"])} = "on"
  end
  enable_mask = opts.Mask == "on";
  if ~enable_mask
    disp("Brain mask is disabled. Analyzing whole brain.");
  end

  disp("Eigendecomposing diffusion tensors...");
  nim.evec = zeros([nim.hdr.ImageSize(1:3) 3 3]);
  nim.eval = zeros([nim.hdr.ImageSize(1:3) 3]);
  i = 0;

  for x=1:nim.xdim
    for y=1:nim.ydim
      for z=1:nim.zdim
        % Apply the mask only when it is enabled.
        if ~enable_mask || nim.mask(x, y, z) == 1
          [Q, L] = eig(nim_reshape_d(nim.DT(x, y, z, :)));
          L = [ L(1, 1) L(2, 2) L(3, 3) ];
          
          % Sort eigenvalues (descending)
          [ML, IL] = maxk(L, 3);
          MQ = Q(:, IL);

          nim.evec(x, y, z, :, :) = MQ;
          nim.eval(x, y, z, :) = ML;
          i = i + 1;
        end
      end
    end
  end

  disp("Eigenvectors and eigenvalues calculated for " + i + " points.");
  nim_out = nim;

  disp("Done");
end
