function nim_out = nim_dt_spd(nim, opts)
  arguments
    nim
    opts.Mask {mustBeMember(opts.Mask, ["on", "off"])} = "on"
    opts.Steps (1,1) = 2000
    opts.FontSize (1,1) = 16
  end
  enable_mask = opts.Mask == "on";
  if ~enable_mask
    disp("Brain mask is disabled. Analyzing whole brain.");
  end

  disp("Calculating diffusion tensors (PD, BFGS)...");
  dt_start = datetime('now');

  % Known values (all voxels)
  b  = nim.bval(nim.bval >= nim.thrsh_b0);
  g  = nim.bvec(nim.bval >= nim.thrsh_b0, :);
  
  size4 = [nim.hdr.ImageSize(1:3) nim.size_bi];
  Y = zeros(size4);
  % A sample is usable only where log(S0/S) is real and finite, i.e. S0 > 0 and
  % S > 0. Externally preprocessed data (e.g. denoised/PreQual output) can hold
  % non-positive intensities; they are left out of that voxel's fit instead of
  % producing a complex tensor. Nothing is substituted for them.
  valid4 = false(size4);
  for t=1:nim.size_bi
    valid4(:,:,:,t) = nim.img_b0 > 0 & nim.img_bi(:,:,:,t) > 0;
    yt = zeros(nim.hdr.ImageSize(1:3));
    bt = nim.img_bi(:,:,:,t);
    v = valid4(:,:,:,t);
    yt(v) = log(nim.img_b0(v) ./ bt(v)) ./ b(t);
    Y(:,:,:,t) = yt;
  end
  min_samples = 7;                 % 6 tensor unknowns + 1
  n_vox_dropped_samples = 0;       % voxels fitted with some samples left out
  n_vox_too_few = 0;               % voxels left at zero: too few usable samples

  gx = g(:,1);  gy = g(:,2);  gz = g(:,3);
  H = [ gx.^2 gy.^2 gz.^2 2.*gx.*gy 2.*gy.*gz 2.*gz.*gx ];

  % Store eigenvectors/eigenvalues
  nim.evec = zeros([nim.hdr.ImageSize(1:3) 3 3]);
  nim.eval = zeros([nim.hdr.ImageSize(1:3) 3]);

  % DEBUG: Plot BFGS performance metrics
  % figure("Name", "BFGS: Performance", "Position", [50 50 500 500]);
  % ax = axes;
  % axis padded;
  % title(ax, "Number of Steps", "FontSize", opts.FontSize * 1.25);
  % subtitle(ax, "Total voxel count: 0", "FontSize", opts.FontSize);
  % xlabel(ax, "Voxel index", "FontSize", opts.FontSize);
  % ylabel(ax, "Number of Steps", "FontSize", opts.FontSize);
  % set(ax, "NextPlot", "replacechildren");

  bfgs_nvox = 0;
  bfgs_nsteps = zeros(nim.size3, 1);

  % For logging progress
  vox_i = 1;          % Voxel index
  vox_n = nim.size3;  % Total voxel count

  % Approximate diffusion tensor
  warning("off", "MATLAB:rankDeficientMatrix");
  for x=1:nim.xdim
    for y=1:nim.ydim
      for z=1:nim.zdim

        % Regions that are not the brain
        if enable_mask && nim.mask(x,y,z) == 0
          nim.DT(x,y,z,:) = zeros(1,6);
          nim.evec(x,y,z,:,:) = zeros(3,3);
          nim.eval(x,y,z,:) = zeros(3,1);

        % b0 is zero or negative, or too few usable samples to fit a tensor
        elseif nim.img_b0(x,y,z) <= 0 || nnz(valid4(x,y,z,:)) < min_samples
          if nim.img_b0(x,y,z) > 0
            n_vox_too_few = n_vox_too_few + 1;
          end
          nim.DT(x,y,z,:) = zeros(1,6);
          nim.evec(x,y,z,:,:) = zeros(3,3);
          nim.eval(x,y,z,:) = zeros(3,1);

         % For each voxel inside the brain,
        else

          % Solve using LSF on the usable samples only
          Y_i = reshape(Y(x,y,z,:), [nim.size_bi 1]);
          use = reshape(valid4(x,y,z,:), [nim.size_bi 1]);
          all_usable = all(use);
          if ~all_usable
            n_vox_dropped_samples = n_vox_dropped_samples + 1;
          end
          D = H(use,:) \ Y_i(use);

          % Ensure positive-definiteness
          [Q,l] = eig(nim_reshape_d(D));
          l = [l(1,1) l(2,2) l(3,3)];
          if any(l<0)
            % DEBUG: Plot performance
            % dt_start = datetime('now', 'Format', 'hh:mm:ss');
            % fprintf("[%s]  BFGS for voxel: x=%d,y=%d,z=%d ",string(dt_start),x,y,z);

            % If not positive definite, solve using BFGS (on the usable samples only)
            if all_usable
              [D,nsteps] = vox_dt_bfgs(nim,x,y,z,"Plot","off","Verbose","off");
            else
              sub = struct('img_b0', nim.img_b0(x,y,z), ...
                           'img_bi', reshape(nim.img_bi(x,y,z,use), [1 1 1 nnz(use)]), ...
                           'bval', b(use), 'bvec', g(use,:), ...
                           'thrsh_b0', nim.thrsh_b0, 'size_bi', nnz(use));
              [D,nsteps] = vox_dt_bfgs(sub,1,1,1,"Plot","off","Verbose","off");
            end

            % DEBUG: Plot performance
            % dt_fin = datetime('now', 'Format', 'hh:mm:ss');
            % dt_dur = duration(dt_fin - dt_start, 'Format', 'hh:mm:ss');
            % fprintf("(Duration: %s) \r", string(dt_dur));

            % Record performance metrics
            bfgs_nvox = bfgs_nvox + 1;
            bfgs_nsteps(vox_i) = nsteps;

            % Re-compute eigenvalues
            [Q,l] = eig(nim_reshape_d(D));
            l = [l(1,1) l(2,2) l(3,3)];
          end

          % Store results
          nim.DT(x,y,z,:) = D;

          % Sort eigenvalues (descending)
          [lM, ilM] = maxk(l,3);
          QM = Q(:, ilM);

          nim.evec(x,y,z,:,:) = QM;
          nim.eval(x,y,z,:)   = lM;

         % Enhanced tensor validation (fix for NaN/Inf issues)
         if(~isempty(find(isnan(D)==1)) || ~isempty(find(isinf(D)==1)))
              nim.DT(x,y,z,:) = zeros(size(D));
              nim.evec(x,y,z,:,:) = zeros(size(QM));
              nim.eval(x,y,z,:)   = zeros(size(lM));
         end
         
         % Additional check for unrealistic eigenvalues
         if any(lM < 0) || any(lM > 0.01) % Typical diffusivity range
              nim.DT(x,y,z,:) = zeros(size(D));
              nim.evec(x,y,z,:,:) = zeros(size(QM));
              nim.eval(x,y,z,:)   = zeros(size(lM));
         end

        end

        % For progress logging
        vox_i = vox_i + 1;

      end %for z
    end % for y

    progress = floor((vox_i / vox_n) * 100);
    dt_progress = datetime('now', 'Format', 'hh:mm:ss');
    fprintf("[%s]  Voxel "+vox_i+"/"+vox_n+" ("+progress+" %%. BFGS: "+bfgs_nvox+") \r", string(dt_progress));
  end % for x
  warning("on", "MATLAB:rankDeficientMatrix");
  fprintf("\n");

  % ns_nz = bfgs_nsteps(bfgs_nsteps ~= 0);
  % ns_nz_sz = size(ns_nz, 1);
  % subtitle(ax, "Total voxel count: " + bfgs_nvox);
  % scatter(linspace(1,ns_nz_sz,ns_nz_sz),ns_nz,...
  %   20,"*",'black');

  dt_fin = datetime('now');
  elapsed = duration(dt_fin - dt_start, 'Format', 'hh:mm:ss');
  disp("Elapsed time: " + string(elapsed));
  if n_vox_dropped_samples > 0 || n_vox_too_few > 0
    fprintf("Non-positive intensities: %d voxels fitted without their unusable samples; %d left at zero (< %d usable samples).\n", ...
      n_vox_dropped_samples, n_vox_too_few, min_samples);
  end
  nim.dt_dropped_sample_voxels = n_vox_dropped_samples;
  nim.dt_too_few_sample_voxels = n_vox_too_few;
  nim_out = nim;
  disp("Done.");
end

