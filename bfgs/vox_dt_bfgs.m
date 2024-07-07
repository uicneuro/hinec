function [D,nsteps] = vox_dt_bfgs(nim, x, y, z, opts)
  arguments
    nim
    x (1,1)
    y (1,1)
    z (1,1)
    opts.Plot string {mustBeMember(opts.Plot, ["on", "off"])} = "off"
    opts.Verbose string {mustBeMember(opts.Verbose, ["on", "off"])} = "off"
    opts.Steps (1,1) = 2000
    opts.FontSize (1,1) = 16
  end

  enable_plot    = strcmp(opts.Plot, "on");
  enable_verbose = strcmp(opts.Verbose, "on");

  % Known values
  S  = reshape(nim.img_bi(x,y,z,:), [nim.size_bi 1]);
  S0 = nim.img_b0(x,y,z);
  b  = nim.bval(nim.bval >= nim.thrsh_b0);
  g  = nim.bvec(nim.bval >= nim.thrsh_b0, :);

  % We skip voxels where S0 == 0.
  if S0 == 0
    D = eye(3);
    return;
  end

  % Normalize. Change unit of b to 10^3 s/mm^2
  b = b ./ 1000;

  % Calculate Apparent Diffusion Coefficient(ADC)
  if enable_plot
    adc = log(S0./S)./b;
    adc_vec = g.*adc;
  end

  % Compare results with LSF
  if enable_plot
    Y = log(S0./S)./b;
    gx = g(:,1);
    gy = g(:,2);
    gz = g(:,3);
    H = [gx.^2 gy.^2 gz.^2 2.*gx.*gy 2.*gy.*gz 2.*gz.*gx];
    D_lsf = H\Y;
    D_lsf = reshapeDM(D_lsf);
    f_lsf = objf(sqrtm(D_lsf),S,S0,b,g);

    figure("Name", "LSF", "Position", [50 50 500 500]);
    ax_lsf = axes;
    plotADC(ax_lsf, adc_vec);
    plotD(ax_lsf, D_lsf);
    title(ax_lsf, "LSF", "FontSize", opts.FontSize*1.25);
    subtitle(ax_lsf, "f = " + f_lsf, "FontSize", opts.FontSize);
    xlabel(ax_lsf, "x", "FontSize", opts.FontSize);
    ylabel(ax_lsf, "y", "FontSize", opts.FontSize);
    zlabel(ax_lsf, "z", "FontSize", opts.FontSize);
    xlim(ax_lsf, [-1 1]);
    ylim(ax_lsf, [-1 1]);
    zlim(ax_lsf, [-1 1]);
  end
  

  % Initial conditions
  D_init = eye(3);
  X = sqrtm(D_init);
  B = eye(6);

  % Stopping conditions
  consec_small_dX = 0;
  consec_small_gf = 0;

  % Plot
  if enable_plot
    figure("Name", "Diffusion Tensor (BFGS)", "Position", [100 50 500 500])
    ax_d = axes;
    plotADC(ax_d, adc_vec);
    plotD(ax_d, D_init);
    title(ax_d, "Diffusion Tensor (BFGS)", "FontSize", opts.FontSize*1.25);
    xlabel(ax_d, "x", "FontSize", opts.FontSize);
    ylabel(ax_d, "y", "FontSize", opts.FontSize);
    zlabel(ax_d, "z", "FontSize", opts.FontSize);
    xlim(ax_d, [-1 1]); ylim(ax_d, [-1 1]); zlim(ax_d, [-1 1]); 
    view(ax_d, 3);
    set(ax_d, "NextPlot", "replacechildren");

    figure("Name", "BFGS Results", "Position", [600 50 800 500]);
    ax_f  = subplot(1,2,1);
    title(ax_f, "Objective Function", "FontSize", opts.FontSize*1.25);
    xlabel(ax_f, "Steps", "FontSize", opts.FontSize);
    ylabel(ax_f, "f", "FontSize", opts.FontSize);
    set(ax_f, "NextPlot", "replacechildren");
    fs  = zeros(opts.Steps, 1);

    ax_gf = subplot(1,2,2);
    title(ax_gf, "Gradient Function", "FontSize", opts.FontSize*1.25);
    xlabel(ax_gf, "Steps", "FontSize", opts.FontSize);
    ylabel(ax_gf, "gf", "FontSize", opts.FontSize);
    set(ax_gf, "NextPlot", "replacechildren");
    gfs = zeros(opts.Steps, 1);
  end

  % BFGS
  warning("off", "MATLAB:nearlySingularMatrix");
  warning("off", "MATLAB:singular");
  warning("off", "MATLAB:singularMatrix");
  for step=1:opts.Steps
    % How should parameters change?
    gf = gradf(X, S, S0, b, g);
    sum_gf = sum(gf);
    dX = -B\gf;

    if any(isnan(dX))
      break;
    end

    % Update parameters
    ndX = vecnorm(dX);
    if ndX > 100
      dX = dX ./ (ndX .* 0.01);
    end
    nX = norm(reshapeDV(X));
    if nX > 100
      X = X ./ (nX .* 0.01);
    end
    X = X + reshapeDM(dX);

    % Stopping conditions
    if -1e-1 < ndX && ndX < 1e-1
      consec_small_dX = consec_small_dX + 1;
    else
      consec_small_dX = 0;
    end
    if -1e1 < sum_gf && sum_gf < 1e1
      consec_small_gf = consec_small_gf + 1;
    else
      consec_small_gf = 0;
    end
    if consec_small_dX > 5 && consec_small_gf > 5 
      break;
    end

    if enable_verbose
      disp("Step: " + step);
      disp("gf  : " + sum_gf);
      disp("ndX : " + ndX);
      disp("objf: " + objf(X, S, S0, b, g));
      disp("sm_dX: " + consec_small_dX + ", sm_gf: " + consec_small_gf);
    end

    if enable_plot
      fs(step)  = objf(X, S, S0, b, g);
      gfs(step) = sum_gf;
      start = max(1,step-100);
      plot(ax_f, linspace(start, step, step-start+1), fs(start:step), '-k', 'LineWidth', 1);
      plot(ax_gf, linspace(start, step, step-start+1), gfs(start:step), '-k', 'LineWidth', 1);
      plotADC(ax_d, adc_vec);
      plotD(ax_d, X*X);
    end

    % Approximate Hessian
    gfn = gradf(X, S, S0, b, g);
    Y = gfn - gf;
    B = B + (Y*Y')/(Y'*dX) - B*(dX*dX')*B'/(dX'*B*dX);
    B = 0.5*(B+B');

    if any(isnan(B))
      break;
    end
    nsteps = step;
  end  % BFGS step
  warning("on", "MATLAB:nearlySingularMatrix");
  warning("on", "MATLAB:singular");
  warning("on", "MATLAB:singularMatrix");

  % ========= Store final result ==========
  D = X*X;
  D = D .* 1000; % Return units of b to 1s/mm^2

  if enable_verbose
    disp("D: ");
    disp(D);
    [Q,l]=eig(D);
    disp("Q: ");
    disp(Q);
    disp("l: ");
    disp(maxk([l(1,1) l(2,2) l(3,3)],3));
  end

  if enable_plot
    fs_end = find(fs, 1, "last");
    gfs_end = find(gfs, 1, "last");
    plot(ax_f, fs(1:fs_end), '-k', 'LineWidth', 1);
    plot(ax_gf, gfs(1:gfs_end), '-k', 'LineWidth', 1);

    plotADC(ax_d, adc_vec);
    plotD(ax_d, D);
    subtitle(ax_d, "f = "+objf(X,S,S0,b,g), "FontSize", opts.FontSize);
  end

  D = reshapeDV(D); % Output as vector for "nim_dt.m"
end % function

% Visualize the apparent diffusion coefficient(ADC)
function plotADC(ax, adc_vec)
  ori = -0.5.*adc_vec; % origin
  dst =       adc_vec; % destination
  quiver3(ax, ori(:,1), ori(:,2), ori(:,3), dst(:,1), dst(:,2), dst(:,3),...
    "-k", "LineWidth", 1, "AutoScale", "off", "ShowArrowHead", "off");
end

% Visualize the diffusion tensor to an ellipsoid
function plotD(ax, D)
  % Eigendecompose
  [Q,l] = eig(D);
  l_max = max(l, [], "all");
  l_nrm = [ l(1,1)./l_max l(2,2)./l_max l(3,3)./l_max ];

  % Diffusion along x, y, z
  [X,Y,Z] = ellipsoid(0, 0, 0, l_nrm(1), l_nrm(2), l_nrm(3), 10);

  % Diffusion along xy, yz, zx
  sz = size(X);
  for a=1:sz(1)
    for b=1:sz(2)
      T = [X(a,b) Y(a,b) Z(a,b)]';
      T = Q * T;
      X(a,b)=T(1);
      Y(a,b)=T(2);
      Z(a,b)=T(3);
    end
  end

  old_nxpt = ax.NextPlot;
  ax.NextPlot = 'add';

  % Draw diffusion ellipsoid
  surf(ax,X,Y,Z,"FaceAlpha",0.3);
  view(ax,3);
  
  % Draw eigenvectors scaled by eigenvalues
  Q = Q .* vecnorm(Q,2) .* l_nrm;
  ori = zeros(3,1);
  quiver3(ax, ori, ori, ori, Q(1,:)', Q(2,:)', Q(3,:)',...
    '-r', 'LineWidth', 2, 'AutoScale', 'off');

  ax.NextPlot = old_nxpt;
end

function f = objf(X, S, S0, b, g)
  Y = X*g';
  YtY = Y(1,:).^2 + Y(2,:).^2 + Y(3,:).^2;
  f = (log(S./S0) + b.* YtY') .^ 2;
  f = sum(f); % sum over each k
end

function gf = gradf(X, S, S0, b, g)
  Y = X*g';
  YtY = Y(1,:).^2 + Y(2,:).^2 + Y(3,:).^2;
  
  gx = g(1); gy = g(2); gz = g(3);
  dY2dX = (2.* [ gx 0 0; 0 gy 0; 0 0 gz; gy gx 0; 0 gz gy; gz 0 gx; ] * Y)';
  gf = (2.*(log(S./S0) + b.*YtY').*dY2dX);

  % Sum over each k
  gf = sum(gf)';
end

% Reshape diffusion tensor to a vector
function d = reshapeDV(D)
  arguments
    D (3,3)
  end
  d = [ D(1,1) D(2,2) D(3,3) D(1,2) D(2,3) D(1,3) ];
end

% Reshape a vector to a diffusion tensor
function D = reshapeDM(d)
  arguments
    d (6,1)
  end
  D = [d(1) d(4) d(6);...
       d(4) d(2) d(5);...
       d(6) d(5) d(3);];
end
