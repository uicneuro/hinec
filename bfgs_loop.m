function bfgs_loop(nim, x, y, z, bfgs_max_step)
  arguments
    nim
    x (1,1)
    y (1,1)
    z (1,1)
    bfgs_max_step (1,1) = 10000
  end
  base_font_size = 16;
  % format compact;
  close all;
  disp("Unit of D: 10^-3 mm^2/s");

  % Known values
  S  = reshape(nim.img_bi(x,y,z,:), [nim.size_bi 1]);
  S0 = nim.img_b0(x,y,z);
  b  = nim.bval(nim.bval >= nim.thrsh_b0);
  g  = nim.bvec(nim.bval >= nim.thrsh_b0, :);
  k_size = nim.size_bi;


  % Change unit to 10^-3 mm^2/s
  b = b./1000;

  % Calculate ADC
  adc = log(S0./S)./b;
  adc_vec = g.*adc;

  % ========== LSF ==========
  gx = g(:,1); gy = g(:,2); gz = g(:,3); % normalized gradient vectors
  H = [gx.^2 gy.^2 gz.^2 2.*gx.*gy 2.*gy.*gz 2.*gz.*gx];
  Y = log(S0./S)./b;
  d = H\Y;
  D_lsf = reshape_dD(d);
  disp("[LSF] Approximated D");
  disp(D_lsf);
  X = sqrtm(D_lsf);
  disp("[LSF] Approximated X=sqrtm(D):");
  disp(X);
  
  err_lsf = 0;
  for k=1:k_size
    e = objf(X, S(k), S0, b(k), g(k,:)');
    err_lsf = err_lsf + e;
  end

  gerr_lsf = 0;
  for k=1:k_size
    gerr = sum(gradf(X, S(k), S0, b(k), g(k,:)'));
    gerr_lsf = gerr_lsf + gerr;
  end
  

  % ======== 1. LSF =========
  fig_lsf = figure("Name", "LSF", "Position", [100 200 600 600]);
  ax_lsf = axes;
  plotD(D_lsf, ax_lsf);
  hold on; plotADC(adc_vec, ax_lsf); hold off;
  xlabel("x", "FontSize", base_font_size*1.125);
  ylabel("y", "FontSize", base_font_size*1.125);
  zlabel("z", "FontSize", base_font_size*1.125);
  title("LSF", "FontSize", base_font_size*1.5);
  subtitle("Matlab \\ Operator. $\phi = " + err_lsf + "$", "Interpreter", "latex", "FontSize", base_font_size);
  pause;
  
  % ======== 2. BFGS ========
  
  % Initial condition
  X = sqrtm(D_lsf);
  if any(~isreal(X))
    X = sqrtm(eye(3));
  end
  disp("[BFGS] Initial X");
  disp(X);
  B = eye(6); % Approx. of Hessian

  % Set up plot
  fig_bfgs = figure("Name", "BFGS", "Position", [750 200 600 600]);
  ax_bfgs = axes;
  fig_objf = figure("Name", "Objective Function", "Position", [100 0 600 800]);
  ax_objf = subplot(2,1,1);
  title("Objective Function", "FontSize", base_font_size*1.5);
  xlabel("Step", "FontSize", base_font_size);
  ylabel("$\phi$", "Interpreter", "latex", "FontSize", base_font_size);
  ax_gradf = subplot(2,1,2);
  title("Gradient Function", "FontSize", base_font_size*1.5);
  xlabel("Step", "FontSize", base_font_size);
  ylabel("$\nabla \phi$", "Interpreter", "latex", "FontSize", base_font_size);


  % For each step
  errs = zeros(bfgs_max_step, 1);
  gerrs = zeros(bfgs_max_step, 1);

  for step=1:bfgs_max_step
    D = X * X;
    cla(ax_bfgs);
    plotADC(adc_vec, ax_bfgs);
    xlim(ax_bfgs, [-1 1]); ylim(ax_bfgs, [-1 1]); zlim(ax_bfgs, [-1 1]);
    view(ax_bfgs, 3);
    title(ax_bfgs, "BFGS", "FontSize", base_font_size*1.5);
    plotD(D, ax_bfgs);

    gerr = 0;
    gradf_ = zeros([k_size 6]);
    for k=1:k_size
      gerr_k = gradf(X, S(k), S0, b(k), g(k,:)');
      gerr = gerr + gerr_k;
      gradf_(k, :) = gerr_k;
    end
    gradf_ = sum(gradf_)';
    dX = -B\gradf_
    x = reshape_Dd(X) + dX;
    size(x)
    X = reshape_dD(x);

    err = 0;
    for k=1:k_size
      err_k = objf(X, S(k), S0, b(k), g(k,:)');
      err = err + err_k;
    end
    
    errs(step)  =  err;
    gerrs(step) = gerr;
    subtitle(ax_bfgs, "$\phi="+err+"$, $\nabla \phi="+gerr+"$", "Interpreter", "latex", "FontSize", base_font_size);

    plot(ax_objf, errs(1:step));
    plot(ax_gradf, errs(1:step));
    
    gerr1 = 0;
    gradf1_ = zeros([k_size 6]);
    for k=1:k_size
      gerr1_k = gradf(X, S(k), S0, b(k), g(k,:)');
      gerr1 = gerr1 + gerr1_k;
      gradf1_(k, :) = gerr1_k;
    end

    Y = gradf1_ - gradf_;
    B = B + (Y*Y')/(Y'*dX) - B*(dX*dX')*B/(dX'*B*dX);
  end
end

function err = objf(X, S, S0, b, g)
  arguments
    X
    S (1,1)
    S0 (1,1)
    b (1,1)
    g (3,1)
  end
  % d: 6x1, g: 3x1, everything else: scalar
  if any(size(X) ~= [3 3])
    X = reshape_dD(X);
  end
  Y = X*g;
  err = (log(S/S0) + b.*Y'*Y)^2;
end

function gerr = gradf(X, S, S0, b, g)
  arguments
    X
    S (1,1)
    S0 (1,1)
    b (1,1)
    g (3,1)
  end
  % d: 6x1, g: 3x1, everything else: scalar
  Y = X*g; % 3x3 * 3x1 = 3x1
  gx = g(1,1); gy = g(2,1); gz = g(3,1);
  dYdX = [ gx   0   0 ;...
            0  gy   0 ;...
            0   0  gz ;...
           gy  gx   0 ;...
            0  gz  gy ;...
            gz  0  gx ; ];

  dY2dX = dYdX * Y;
  gerr = 2 .* (log(S/S0) + b.*(Y'*Y)) .* b .* dY2dX;
  % gerr = sum(gerr_X);
end

function plotADC(adc_vec, ax)
  o = -1.*adc_vec; % origin
  d =  2.*adc_vec; % direction
  quiver3(ax, o(:,1),o(:,2),o(:,3),d(:,1),d(:,2),d(:,3),...
    "-k", "LineWidth", 1, "AutoScale", "off", "ShowArrowHead", "off");
end

function plotD(D, ax)
  arguments
    D (3,3)
    ax
  end

  [Q,l] = eig(D);
  l_max = max(l, [], "all");
  l_nrm = [ l(1,1)./l_max l(2,2)./l_max l(3,3)./l_max ];

  % Diffusion along x, y, z
  [X,Y,Z] = ellipsoid(0, 0, 0, l_nrm(1), l_nrm(2), l_nrm(3), 10);

  % Diffusion along xy, yz, zx
  sz = size(X);
  for x=1:sz(1)
    for y=1:sz(2)
      A = [X(x,y) Y(x,y) Z(x,y)]';
      A = Q * A;
      X(x,y)=A(1);
      Y(x,y)=A(2);
      Z(x,y)=A(3);
    end
  end

  % Draw diffusion ellipsoid
  hold on;
  surf(ax, X, Y, Z, "FaceAlpha", 0.3);
  view(3);
  grid on;

  % Draw eigenvectors, scaled by eigenvalues
  Q = Q .* vecnorm(Q,2) .* l_nrm;
  origin = zeros(3,1);
  quiver3(ax, origin, origin, origin, Q(1,:)', Q(2,:)', Q(3,:)', '-r', 'LineWidth', 2, 'AutoScale', 'off');
  hold off;
end

function D = reshape_dD(d)
  arguments
    d (6,1)
  end
  D = [ d(1) d(4) d(6);...
        d(4) d(2) d(5);...
        d(6) d(5) d(3);];
end

function d = reshape_Dd(D)
  arguments
    D (3,3)
  end
  d = [D(1,1) D(2,2) D(3,3) D(1,2) D(2,3) D(1,3)]';
end
