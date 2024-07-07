function bfgs_test(nim, x, y, z)
  format compact;
  % Known values
  S =reshape(nim.img_bi(x, y, z, :), [nim.size_bi 1]);
  S0=nim.img_b0(x, y, z);
  min(S0, [], "all")
  b =nim.bval(nim.bval >= nim.thrsh_b0);
  % Set units to 10^-3 mm^2/s (b = b/1000);
  g =nim.bvec(nim.bval >= nim.thrsh_b0, :);

  close all;

  fig_lsf = figure("Name", "LSF", "Position", [100 300 600 600]);
  ax = axes();
  bv=g.*log(S0./S);
  bv=bv./vecnorm(g,2,2);
  quiver3(ax, -1.*bv(:,1), -1.*bv(:,2), -1.*bv(:,3), 2.*bv(:,1), 2.*bv(:,2), 2.*bv(:,3), '-k', 'LineWidth', 1, 'AutoScale', 'off', "ShowArrowHead", "off");
  title("LSF");
  xlim(ax, [-1 1]);
  ylim(ax, [-1 1]);
  zlim(ax, [-1 1]);
  ax.XTick = [-1 0 1];
  ax.YTick = [-1 0 1];
  ax.ZTick = [-1 0 1];
  view(ax, 3);


  % ================ LSF ================ 
  % figure("Name", "LSF", "Position", [800 200 600 600]);
  gx = g(:,1);
  gy = g(:,2);
  gz = g(:,3);
  H = [gx.^2 gy.^2 gz.^2 2.*gx.*gy 2.*gx.*gz 2.*gy.*gz];
  % Y = log(S0./S)./b;
  Y = log(S0./S);
  D = H\Y;
  D = reshapeD(D)
  plotD(D, ax);

  exportgraphics(fig_lsf, "figures/DT_LSF.png");

  if exist("figures/DT_BFGS.gif", "file") ~= 0
    delete "figures/DT_BFGS.gif";
  end

  % ================ BFGS ================ 
  fig_bfgs = figure("Name", "BFGS", "Position", [800 300 600 600]);
  ax = axes();

  % Initial conditions
  initX = sqrtm(D);
  if any(~isreal(initX))
    initX = eye(3);
    disp("Initial D is Identity.");
  end
  % initX = eye(3) ./ 1000;
  % initX = sqrtm(initX)
  initX = sqrtm(eye(3));
  [~, initl] = eig(initX);
  X=reshapeX(initX);
  B=eye(6);
  initX
  pause;

  objfs = zeros(1000, 1);
  gradfs = zeros(1000, 1);

  for k=1:2000
    % Plot steps
    X_ = reshapeD(X);
    D_ = X_ * X_;
    cla(ax);
    quiver3(ax, -1.*bv(:,1), -1.*bv(:,2), -1.*bv(:,3), 2.*bv(:,1), 2.*bv(:,2), 2.*bv(:,3), '-k', 'LineWidth', 1, 'AutoScale', 'off', "ShowArrowHead", "off");
    xlim(ax, [-1 1]); ylim(ax, [-1 1]); zlim(ax, [-1 1]);
    ax.XTick = [-1 0 1]; ax.YTick = [-1 0 1]; ax.ZTick = [-1 0 1];
    view(ax, 3);
    title(ax, "BFGS");
    subtitle("Initial: Identity");

    % Plot diffusion tensor
    plotD(D_, ax);

    % exportgraphics(fig_bfgs, "figures/DT_BFGS.gif", "Append", true);

    % Check eigenvalues
    [~, l] = eig(D_);
    l_ = [l(1,1) l(2,2) l(3,3)];
    if ~any(l < 0)
      % Stop if PD
      % D_
      % l
      % k
      % break;
    end

    % Update params
    gradf_ = gradf(X, S, S0, b, g);
    gradfs(k) = sum(gradf_);
    sum_gradf_ = sum(gradf_)
    dX = -(B\gradf_);
    vndX = vecnorm(dX);
    if vndX > 1
      dX = (dX ./ vndX) .* 1;
    end
    X = X + dX;
    dX
    Y  = gradf(X, S, S0, b, g) - gradf_;
    % B  = B + (Y*Y')/(Y'*dX) - B*(dX*dX')*B/(dX'*B*dX);
    B  = B + (Y*Y')/(Y'*dX) - B*(dX*dX')*B'/(dX'*B*dX);
    B = 0.5*(B+B');
    

    % BFGS with trust regions


    % Check error (difference)
    objf_ = objf(X, S, S0, b, g);
    objfs(k) = objf_;
    disp("objf("+k+") = "+objf_);

    % Stopping condition
    if (any(isnan(B)))
      disp("Done. ("+k+" steps)");
      break
    end

    % if k > 10 && -10 < sum_gradf_ && sum_gradf_ < 10
    %   sum_gradf_
    %   disp("Done. ("+k+" steps)");
    %   break
    % end
    % pause(.05);
  end

  figure;
  subplot(1, 2, 1);
  plot(objfs(1:find(objfs, 1, "last")));
  title("Objective Function"); xlabel("Steps"); ylabel("f");

  subplot(1, 2, 2);
  plot(gradfs(1:find(gradfs, 1, "last")));
  title("Gradient"); xlabel("Steps"); ylabel('$|| \nabla f ||$', 'Interpreter', 'latex');
end

% Objective function (squared error of difference)
function e = objf(X, S, S0, b, g)
  Y=reshapeD(X)*g';
  YtY=Y(1,:).^2 + Y(2,:).^2 + Y(3,:).^2;
  % e=(log(S./S0) + b.*YtY') .^ 2;
  e=(log(S./S0) + YtY').^2;
  e=sum(e);

  % X0 = zeros(6, 1);
  % e = (X - X0)' * (X - X0);
end

% Gradient of objective function
function ge = gradf(X, S, S0, b, g)
  % S: kx1, S0: -, b: kx1, g: 3xk, X: 3x3
  Y=reshapeD(X)*g';
  YtY=Y(1,:).^2 + Y(2,:).^2 + Y(3,:).^2;
  gx=g(1); gy=g(2); gz=g(3);
  dY2dX=(2.*[gx 0 0; 0 gy 0; 0 0 gz; gy gx 0; gz 0 gx; 0 gz gy;]*Y)';
  % ge=sum(2.*(log(S./S0) + b.*YtY').*b.*dY2dX)';
  ge = (2.*(log(S./S0) + YtY').*dY2dX);
  ge=sum(ge)';

  % X0 = zeros(6, 1);
  % ge = 2 * (X - X0);
end

function D = reshapeD(X)
  arguments
    X (6,1)
  end
  D = [X(1) X(4) X(5); ...
       X(4) X(2) X(6); ...
       X(5) X(6) X(3);];
end

function X = reshapeX(D)
  arguments
    D (3,3)
  end
  X = [D(1,1) D(2,2) D(3,3) D(1,2) D(1,3) D(2,3)]';
end

function plotD(D, ax)
  arguments
    D (3,3)
    ax
  end

  % Draw the diffusion tensor in an ellipsoid
  [Q,l]=eig(D);
  maxl=max([l(1,1) l(2,2) l(3,3)]);
  nl=[l(1,1)/maxl l(2,2)/maxl l(3,3)/maxl];
  [X,Y,Z]=ellipsoid(0, 0, 0, nl(1), nl(2), nl(3), 10);
  sz = size(X);
  for x=1:sz(1)
    for y=1:sz(2)
      A=[X(x,y) Y(x,y) Z(x,y)]';
      A=Q*A;
      X(x,y)=A(1);
      Y(x,y)=A(2);
      Z(x,y)=A(3);
    end
  end
  hold on;
  surf(ax, X,Y,Z,"FaceAlpha", 0.3);
  view(ax, 3);
  grid on;
  xlabel("x"); ylabel("y"); zlabel("z");

  % Plot eigenvectors, scaled by eigenvalues
  nQ = Q .* vecnorm(Q,2);
  nQ = nQ .* nl;
  quiver3(ax, zeros(3,1), zeros(3,1), zeros(3,1), nQ(1,:)', nQ(2,:)', nQ(3,:)', '-r', 'LineWidth', 2, 'AutoScale', 'off');
  hold off;
end
