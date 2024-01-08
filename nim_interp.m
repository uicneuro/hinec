function nim_interp(nim,p)
  arguments
    % nim struct
    nim

    % interpolation order
    p = 3;
  end

  hdr

  Nvox_x = nim.hdr.ImageSize(1);
  Nvox_y = nim.hdr.ImageSize(2);
  Nvox_z = nim.hdr.ImageSize(3);

  [X,~] = zwuni(Nvox_x);     % [-1 1]
  [Y,~] = zwuni(Nvox_y);
  [Z,~] = zwuni(Nvox_z);

  X = X .* floor(Nvox_x/2);  % Scale. [-Nvox/2 Nvox/2]
  Y = Y .* floor(Nvox_y/2);
  Z = Z .* floor(Nvox_z/2);

  % N = nim.hdr.ImageSize(1:3) + 1;
  % Nx = N(1); Ny = N(2); Nz = N(3);

  % Vectors for each voxel
  Nvox = nim.hdr.ImageSize(1:3);
  Vx = reshape(nim.evec(:, :, :, 1, 1), Nvox);
  Vy = reshape(nim.evec(:, :, :, 1, 2), Nvox);
  Vz = reshape(nim.evec(:, :, :, 1, 3), Nvox);

  % Voxel Vertices
  [XX,YY,ZZ] = meshgrid(X,Y,Z);

  % Get dimensions
  NX = size(XX,1);
  NY = size(XX,2);
  NZ = size(XX,3);

  % Voxel center points
  xzcenter = 0.5.*(X(1:end-1) + X(2:end));
  yzcenter = 0.5.*(Y(1:end-1) + Y(2:end));
  zzcenter = 0.5.*(Z(1:end-1) + Z(2:end));
  [Xcenter,Ycenter,Zcenter] = meshgrid(xzcenter,yzcenter,zzcenter);

  % hold off;

  % Vectors for each grid point
  % Interpolate vectors
  Vxp = InterpCtoV3D(Vx);
  Vyp = InterpCtoV3D(Vy);
  Vzp = InterpCtoV3D(Vz);

  % Diffusion tensor at the center of each voxel
  % figure(1);
  % plotcubicgrid(NX,NY,NZ,XX,YY,ZZ);
  % quiver3(Xcenter,Ycenter,Zcenter, Vx, Vy, Vz, 'k-', 'LineWidth', 2);
  % title("DT vector at voxel center");
  % subtitle("Original data");
  % 
  % % Diffusion tensor at vertices
  % figure(2);
  % plotcubicgrid(NX,NY,NZ,XX,YY,ZZ);
  % quiver3(Xcenter,Ycenter,Zcenter, Vx, Vy, Vz, 'b-.', 'LineWidth', 2);
  % title("DT vectors at Voxel Vertices");
  % subtitle("Preprocessed");
  % hold off;
  % 
  % % Diffusion tensor interpolated to each collocation points
  % figure(3);
  % plotcubicgrid(NX,NY,NZ,XX,YY,ZZ); hold on;
  % title("DT vectors at collocation grid points (interpolated)");
  % subtitle("Interpolated (spline)");

  [zi,~] = zwuni(p);
  xi = 0.5*(1+zi); % on [0,1]
 
  % D = deriv_mat(xi);
  % Np = size(D(:,1),1);

  
  for ez=1:Nvox_z
    for ey=1:Nvox_y
      for ex=1:Nvox_x

        dx = X(ex+1)-X(ex);
        dy = Y(ey+1)-Y(ey);
        dz = Z(ez+1)-Z(ez);
        xf = X(ex) + xi*dx;  % 4x4x4 array of
        yf = Y(ey) + xi*dy;  % output points
        zf = Z(ez) + xi*dz;

        [Xf,Yf,Zf] = meshgrid(xf,yf,zf);  % Meshgrid needed for interp3() :(
        Vxf = interp3(X,Y,Z,Vxp,Xf,Yf,Zf,'spline');
        Vyf = interp3(X,Y,Z,Vyp,Xf,Yf,Zf,'spline');
        Vzf = interp3(X,Y,Z,Vzp,Xf,Yf,Zf,'spline');

        quiver3(Xf,Yf,Zf,Vxf,Vyf,Vzf, '-k', 'LineWidth', 2);
        axis equal

      end  % for ez
    end  % for ey
  end  % for ex

  hold off;

  savefile = 'nim_interp';
  save(savefile,'Xf','Yf','Zf','Vxf','Vyf','Vzf')

end


% 1D expansion
function v = InterpCtoV1D(vc)
  N0 = size(vc(:),1); N = N0+1;
  v = zeros(N,1);
  
  v(1)=vc(1);
  v(N)=vc(N-1);

  for i=2:N-1
    v(i) = 0.5*(vc(i-1)+vc(i));
  end
end


% 2D Expansion
function V = InterpCtoV2D(Vc)
  Nx0 = size(Vc,1);  Nx = Nx0+1;
  Ny0 = size(Vc,2);  Ny = Ny0+1;
  V = zeros(Nx,Ny);

  V(1,:)  = InterpCtoV1D(Vc(1,:));     % Y = -1.0
  V(Nx,:) = InterpCtoV1D(Vc(Nx-1,:));  % Y = 1.0
  V(:,1)  = InterpCtoV1D(Vc(:,1));     % X = -1.0
  V(:,Ny) = InterpCtoV1D(Vc(:,Ny-1));  % X = 1.0
  
  for i=2:Nx-1
    for j=2:Ny-1
      V(i,j) = 0.25 *( Vc(i-1,j-1) + Vc(i,j-1) + Vc(i-1,j) + Vc(i,j) );
    end
  end
end


% 3D Expansion
function V = InterpCtoV3D(Vc)
  Nx0 = size(Vc,1);  Nx = Nx0 + 1;
  Ny0 = size(Vc,2);  Ny = Ny0 + 1;
  Nz0 = size(Vc,3);  Nz = Nz0 + 1;
  
  % XY - plane
  V(:,:,1) = InterpCtoV2D(Vc(:,:,1));
  V(:,:,Nz) = InterpCtoV2D(Vc(:,:,Nz-1));

  % YZ - plane
  tmp1 = zeros(Ny0,Nz0);
  tmpN = zeros(Ny0,Nz0);
  for i=1:Ny0
    tmp1(i,:) = Vc(1,i,:);
    tmpN(i,:) = Vc(Nx0,i,:);
  end
  R1 = InterpCtoV2D(tmp1); 
  RN = InterpCtoV2D(tmpN); 
  for i=1:Ny
    V(1,i,:) = R1(i,:);
    V(Nx,i,:) = RN(i,:);
  end

  % XZ - plane
  tmp1 = zeros(Nx0,Nz0);
  tmpN = zeros(Nx0,Nz0);
  for i=1:Nx0
    tmp1(i,:) = Vc(i,1,:);
    tmpN(i,:) = Vc(i,Ny0,:);
  end
  R1 = InterpCtoV2D(tmp1); 
  RN = InterpCtoV2D(tmpN); 
  for i=1:Nx
    V(i,1,:) = R1(i,:);
    V(i,Nx,:) = RN(i,:);
  end

  % Linear interpolation
  for i=2:Nx-1
    for j=2:Ny-1
      for k=2:Nz-1
        V(i,j,k) = 0.125 *( Vc(i-1,j-1,k) + Vc(i,j-1,k)...
                           + Vc(i-1,j,k) + Vc(i,j,k) ...
                           + Vc(i-1,j-1,k-1) + Vc(i,j-1,k-1)...
                           + Vc(i-1,j,k-1) + Vc(i,j,k-1) );
      end  %k
    end  %j
  end  %i

end % InterpCtoV3D


function plotcubicgrid(NX,NY,NZ,XX,YY,ZZ)
  % XY - plane
  plot3(XX(:,:,1),YY(:,:,1),ZZ(:,:,1),'r-',XX(:,:,1)',YY(:,:,1)',ZZ(:,:,1)','r-');
  hold on;

  % YZ - plane
  for j=1:NY
    for k=1:NZ
        X1(:,k) = XX(:,j,k);
        Y1(:,k) = YY(:,j,k);
        Z1(:,k) = ZZ(:,j,k);
    end
    plot3(X1,Y1,Z1,'r-', X1',Y1',Z1','r-' );
    hold on;
  end
  axis equal;
  drawnow;

  % ZX - plane
  for j=1:NX
    for k=1:NZ
        X1(:,k) = XX(j,:,k);
        Y1(:,k) = YY(j,:,k);
        Z1(:,k) = ZZ(j,:,k);
    end
    plot3(X1,Y1,Z1,'r-', X1',Y1',Z1','r-' );
    hold on;
  end  % for j
  axis equal;

end  % function plotcubicgrid
