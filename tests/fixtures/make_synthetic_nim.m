function nim = make_synthetic_nim(kind, varargin)
% make_synthetic_nim: build a nim struct DIRECTLY (no preprocessing) with an
% analytic principal-direction field, for VERIFICATION tests where the exact
% streamline / curvature is known. See MMF_Verification_Validation.md §2.1.
%
%   nim = make_synthetic_nim('straight')   % e1=[0 0 1], kappa=0 (exact = a line)
%   nim = make_synthetic_nim('arc')        % concentric circles, kappa=1/rho
%   nim = make_synthetic_nim('helix')      % helix: kappa,tau constant
%   nim = make_synthetic_nim('crossing')   % two planar-crossing bundles (x & y)
%
% Options (name/value): 'dims', 'center', 'Rin', 'Rout', 'lambda' ([l1 l2 l3]),
%   'helix_a', 'helix_b'.
%
% Returns nim with .evec [D1 D2 D3 3 3], .eval [D1 D2 D3 3], .FA, .mask (and, for
% 'arc'/'helix', the analytic parameters in nim.truth for the tests to compare to).

p = inputParser;
p.addParameter('dims', []);
p.addParameter('center', []);
p.addParameter('Rin', 10);
p.addParameter('Rout', 34);
p.addParameter('lambda', [1.7e-3 3e-4 3e-4]);
p.addParameter('helix_a', 15);
p.addParameter('helix_b', 6);
p.parse(varargin{:});
opt = p.Results;
lam = opt.lambda(:)';
FAs = fa_from_eig(lam);

switch lower(kind)
    case 'straight'
        dims = default(opt.dims, [24 24 60]);
        nim = blank_nim(dims);
        nim.evec(:,:,:,3,1) = 1;                 % e1 = [0 0 1]
        nim.evec(:,:,:,1,2) = 1;                 % e2 = [1 0 0]
        nim.evec(:,:,:,2,3) = 1;                 % e3 = [0 1 0]
        nim.eval = repmat(reshape(lam,1,1,1,3),[dims 1]);
        nim.FA(:) = FAs; nim.mask(:) = 1;
        nim.truth = struct('kind','straight','dir',[0 0 1]);

    case 'arc'
        dims = default(opt.dims, [80 80 8]);
        c = default(opt.center, [dims(1)/2, dims(2)/2]);
        nim = blank_nim(dims);
        [X,Y,~] = ndgrid(1:dims(1),1:dims(2),1:dims(3));
        rx = X - c(1); ry = Y - c(2); rho = sqrt(rx.^2+ry.^2); rg = rho; rg(rg<1e-6)=1;
        % e1 tangent to the circle = normalize([0 0 1] x r) = [-ry, rx, 0]/rho
        e1x = -ry./rg; e1y = rx./rg;
        er = rx./rg; ery = ry./rg;               % radial = e2
        nim.evec(:,:,:,1,1)=e1x; nim.evec(:,:,:,2,1)=e1y;         % e1
        nim.evec(:,:,:,1,2)=er;  nim.evec(:,:,:,2,2)=ery;        % e2 radial
        nim.evec(:,:,:,3,3)=1;                                    % e3 = z
        nim.eval = repmat(reshape(lam,1,1,1,3),[dims 1]);
        ann = (rho>=opt.Rin) & (rho<=opt.Rout);
        nim.mask = double(ann); nim.FA = FAs*double(ann);
        nim.truth = struct('kind','arc','center',c,'Rin',opt.Rin,'Rout',opt.Rout);

    case 'helix'
        dims = default(opt.dims, [60 60 60]);
        a = opt.helix_a; b = opt.helix_b;
        c = default(opt.center, [dims(1)/2, dims(2)/2]);
        nim = blank_nim(dims);
        [X,Y,Z] = ndgrid(1:dims(1),1:dims(2),1:dims(3));
        % nearest helix parameter theta for each voxel (helix axis = z, radius a)
        ang = atan2(Y-c(2), X-c(1));
        % tangent of helix (a cos t, a sin t, b t): T = [-a sin t, a cos t, b]/sqrt(a^2+b^2)
        s = sqrt(a^2+b^2);
        e1x = -a*sin(ang)/s; e1y = a*cos(ang)/s; e1z = b/s*ones(dims);
        nim.evec(:,:,:,1,1)=e1x; nim.evec(:,:,:,2,1)=e1y; nim.evec(:,:,:,3,1)=e1z;
        nim.eval = repmat(reshape(lam,1,1,1,3),[dims 1]);
        rho = sqrt((X-c(1)).^2+(Y-c(2)).^2);
        tube = abs(rho - a) <= 3;
        nim.mask = double(tube); nim.FA = FAs*double(tube);
        nim.truth = struct('kind','helix','a',a,'b',b, ...
            'kappa', a/(a^2+b^2), 'tau', b/(a^2+b^2));

    case 'crossing'
        % two planar-crossing bundles (A along x, B along y); crossing region is
        % planar with an adversarial principal eigenvector (undirected recovery).
        dims = default(opt.dims, [80 60 6]);
        bA = [27 33]; bB = [40 46];
        nim = blank_nim(dims);
        hi = lam(1); lo = lam(2); FAp = fa_from_eig([hi+lo hi+lo 2*lo]);
        xh=[1 0 0]; yh=[0 1 0]; zh=[0 0 1];
        for x=1:dims(1), for y=1:dims(2)
            inA = (y>=bA(1)&&y<=bA(2)); inB = (x>=bB(1)&&x<=bB(2));
            if inA&&inB, e1=yh;e2=xh;e3=zh; ev=[hi+lo hi+lo 2*lo]; fa=FAp;
            elseif inA, e1=xh;e2=yh;e3=zh; ev=[hi lo lo]; fa=FAs;
            elseif inB, e1=yh;e2=xh;e3=zh; ev=[hi lo lo]; fa=FAs;
            else, continue; end
            for z=1:dims(3)
                nim.FA(x,y,z)=fa; nim.mask(x,y,z)=1;
                nim.evec(x,y,z,:,1)=e1; nim.evec(x,y,z,:,2)=e2; nim.evec(x,y,z,:,3)=e3;
                nim.eval(x,y,z,:)=ev;
            end
        end, end
        nim.truth = struct('kind','crossing','bandA',bA,'bandB',bB);

    otherwise
        error('make_synthetic_nim: unknown kind "%s"', kind);
end
end

% ---- helpers ---------------------------------------------------------------
function nim = blank_nim(dims)
nim = struct();
nim.FA = zeros(dims); nim.mask = zeros(dims);
nim.evec = zeros([dims 3 3]); nim.eval = zeros([dims 3]);
nim.xdim=dims(1); nim.ydim=dims(2); nim.zdim=dims(3);
end
function v = default(v, d), if isempty(v), v = d; end, end
function fa = fa_from_eig(l), l=l(:); md=mean(l); fa=sqrt(1.5)*norm(l-md)/max(norm(l),1e-12); end
