function nim = nim_mmf_from_dwi(nim, options)
%NIM_MMF_FROM_DWI  Moving frame + connection curvature estimated from the RAW DWI.
%
%   WHY THIS EXISTS. The normal route is DWI -> tensor -> principal eigenvector
%   -> finite-difference that field -> curvature. Measured against the ISMRM-2015
%   ground-truth geometry, that route obeys
%
%       kappa_field = a * kappa_true + c,     a ~ 0.4-0.5,  c ~ 0.015 /mm
%
%   i.e. it returns about HALF the true bending and adds a floor of spurious
%   bending where the truth is straight. Raising the derivative order does not
%   help: 2nd, 4th and 6th order stencils, and trilinear vs cubic interpolation,
%   all land within +-4% of the same gain. The information is not lost by the
%   derivative -- it is lost when 33 signal volumes are collapsed into one rank-1
%   direction per voxel, independently of the neighbours.
%
%   THE MODEL. For one fibre population the DW signal along gradient g is
%       log S(g) = C - A (g.u)^2 ,      A = b (lam_par - lam_perp) >= 0
%   and across a neighbourhood the direction is TRANSPORTED by the connection
%   rather than being a free parameter per voxel:
%       u(d) = normalise( e1 + (d.e1) * kappa ),      kappa _|_ e1
%   Centering the log-signal per voxel over gradients removes the unknown C(d)
%   (which absorbs S0 and the isotropic part), leaving five parameters
%       e1 (2 angles), kappa (2 components in the e2/e3 plane), A (1)
%   fitted to (#gradients x #neighbours) measurements -- 32 x 27 = 864 here.
%   Curvature is therefore ESTIMATED FROM THE SIGNAL, not differenced out of a
%   field that already threw it away.
%
%   RELATION TO ASYMMETRIC FODs. The DW signal is antipodally symmetric,
%   S(g) = S(-g), so any VOXEL-WISE fibre orientation distribution estimated from
%   it is a symmetric function on the sphere -- and a symmetric spherical function
%   cannot separate bending from fanning from crossing (fanning in and fanning out
%   give the identical FOD). Reisert et al. showed the ASYMMETRIC part of the FOD
%   is precisely the fibre curvature, and that it is recovered by imposing FIBRE
%   CONTINUITY across neighbours: a fibre leaving x along u must enter x+du along
%   u. This function is that construction written in connection coordinates --
%   (e1, kappa) is a 5-number parameterisation of the asymmetric part, against the
%   (lmax+1)^2 spherical-harmonic coefficients a general asymmetric FOD needs.
%
%   Continuity is imposed by relaxation (CONTINUITY_SWEEPS below): after the
%   independent per-voxel fits, each voxel is re-solved with a penalty tying the
%   direction its own connection PREDICTS at each neighbour to that neighbour's
%   estimated direction. Without it the per-voxel fits are mutually inconsistent:
%   they recover more curvature (0.89 of the true fibre curvature vs 0.65 for the
%   tensor route) but the excess is uncorrelated between voxels, and a tracker
%   integrating that field follows the noise -- measured as F1 0.395 -> 0.366 on
%   the ISMRM-2015 whole brain.
%
%   Stored as MMF the result is neither a tensor nor an FOD: it is the frame
%   field {e1,e2,e3} plus the connection (curvature vector, torsion).

if nargin < 2, options = struct(); end
fa_floor = getfielddef(options,'termination_fa',0.08);

dims = [nim.xdim nim.ydim nim.zdim];
bval = nim.bval(:); bvec = nim.bvec;
dwi  = double(nim.img);

isb0 = bval < 50;
g    = bvec(~isb0,:);                       % [G x 3] diffusion directions
g    = g ./ max(sqrt(sum(g.^2,2)),1e-12);
G    = size(g,1);
sig  = dwi(:,:,:,~isb0);
b0   = double(nim.img_b0);

% log-signal, centred per voxel over gradients (removes C(d))
L    = log(max(sig, 1e-6));
Lc   = L - mean(L,4);

% which voxels to fit
fit_mask = (nim.FA >= fa_floor) & (nim.mask > 0) & (b0 > 0);
idx = find(fit_mask);
fprintf('nim_mmf_from_dwi: fitting %d voxels (%d gradients)\n', numel(idx), G);

% neighbourhood offsets: a 3x3x3 block, weighted toward the fibre axis at fit time
[ox,oy,oz] = ndgrid(-1:1,-1:1,-1:1);
off = [ox(:) oy(:) oz(:)];
N   = size(off,1);

% initialisation from the existing tensor solution
E1i = nim.evec(:,:,:,:,1);
lam = nim.eval;
Ai  = max(mean(bval(~isb0)) * (lam(:,:,:,1)-lam(:,:,:,3)), 1e-3);

kap_out = zeros([dims 3]);
e1_out  = zeros([dims 3]);
resid   = nan(dims);

[sx,sy,sz] = size(nim.FA(:,:,:,1));
Lc_r = reshape(Lc, [], G);
lin  = @(i,j,k) i + (j-1)*sx + (k-1)*sx*sy;

nidx = numel(idx);
kap_l = zeros(nidx,3); e1_l = zeros(nidx,3); res_l = nan(nidx,1);
t0 = tic;
parfor t = 1:nidx
    [i,j,k] = ind2sub(dims, idx(t));
    if i<2||j<2||k<2||i>dims(1)-1||j>dims(2)-1||k>dims(3)-1
        e1_l(t,:) = squeeze(E1i(i,j,k,:))'; continue;
    end
    % gather the neighbourhood's centred log-signal
    Y = zeros(N,G); ok = false(N,1);
    for n = 1:N
        ii=i+off(n,1); jj=j+off(n,2); kk=k+off(n,3);
        if fit_mask(ii,jj,kk)
            Y(n,:) = Lc_r(lin(ii,jj,kk),:); ok(n)=true;
        end
    end
    if nnz(ok) < 8, e1_l(t,:) = squeeze(E1i(i,j,k,:))'; continue; end
    Yv = Y(ok,:); dv = off(ok,:);

    e10 = squeeze(E1i(i,j,k,:))'; e10 = e10/max(norm(e10),1e-12);
    p0  = [dir2ang(e10), 0, 0, Ai(i,j,k)];
    [p, r] = fit_frame(p0, dv, g, Yv);
    [e1f, kf] = unpack(p);
    e1_l(t,:) = e1f; kap_l(t,:) = kf; res_l(t) = r;
end
fprintf('nim_mmf_from_dwi: fit loop %.1f s\n', toc(t0));

% ---- FIBRE CONTINUITY RELAXATION (the A-FOD constraint, in connection form) ---
CONTINUITY_SWEEPS = getfielddef(options,'continuity_sweeps',3);
CONTINUITY_WEIGHT = getfielddef(options,'continuity_weight',1.0);
if CONTINUITY_SWEEPS > 0
    tc = tic;
    E1f = zeros(dims(1)*dims(2)*dims(3),3); E1f(idx,:) = e1_l;
    Kf  = zeros(dims(1)*dims(2)*dims(3),3); Kf(idx,:)  = kap_l;
    for sweep = 1:CONTINUITY_SWEEPS
        e1_prev = E1f; k_prev = Kf;
        e1_new = e1_l; k_new = kap_l;
        parfor t = 1:nidx
            [i,j,k] = ind2sub(dims, idx(t));
            if i<2||j<2||k<2||i>dims(1)-1||j>dims(2)-1||k>dims(3)-1, continue; end
            Y = zeros(N,G); okn = false(N,1); En = zeros(N,3);
            for n = 1:N
                ii=i+off(n,1); jj=j+off(n,2); kk=k+off(n,3);
                if fit_mask(ii,jj,kk)
                    li = lin(ii,jj,kk);
                    Y(n,:) = Lc_r(li,:); En(n,:) = e1_prev(li,:); okn(n)=true;
                end
            end
            if nnz(okn) < 8, continue; end
            e1c = e1_prev(lin(i,j,k),:); kc = k_prev(lin(i,j,k),:);
            if norm(e1c) < 1e-9, continue; end
            e1c = e1c/norm(e1c);
            [e2c,e3c] = perp_basis(e1c);
            p0 = [dir2ang(e1c), dot(kc,e2c), dot(kc,e3c), Ai(i,j,k)];
            p  = fit_frame_cont(p0, off(okn,:), g, Y(okn,:), En(okn,:), CONTINUITY_WEIGHT);
            [ef,kf2] = unpack(p);
            e1_new(t,:) = ef; k_new(t,:) = kf2;
        end
        e1_l = e1_new; kap_l = k_new;
        E1f(idx,:) = e1_l; Kf(idx,:) = kap_l;
    end
    fprintf('nim_mmf_from_dwi: %d continuity sweeps %.1f s\n', CONTINUITY_SWEEPS, toc(tc));
end

for c=1:3
    tmp = zeros(dims); tmp(idx) = kap_l(:,c); kap_out(:,:,:,c) = tmp;
    tmp = zeros(dims); tmp(idx) = e1_l(:,c);  e1_out(:,:,:,c)  = tmp;
end
resid(idx) = res_l;

% outside the fitted set fall back to the tensor eigenvector, zero curvature
fb = ~fit_mask;
for c=1:3
    t1 = e1_out(:,:,:,c); t2 = E1i(:,:,:,c); t1(fb) = t2(fb); e1_out(:,:,:,c) = t1;
end

nim.mmf_e1_dwi    = e1_out;
nim.mmf_kappa_dwi = kap_out;
nim.mmf_dwi_resid = resid;
km = sqrt(sum(kap_out.^2,4));
fprintf('nim_mmf_from_dwi: |kappa| median %.4f /voxel (%.4f /mm at %.1f mm)\n', ...
        median(km(fit_mask)), median(km(fit_mask))/2, 2.0);
end

% ---------------------------------------------------------------------------
function [p, rn] = fit_frame(p0, dv, g, Y)
% Levenberg-Marquardt on 5 parameters: [theta phi k2 k3 A]
p = p0(:)'; lam = 1e-3;
r = resid_of(p, dv, g, Y); f = sum(r.^2);
for it = 1:12
    J = zeros(numel(r), 5);
    for q = 1:5
        h = max(1e-4, abs(p(q))*1e-4); pq = p; pq(q) = pq(q)+h;
        J(:,q) = (resid_of(pq, dv, g, Y) - r)/h;
    end
    H = J'*J; grad = J'*r;
    for tryi = 1:6
        dg = diag(H); dg = max(dg, max(dg)*1e-8 + 1e-12);
        dp = -(H + lam*diag(dg) + 1e-9*(1+max(dg))*eye(5)) \ grad;
        pn = p + dp';
        pn(5) = max(pn(5), 1e-4);
        rn_ = resid_of(pn, dv, g, Y); fn = sum(rn_.^2);
        if fn < f
            p = pn; r = rn_; f = fn; lam = max(lam*0.4, 1e-9); break;
        else
            lam = lam*4;
        end
    end
    if norm(dp) < 1e-7, break; end
end
rn = sqrt(f/numel(r));
end

function r = resid_of(p, dv, g, Y)
[e1, kap] = unpack(p); A = p(5);
s   = dv * e1(:);                      % arc offset along the fibre, in voxels
U   = e1 + s .* kap;                   % transported direction per neighbour
U   = U ./ max(sqrt(sum(U.^2,2)),1e-12);
P   = (U * g').^2;                     % [N x G]
Yh  = -A * (P - mean(P,2));
r   = reshape(Yh - Y, [], 1);
end

function [e1, kap] = unpack(p)
e1 = ang2dir(p(1), p(2));
[e2, e3] = perp_basis(e1);
kap = p(3)*e2 + p(4)*e3;               % curvature vector, orthogonal to e1
end

function a = dir2ang(v)
a = [acos(max(min(v(3),1),-1)), atan2(v(2), v(1))];
end
function v = ang2dir(th, ph)
v = [sin(th)*cos(ph), sin(th)*sin(ph), cos(th)];
end
function [e2, e3] = perp_basis(e1)
a = [1 0 0]; if abs(e1(1)) > 0.9, a = [0 1 0]; end
e2 = a - dot(a,e1)*e1; e2 = e2/max(norm(e2),1e-12);
e3 = cross(e1, e2);
end
function v = getfielddef(s, f, d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% ---------------------------------------------------------------------------
function p = fit_frame_cont(p0, dv, g, Y, En, w)
% As fit_frame, plus the FIBRE CONTINUITY penalty: the direction this voxel's
% connection predicts at neighbour d must agree with that neighbour's own
% estimated direction. Sign-invariant (a line field), so the residual is
% sqrt(1 - (u.e)^2) rather than |u - e|.
p = p0(:)'; lam = 1e-3; dp = zeros(1,5);
% BALANCE THE TWO BLOCKS. The data block has (#neighbours x #gradients)
% residuals and the continuity block only (#neighbours) -- 864 vs 27 here. Left
% unnormalised the data term outweighs continuity 32:1 and the constraint is
% inert (measured: neighbour agreement 15.7 -> 14.3 deg for 5.8x the cost).
% Divide each block by sqrt(its count) so both contribute equally at w = 1.
nd = numel(Y); nc = size(En,1);
sd = 1/sqrt(max(nd,1)); sc = 1/sqrt(max(nc,1));
r = [sd*resid_of(p,dv,g,Y); w*sc*cont_resid(p,dv,En)];
f = sum(r.^2);
for it = 1:8
    J = zeros(numel(r),5);
    for q = 1:5
        h = max(1e-4, abs(p(q))*1e-4); pq = p; pq(q) = pq(q)+h;
        rq = [sd*resid_of(pq,dv,g,Y); w*sc*cont_resid(pq,dv,En)];
        J(:,q) = (rq - r)/h;
    end
    H = J'*J; grad = J'*r;
    for tryi = 1:6
        dg = diag(H); dg = max(dg, max(dg)*1e-8 + 1e-12);
        dp = -((H + lam*diag(dg) + 1e-9*(1+max(dg))*eye(5)) \ grad)';
        pn = p + dp; pn(5) = max(pn(5),1e-4);
        rn = [sd*resid_of(pn,dv,g,Y); w*sc*cont_resid(pn,dv,En)];
        fn = sum(rn.^2);
        if fn < f, p = pn; r = rn; f = fn; lam = max(lam*0.4,1e-9); break;
        else, lam = lam*4; end
    end
    if norm(dp) < 1e-7, break; end
end
end

function rc = cont_resid(p, dv, En)
[e1, kap] = unpack(p);
s = dv * e1(:);
U = e1 + s .* kap;
U = U ./ max(sqrt(sum(U.^2,2)),1e-12);
c = sum(U .* En, 2);
nn = sqrt(sum(En.^2,2));
ok = nn > 1e-9;
rc = zeros(size(c));
rc(ok) = sqrt(max(1 - (c(ok)./nn(ok)).^2, 0));   % sin of the disagreement angle
end
