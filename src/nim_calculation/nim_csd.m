function nim = nim_csd(nim, options)
% nim_csd: Constrained Spherical Deconvolution (Tournier 2007), single-shell,
% implemented natively in MATLAB. Produces a fiber orientation distribution (FOD)
% per voxel and extracts peaks into the SAME interface the multi-directional
% tracker consumes (nim.peaks / nim.npeaks), so CSD and the DTI-field
% reconstruction are interchangeable upstream of tractography.
%
% Pipeline (single-shell, single-tissue):
%   1. real SH fit of the DW signal (design matrix from bvecs, lmax)
%   2. single-fiber response function from high-FA/high-cl voxels
%   3. constrained spherical deconvolution (iterative non-negativity, Tournier07)
%   4. peak extraction from the FOD
%
% Arguments:
%   nim     - struct with .img (or .img_bi + .img_b0), .bval, .bvec, .evec,
%             .eval, .FA, .mask.
%   options - struct; uses:
%       .lmax          SH order (default 6; capped by #directions)
%       .sf_fa         FA threshold for single-fiber response voxels (0.7)
%       .fod_lambda    non-negativity regularization weight (1.0)
%       .fod_thresh    FOD amplitude treated as "negative" (0.1 * mean)
%       .n_iter        CSD iterations (50)
%       .peak_thresh   min peak amplitude, fraction of max FOD (0.1)
%       .peak_min_sep  min angular peak separation, degrees (25)
%       .max_peaks     max peaks kept per voxel (3)
%       .response      [optional] precomputed zonal response r_l (row) to skip est.
%
% Returns nim with .peaks, .npeaks, .peak_w, .fod_sh (SH FOD coeffs), .response.

if nargin < 2, options = struct(); end
o = @(f,d) getdef(options,f,d);
sf_fa      = o('sf_fa',0.7);
fod_lambda = o('fod_lambda',1.0);
n_iter     = o('n_iter',50);
peak_thresh= o('peak_thresh',0.1);
peak_sep   = o('peak_min_sep',25);
max_peaks  = o('max_peaks',3);

dims = size(nim.FA);

% ---- gradient table & DW signal --------------------------------------------
bval = nim.bval(:);
bvec = nim.bvec;                          % [Nvol x 3]
dw = bval > max(bval)*0.5;                % DW volumes (single shell)
b0 = bval <= max(bval)*0.1;
g = bvec(dw,:);                           % DW directions
g = g ./ max(sqrt(sum(g.^2,2)),1e-9);
ndir = size(g,1);

lmax = o('lmax', 6);
while (lmax+1)*(lmax+2)/2 > ndir && lmax > 2, lmax = lmax - 2; end   % cap by #dirs
fprintf('CSD: %d DW dirs, lmax=%d (%d SH coeffs)\n', ndir, lmax, (lmax+1)*(lmax+2)/2);

% signal attenuation S/S0
img = double(nim.img);
S0 = mean(img(:,:,:,b0), 4); S0(S0<=0) = NaN;
Sdw = img(:,:,:,dw);
A = Sdw ./ S0;                            % [D1 D2 D3 ndir] attenuation

% ---- real SH design matrix for the DW directions ---------------------------
B = real_sh_basis(g, lmax);              % [ndir x ncoef]
ncoef = size(B,2);
% least-squares SH fit operator (with light Tikhonov for stability)
Binv = (B'*B + 1e-3*eye(ncoef)) \ B';

% fit signal SH per voxel (vectorized over voxels)
if isfield(nim,'mask') && ~isempty(nim.mask), brain = nim.mask>0.5; else, brain = nim.FA>0.05; end
Avox = reshape(A, [], ndir)';            % [ndir x Nvox]
bad = any(isnan(Avox),1) | ~brain(:)';
Avox(:,bad) = 0;
Ssh = Binv * Avox;                        % [ncoef x Nvox] signal SH

% ---- response function (single-fiber) --------------------------------------
if ~isempty(o('response',[]))
    r = o('response',[]); r = r(:)';
else
    r = estimate_response(nim, A, g, lmax, sf_fa, brain);
end
fprintf('CSD: response zonal r_l = [%s]\n', num2str(r, '%.3g '));

% ---- forward convolution operator C (FOD SH -> signal SH) -------------------
% Axially-symmetric kernel: per order l, s_lm = f_lm * sqrt(4pi/(2l+1)) * r_l.
cscale = zeros(1, ncoef);
ci = 1;
for l = 0:2:lmax
    fac = sqrt(4*pi/(2*l+1)) * r(l/2+1);
    for m = -l:l, cscale(ci) = fac; ci = ci + 1; end
end
C = diag(cscale);                         % signal_sh = C * fod_sh

% ---- constrained spherical deconvolution -----------------------------------
% Directions for the non-negativity constraint (dense, ~300 on the hemisphere).
Uc = sphere_dirs(300);
Bc = real_sh_basis(Uc, lmax);            % [Ncon x ncoef] FOD-amp evaluator
tau = 0.1;                                % negativity threshold (fraction of mean amp)

CtC = C'*C;
Cts = C' * Ssh;                           % [ncoef x Nvox]
% initial FOD: unconstrained deconvolution
Fsh = (CtC + 1e-2*eye(ncoef)) \ Cts;

vidx = find(brain(:) & ~bad');
fprintf('CSD: deconvolving %d voxels (%d iters)...\n', numel(vidx), n_iter);
L2fixed = fod_lambda^2;
for it = 1:n_iter
    amp = Bc * Fsh(:, vidx);              % [Ncon x Nv] FOD amplitudes
    meanamp = mean(max(amp,0), 1) + 1e-9;
    for jj = 1:numel(vidx)
        v = vidx(jj);
        a = amp(:, jj);
        neg = a < tau * meanamp(jj);
        if ~any(neg), continue; end
        Ln = Bc(neg, :);
        M = CtC + L2fixed * (Ln'*Ln);
        Fsh(:, v) = M \ Cts(:, v);
    end
end

% ---- peak extraction --------------------------------------------------------
Up = sphere_dirs(600);                    % finer grid for peaks
Bp = real_sh_basis(Up, lmax);
peaks  = zeros([dims, max_peaks, 3]);
npeaks = zeros(dims);
peak_w = zeros([dims, max_peaks]);
for jj = 1:numel(vidx)
    v = vidx(jj);
    amp = Bp * Fsh(:, v);
    [pk_dirs, pk_amp] = find_peaks(Up, amp, peak_thresh, peak_sep, max_peaks);
    m = size(pk_dirs,1);
    [x,y,z] = ind2sub(dims, v);
    for k = 1:m
        peaks(x,y,z,k,:) = pk_dirs(k,:);
        peak_w(x,y,z,k)  = pk_amp(k);
    end
    npeaks(x,y,z) = m;
end

nim.fod_sh = reshape(Fsh', [dims, ncoef]);
nim.response = r;
nim.peaks = peaks; nim.npeaks = npeaks; nim.peak_w = peak_w;
fprintf('CSD: multi-peak voxels: %d (%.1f%% of brain)\n', nnz(npeaks>=2), 100*nnz(npeaks>=2)/nnz(brain));
end


% ===========================================================================
function r = estimate_response(nim, A, g, lmax, sf_fa, brain)
% Single-fiber response as zonal harmonics r_l (l=0,2,..,lmax): for each
% single-fiber voxel, rotate its DW directions into the fiber frame (e1 -> +z),
% fit the attenuation to the real SH basis, and take the zonal (m=0) coefficients
% -- which are the only non-zero ones for an axially symmetric response. Average
% over voxels. This is the correct, orthonormal projection (my earlier Legendre-
% on-|cos| fit was non-orthogonal and mis-scaled r_2 vs r_0).
l1 = nim.eval(:,:,:,1); l2 = nim.eval(:,:,:,2); l3 = nim.eval(:,:,:,3);
S = l1+l2+l3; Sg=S; Sg(Sg<=0)=1; cl = (l1-l2)./Sg;
sf = brain & (nim.FA > sf_fa) & (cl > 0.4);
idx = find(sf);
if numel(idx) < 20, sf = brain & (nim.FA > sf_fa-0.15); idx = find(sf); end
dims = size(nim.FA);
e1 = nim.evec(:,:,:,:,1);

nl = lmax/2 + 1;
zidx = zeros(1,nl); col=0; li=0;
for l = 0:2:lmax, li=li+1; zidx(li) = col + (l+1); col = col + (2*l+1); end  % m=0 positions
acc = zeros(1, nl); cnt = 0;
for ii = 1:numel(idx)
    [x,y,z] = ind2sub(dims, idx(ii));
    f = squeeze(e1(x,y,z,:))'; nf=norm(f); if nf<1e-6, continue; end; f=f/nf;
    a = squeeze(A(x,y,z,:)); if any(isnan(a)), continue; end
    R = rot_to_z(f);                       % rotates f -> +z
    glocal = (R * g')';                    % DW dirs in fiber frame
    Bl = real_sh_basis(glocal, lmax);
    c = (Bl'*Bl + 1e-3*eye(size(Bl,2))) \ (Bl' * a);
    acc = acc + c(zidx)'; cnt = cnt + 1;
end
r = acc / max(cnt,1);
end

function R = rot_to_z(f)
% rotation matrix taking unit vector f to +z.
z = [0 0 1];
v = cross(f, z); s = norm(v); c = dot(f, z);
if s < 1e-8
    if c > 0, R = eye(3); else, R = diag([1 -1 -1]); end
    return;
end
vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
R = eye(3) + vx + vx*vx*((1-c)/s^2);
end

% ===========================================================================
function B = real_sh_basis(dirs, lmax)
% Real, symmetric (even-order) SH basis, MRtrix ordering: for each even l,
% m = -l..l; m<0 -> sqrt2*Im(Ylm), m=0 -> Yl0, m>0 -> sqrt2*Re(Ylm).
n = size(dirs,1);
[phi, theta] = cart2sphr(dirs);           % theta=polar from +z, phi azimuth
ncoef = 0; for l=0:2:lmax, ncoef = ncoef + (2*l+1); end
B = zeros(n, ncoef);
col = 1;
for l = 0:2:lmax
    Pl = legendre(l, cos(theta(:)'));     % [(l+1) x n], rows m=0..l
    for m = -l:l
        am = abs(m);
        % normalization for real SH
        Nlm = sqrt((2*l+1)/(4*pi) * factorial(l-am)/factorial(l+am));
        P = Pl(am+1, :)';
        if m < 0
            B(:,col) = sqrt(2) * Nlm * P .* sin(am*phi(:));
        elseif m == 0
            B(:,col) = Nlm * P;
        else
            B(:,col) = sqrt(2) * Nlm * P .* cos(am*phi(:));
        end
        col = col + 1;
    end
end
end

function [phi, theta] = cart2sphr(d)
d = d ./ max(sqrt(sum(d.^2,2)),1e-9);
theta = acos(max(-1,min(1,d(:,3))));
phi = atan2(d(:,2), d(:,1));
end

function U = sphere_dirs(n)
% ~n roughly-uniform directions on the sphere (Fibonacci lattice).
k = (0:n-1)';
ga = pi*(3-sqrt(5));
z = 1 - 2*(k+0.5)/n;
rho = sqrt(max(0,1-z.^2));
ph = ga*k;
U = [rho.*cos(ph), rho.*sin(ph), z];
end

function [dirs, amps] = find_peaks(U, amp, thr, sep_deg, maxp)
% local maxima of amp over directions U (undirected: fold to a hemisphere).
mx = max(amp);
if mx <= 0, dirs = zeros(0,3); amps = []; return; end
cand = find(amp > thr*mx);
[~, ord] = sort(amp(cand), 'descend');
cand = cand(ord);
dirs = zeros(0,3); amps = [];
cs = cos(deg2rad(sep_deg));
for i = 1:numel(cand)
    d = U(cand(i), :);
    keep = true;
    for j = 1:size(dirs,1)
        if abs(d * dirs(j,:)') > cs, keep = false; break; end   % too close (undirected)
    end
    if keep
        dirs(end+1,:) = d; amps(end+1,1) = amp(cand(i)); %#ok<AGROW>
        if size(dirs,1) >= maxp, break; end
    end
end
end

function val = getdef(s,f,d), if isfield(s,f)&&~isempty(s.(f)), val=s.(f); else, val=d; end, end
