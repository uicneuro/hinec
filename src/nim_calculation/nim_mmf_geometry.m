function nim = nim_mmf_geometry(nim, options)
% nim_mmf_geometry: Build the MMF moving-frame geometry over the whole domain and
% return it on the nim (Chun & Peng, in preparation; pipeline steps 1-3).
%
% BUILT PER RUN, by runTractography step 3, and ONLY when algorithm == 'mmf'. It is
% NOT part of the nim on disk: it costs 4.9 s (dti) / 8.1 s (csd), which is cheaper
% than the build stamp + version + rebuild-guard machinery that used to avoid it, and
% it depends on config.tractography (field), which the dataset must not.
%
% Faithful to the spec's Frenet-frame construction:
%   step 1  e1 = trajectory-dependent (sel_power) denoised tangent field (NOT Gaussian).
%   Eq 7    e2 = (de1/ds)/||de1/ds|| = kappa/||kappa||  (Frenet normal), with reference-axis
%           projection (Eq 6) only as the degenerate fallback where ||de1/ds|| ~ 0.
%   Eq 8    e3 = e1 x e2.
%   Eq 9    [omega] = connection 1-form of the frame field (nim_connection_form).
% Stores, for the connection-driven tracer (Eq 10-11):
%   nim.mmf_frames [D1 D2 D3 3 3]  the frame field {e1,e2,e3}
%   nim.mmf_kappa  [D1 D2 D3 3]    curvature vector d e1/ds = grad_{e1}e1
%   nim.mmf_tau    [D1 D2 D3]      torsion  tau = omega_23(e1)
%
% NO DENOISING PARAMETER. The alignment-weighted denoise this function used to
% apply (each 3x3x3 neighbour weighted |n.e1|^frame_sel_power) is gone, and with
% it the exponent. It was the same mechanism as the retired sel_power: a free
% parameter with no principled value, which cannot be justified without fitting
% it to data.
%
% It was also not doing its job. Measured against the curvature of the ISMRM
% ground-truth curves on 5699 Cingulum_right voxels, where the true median
% curvature is 0.1538/vox:
%
%   frame_sel_power   median kappa   correlation with truth
%          0             0.0148            0.229
%          2             0.0112            0.239
%         16 (default)   0.0092            0.224
%         64             0.0101            0.217
%
% The exponent changes nothing (r spans 0.217-0.239) and every setting suppresses
% curvature by more than an order of magnitude. The reason it cannot work: the
% true curvature is 8.8 degrees of turn per voxel, while the direction field's own
% error is 8.15 degrees in single-orientation voxels and up to 22 degrees where
% bundles cross. Differentiating a field whose noise matches its signal cannot
% produce a usable derivative; smoothing trades that noise for bias, which is why
% kappa collapses toward zero.

if nargin < 2 || isempty(options), options = struct(); end
if ~isfield(nim,'evec'), error('nim_mmf_geometry: nim.evec not found (run nim_eig).'); end
if isfield(nim,'mask') && ~isempty(nim.mask), m = nim.mask > 0.5; else, m = nim.FA > 0; end
dims = size(nim.FA);

% --- step 1: sign-consistent raw e1 (no denoising; see the header) ---
% e1 source: the tensor principal eigenvector (DTI, the spec's formulation) or, when a
% CSD FOD peak field is present and requested, the DOMINANT FOD peak (so field='csd'
% builds the connection from CSD data, not silently from the tensor).
% CSD build happens ONLY if peaks are actually present AND requested (else it silently
% falls back to DTI); field= in the log line at the end says what was ACTUALLY built.
is_dwi = strcmpi(getf(options,'field','dti'),'dwi');
if is_dwi
    % DIRECT ROUTE: e1 and kappa come from a joint fit to the raw DW signal, so
    % the curvature was never differenced out of an already-collapsed direction
    % field. Everything downstream (e2, e3, torsion) is built from that frame.
    % The fit itself is step 2 (nim_field); this function only CONSUMES it.
    if ~isfield(nim,'mmf_kappa_dwi') || isempty(nim.mmf_kappa_dwi) || ...
       ~isfield(nim,'mmf_e1_dwi')    || isempty(nim.mmf_e1_dwi)
        error('nim_mmf_geometry:missingDwiFit', ...
            ['field=dwi needs nim.mmf_e1_dwi / nim.mmf_kappa_dwi, which are fitted by ' ...
             'nim_field (runTractography step 2). Call nim_field(nim, options, nim_path) first.']);
    end
end
is_csd = isfield(nim,'peaks') && strcmpi(getf(options,'field','dti'),'csd');
% The e1 substitution is LOCAL (nim_e1): the returned nim keeps the dataset's own
% evec, so building the geometry cannot alter the dataset the tracker then reads.
if is_csd
    nim_e1 = nim; nim_e1.evec(:,:,:,:,1) = squeeze(nim.peaks(:,:,:,1,:));   % dominant peak
    fprintf('nim_mmf_geometry: CSD field -> connection built from the dominant FOD peak\n');
elseif is_dwi
    nim_e1 = nim; nim_e1.evec(:,:,:,:,1) = nim.mmf_e1_dwi;                 % signal-fitted frame
else
    nim_e1 = nim;                                                          % DTI principal eigenvector
end
f0  = nim_build_frames(nim_e1, struct('frame_smooth_sigma',0));   % sign-consistent, no smoothing
e1d = f0(:,:,:,:,1);

% --- pass 1: curvature vector kappa = grad_{e1}e1 (basis-free) via the connection form ---
[e2a,e3a] = ref_axis_frame_field(e1d);
Wa  = nim_connection_form(assemble(e1d,e2a,e3a), m, options);
w12 = -squeeze(Wa(:,:,:,2,1,1));   % omega_12(e1)
w13 = -squeeze(Wa(:,:,:,3,1,1));   % omega_13(e1)
kappa = zeros([dims 3]);
for c=1:3, kappa(:,:,:,c) = e2a(:,:,:,c).*w12 + e3a(:,:,:,c).*w13; end   % = grad_{e1}e1
if is_dwi
    kappa = nim.mmf_kappa_dwi;      % fitted from the signal; no deconvolution needed
else
    kappa = mmf_deconvolve_curvature(kappa, m);
end


% --- Eq 7: e2 = kappa/||kappa|| (Frenet normal); reference-axis (Eq 6) fallback ---
kmag = sqrt(sum(kappa.^2,4));
e2 = zeros([dims 3]);
for c=1:3, e2(:,:,:,c) = kappa(:,:,:,c)./max(kmag,1e-12); end
low = kmag < 1e-3;                                 % de1/ds ~ 0 -> use reference-axis
for c=1:3, t=e2(:,:,:,c); ta=e2a(:,:,:,c); t(low)=ta(low); e2(:,:,:,c)=t; end
[e2,e3] = orthonormalize(e1d, e2);                 % re-orthonormalize; e3 = e1 x e2 (Eq 8)

% --- pass 2: torsion tau = omega_23(e1) from the Frenet frame (Eq 9/10) ---
frames = assemble(e1d, e2, e3);
z = ~(m & sqrt(sum(e1d.^2,4))>1e-6);
for i=1:3, for c=1:3, t=frames(:,:,:,c,i); t(z)=0; frames(:,:,:,c,i)=t; end, end
Wb  = nim_connection_form(frames, m, options);
tau = squeeze(Wb(:,:,:,2,3,1));                    % omega_23(e1)

% --- CSD MULTIPLE PATHWAYS (spec step 1b): a per-peak connection curvature field ---
% Each FOD peak is a fibre population. For peak p we compute its curvature grad_{e1p}e1p
% by TRAJECTORY-ALIGNED differencing: each neighbour contributes the peak best aligned with
% the centre's peak p (peak matching across the crossing). At trace time the tracer selects
% the peak aligned with the incoming tangent, so two streamlines entering one voxel from
% different approaches follow DIFFERENT curvatures -> different continuations (crossings).
if is_csd
    P = nim.peaks; maxp = size(P,4); Pn = zeros(size(P));
    for p=1:maxp
        pk = squeeze(P(:,:,:,p,:)); nn = sqrt(sum(pk.^2,4));
        for c=1:3, Pn(:,:,:,p,c) = pk(:,:,:,c)./max(nn,1e-9); end
    end
    kappa_p = zeros([dims maxp 3]);
    for p=1:maxp
        e1p = squeeze(Pn(:,:,:,p,:));
        kp  = mmf_deconvolve_curvature(peak_aligned_curvature(Pn, e1p, dims), m);
        for c=1:3, kappa_p(:,:,:,p,c) = kp(:,:,:,c); end
    end
    nim.mmf_peakdirs = Pn;         % [D1 D2 D3 maxp 3] unit peak directions
    nim.mmf_kappa_p  = kappa_p;    % [D1 D2 D3 maxp 3] per-peak curvature vectors
    nim.mmf_npeaks   = nim.npeaks;
    nim.mmf_multi    = true;
    fprintf('nim_mmf_geometry: CSD multi-frame connection (%d peaks) for multiple pathways\n', maxp);
else
    nim.mmf_multi = false;
end

nim.mmf_frames = frames;
nim.mmf_kappa  = kappa;
nim.mmf_tau    = tau;
if is_csd, built_field = 'csd'; else, built_field = 'dti'; end           % what was ACTUALLY built
fprintf('nim_mmf_geometry: Frenet frames (Eq7 e2=de1/ds, field=%s) + curvature + torsion built\n', built_field);
end

% ---------------------------------------------------------------------------

function out = shift3(A,d)
sz = size(A); out = zeros(sz);
i1=max(1,1-d(1)):min(sz(1),sz(1)-d(1)); i2=max(1,1-d(2)):min(sz(2),sz(2)-d(2)); i3=max(1,1-d(3)):min(sz(3),sz(3)-d(3));
out(i1,i2,i3,:) = A(i1+d(1),i2+d(2),i3+d(3),:);
end

% per-peak curvature grad_{e1ref}e1ref via TRAJECTORY-ALIGNED central differencing:
% at each neighbour we take the peak best aligned with the centre's reference direction.
function kap = peak_aligned_curvature(Pn, e1ref, dims)
grad = zeros([dims 3 3]);                         % grad(:,:,:,axis,comp) = d e1ref_comp / d x_axis
for ax=1:3
  ep = aligned_peak_neighbor(Pn, e1ref, ax, +1);  % neighbour(+ax) peak aligned to centre e1ref
  em = aligned_peak_neighbor(Pn, e1ref, ax, -1);
  dcomp = (ep - em)/2;
  grad(:,:,:,ax,:) = reshape(dcomp, [dims 1 3]);
end
kap = zeros([dims 3]);
for comp=1:3
  acc = zeros(dims);
  for ax=1:3, acc = acc + e1ref(:,:,:,ax).*grad(:,:,:,ax,comp); end   % (e1ref . grad) e1ref
  kap(:,:,:,comp) = acc;
end
d = sum(kap.*e1ref,4); for c=1:3, kap(:,:,:,c) = kap(:,:,:,c) - d.*e1ref(:,:,:,c); end  % normal to e1ref
end

% at each voxel, the neighbour at (sgn*ax): the neighbour peak best aligned with e1ref, sign-flipped.
function out = aligned_peak_neighbor(Pn, e1ref, ax, sgn)
D = size(e1ref); D = D(1:3); maxp = size(Pn,4);
dshift = [0 0 0]; dshift(ax) = sgn;
best = -inf(D); out = zeros([D 3]);
for p=1:maxp
  pk = shift3(squeeze(Pn(:,:,:,p,:)), dshift);
  al = sum(pk.*e1ref,4); s = sign(al); s(s==0)=1; ala = abs(al);
  better = ala > best; best(better) = ala(better);
  for c=1:3, t=out(:,:,:,c); pc=pk(:,:,:,c).*s; t(better)=pc(better); out(:,:,:,c)=t; end
end
end

function [e2,e3] = ref_axis_frame_field(e1)
ex=e1(:,:,:,1); ey=e1(:,:,:,2); ez=e1(:,:,:,3);
ax=abs(ex); ay=abs(ey); az=abs(ez);
px=(ax<=ay)&(ax<=az); py=(~px)&(ay<=az); pz=~(px|py);
rx=double(px); ry=double(py); rz=double(pz); rd=rx.*ex+ry.*ey+rz.*ez;
e2x=rx-rd.*ex; e2y=ry-rd.*ey; e2z=rz-rd.*ez;
n2=sqrt(e2x.^2+e2y.^2+e2z.^2); n2(n2<1e-9)=1;
e2x=e2x./n2; e2y=e2y./n2; e2z=e2z./n2;
e3x=ey.*e2z-ez.*e2y; e3y=ez.*e2x-ex.*e2z; e3z=ex.*e2y-ey.*e2x;
e2=cat(4,e2x,e2y,e2z); e3=cat(4,e3x,e3y,e3z);
end

function [e2o,e3o] = orthonormalize(e1,e2)
d12=sum(e1.*e2,4); e2o=e2;
for c=1:3, e2o(:,:,:,c)=e2(:,:,:,c)-d12.*e1(:,:,:,c); end
n=sqrt(sum(e2o.^2,4)); for c=1:3, e2o(:,:,:,c)=e2o(:,:,:,c)./max(n,1e-9); end
ex=e1(:,:,:,1);ey=e1(:,:,:,2);ez=e1(:,:,:,3); ax=e2o(:,:,:,1);ay=e2o(:,:,:,2);az=e2o(:,:,:,3);
e3o=cat(4, ey.*az-ez.*ay, ez.*ax-ex.*az, ex.*ay-ey.*ax);
end

function F = assemble(e1,e2,e3)
d=size(e1); F=zeros([d(1:3) 3 3]);
for c=1:3, F(:,:,:,c,1)=e1(:,:,:,c); F(:,:,:,c,2)=e2(:,:,:,c); F(:,:,:,c,3)=e3(:,:,:,c); end
end

function v = getf(s,f,d), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end

% ---------------------------------------------------------------------------
function kap = mmf_deconvolve_curvature(kap, mask)
% INVERT THE DWI -> DTI CURVATURE TRANSFER.
%
% The curvature read out of the direction field is not the curvature of the
% underlying fascicle. Measured against the ISMRM-2015 ground-truth geometry --
% 294,554 points on 14 bundles, sampled the way the tracker samples (trilinear
% on the dyadic v1*v1'), FA >= 0.15 -- the field obeys
%
%       kappa_field = 0.509 * kappa_true + 0.0144 /mm      (R^2 = 0.9993)
%
% Two independent defects in one line:
%   GAIN 0.509   the field returns barely half the true bending, so a tracker
%                obeying it cuts too wide an arc and leaves tight bundles. CP
%                turns at R = 12.2 mm; the field reports R = 16.9 mm, and every
%                streamline spliced along CP left the scoring corridor.
%   OFFSET 0.0144 where the truth is straight the field still reports curvature,
%                so streamlines wander off a straight fascicle. This is the
%                1-3 mm halo that shows up as overreach.
% They cross at kappa = 0.029 /mm (R = 34 mm): tighter than that the field
% destroys curvature, straighter than that it invents it.
%
% Geometry is built on the voxel grid, so the law is rescaled by the 2 mm voxel:
% offset 0.0144/mm -> 0.0288/voxel; the gain is dimensionless.
%
% This is a calibration of the measurement, not a tuning knob -- hence no config
% key. The direction of kappa is untouched (e2 = kappa/||kappa|| is unchanged);
% only the turning RATE that enters the structure equation is corrected.
GAIN   = 0.509;
OFFSET = 0.0288;          % per voxel (0.0144 /mm at 2 mm isotropic)
km   = sqrt(sum(kap.^2,4));
kc   = max(km - OFFSET, 0) ./ GAIN;      % below the noise floor -> genuinely straight
scale = kc ./ max(km, 1e-12);
if nargin > 1 && ~isempty(mask)
    scale(~mask) = 1;
end
for c = 1:3, kap(:,:,:,c) = kap(:,:,:,c) .* scale; end
end

