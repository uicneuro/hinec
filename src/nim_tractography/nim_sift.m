function [tracks_out, info] = nim_sift(tracks, nim, options)
% nim_sift: SIFT-style filtering of a tractogram (Smith 2013 analog), with
% PER-LOBE attribution.
%
% Prunes streamlines so the streamline density on each FOD LOBE (fiber population)
% matches that lobe's FOD amplitude. Each streamline segment is attributed to the
% FOD peak it best aligns with; a lobe carrying more streamline-density than its
% amplitude warrants is over-represented, and the streamlines running through
% over-represented lobes are removed. This distinguishes a legitimately dense
% bundle (matched to a strong peak) from invalid spray (over-loading weak/wrong
% lobes) — which a per-voxel density match cannot.
%
% Requires nim.peaks [X Y Z maxK 3] and nim.peak_w [X Y Z maxK] (from CSD).
%
% options: sift_batch_frac (0.03), sift_n_iter (60), sift_min_keep (0.20),
%          sift_align_pow (weight excess by misalignment; 0=off, default 1).

if nargin < 3, options = struct(); end
o = @(f,d) getdef(options,f,d);
batch_frac = o('sift_batch_frac', 0.03);
n_iter     = o('sift_n_iter', 60);
min_keep   = o('sift_min_keep', 0.20);
align_pow  = o('sift_align_pow', 1);

dims = size(nim.FA); nv = prod(dims);
nt = numel(tracks);
assert(isfield(nim,'peaks') && isfield(nim,'peak_w'), 'nim_sift needs nim.peaks + nim.peak_w (run CSD)');
maxK = size(nim.peak_w, 4);
brain = true(nv,1); if isfield(nim,'mask')&&~isempty(nim.mask), brain = nim.mask(:)>0.5; end

PW = reshape(nim.peak_w, nv, maxK); PW(~brain,:) = 0;                 % target amplitudes
PK = reshape(nim.peaks,  nv, maxK, 3);                                % peak directions

% ---- flatten segments with unit directions --------------------------------
[segvox, seglen, segdir, segtrk] = flatten_dir(tracks, dims);
Ns = numel(segvox);
if Ns == 0, tracks_out = tracks; info = struct('n_in',nt,'n_out',nt,'keep_frac',1); return; end

% ---- attribute each segment to its best-aligned FOD lobe ------------------
best = -inf(Ns,1); seglobe = ones(Ns,1);
for kk = 1:maxK
    pks = squeeze(PK(:,kk,:)); pks = pks(segvox,:);                  % [Ns 3]
    npk = sqrt(sum(pks.^2,2));
    a = abs(sum(segdir.*pks, 2)) ./ max(npk,1e-9);
    a(npk < 1e-6) = -inf;
    b = a > best; seglobe(b) = kk; best(b) = a(b);
end
best(~isfinite(best)) = 0;
% "excess" weight: misaligned segments (best alignment low) count as more removable
segw = seglen .* (1 + align_pow*(1 - max(0,best)));

nvl  = nv*maxK;
segidx = (segvox-1)*maxK + seglobe;                                  % (voxel,lobe) index
TD = zeros(nvl,1);
for kk = 1:maxK, TD((0:nv-1)*maxK + kk) = PW(:,kk); end

% ---- greedy density-matching removal --------------------------------------
keep = true(nt,1); active = true(Ns,1);
SD = accumarray(segidx, segw, [nvl 1]);
k  = sum(SD)/max(sum(TD),eps);
prev = inf; cost_curve = zeros(n_iter,1);
for it = 1:n_iter
    M = SD - k*TD; tot = sum(abs(M)); cost_curve(it) = tot;
    if tot >= prev - 1e-9, cost_curve = cost_curve(1:it); break; end
    prev = tot;
    pos = max(M,0);
    tcost = accumarray(segtrk(active), segw(active).*pos(segidx(active)), [nt 1]);
    tcost(~keep) = -inf;
    cand = find(tcost > 0); if isempty(cand), break; end
    nrem = max(1, round(batch_frac*nt));
    if (sum(keep)-nrem) < min_keep*nt, nrem = max(0, sum(keep)-round(min_keep*nt)); end
    if nrem <= 0, break; end
    [~,ord] = sort(tcost(cand),'descend');
    keep(cand(ord(1:min(nrem,numel(cand))))) = false;
    active = keep(segtrk);
    SD = accumarray(segidx(active), segw(active), [nvl 1]);
    k  = sum(SD)/max(sum(TD),eps);
end

tracks_out = tracks(keep);
info = struct('n_in',nt,'n_out',numel(tracks_out), ...
              'keep_frac',numel(tracks_out)/max(nt,1),'cost_curve',cost_curve);
fprintf('SIFT(lobe): kept %d/%d (%.1f%%), mismatch %.3g -> %.3g\n', ...
        info.n_out, nt, 100*info.keep_frac, cost_curve(1), cost_curve(end));
end

% ---------------------------------------------------------------------------
function [segvox, seglen, segdir, segtrk] = flatten_dir(tracks, dims)
nt = numel(tracks); nseg = zeros(nt,1);
for i=1:nt, nseg(i)=max(0,size(tracks{i},1)-1); end
N = sum(nseg);
segvox=zeros(N,1); seglen=zeros(N,1); segdir=zeros(N,3); segtrk=zeros(N,1); p=0;
for i=1:nt
    t=tracks{i}; if size(t,1)<2, continue; end
    d=diff(t,1,1); L=sqrt(sum(d.^2,2)); mid=(t(1:end-1,:)+t(2:end,:))/2;
    v=round(mid);
    v(:,1)=min(max(v(:,1),1),dims(1)); v(:,2)=min(max(v(:,2),1),dims(2)); v(:,3)=min(max(v(:,3),1),dims(3));
    li=v(:,1)+(v(:,2)-1)*dims(1)+(v(:,3)-1)*dims(1)*dims(2);
    dir=d./max(L,1e-9);
    m=numel(L); rng=p+1:p+m; segvox(rng)=li; seglen(rng)=L; segdir(rng,:)=dir; segtrk(rng)=i; p=p+m;
end
end

function val = getdef(s,f,d), if isfield(s,f)&&~isempty(s.(f)), val=s.(f); else, val=d; end, end
