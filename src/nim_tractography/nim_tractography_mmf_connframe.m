function [tracks, info] = nim_tractography_mmf_connframe(nim, options)
% nim_tractography_mmf_connframe: GENUINE Method-of-Moving-Frames tractography
% (Chun & Peng, in preparation). PURE connection-form tracer -- the
% interpolated streamline tracker (sel_power / RK4 / RKF45 / CSD-peak resampling) lives
% in nim_tractography_hinec.m, NOT here.
%
% Pipeline: runTractography step 3 calls nim_mmf_geometry, which builds the moving-frame
% field {e1,e2,e3} + connection 1-form (curvature + torsion) for THIS run (Eq 6-9) and
% hands it in on the nim. This function traces by evolving the carried frame with the
% structure equation (Eq 10) while advancing dx/ds=e1 (Eq 11). It never builds the
% geometry itself. CSD field -> a per-peak connection gives MULTIPLE PATHWAYS
% (crossing resolution).
%
% options.field      : 'dti' (tensor principal direction) | 'csd' (per-peak, multi-frame)
% options.mmf_anchor : 0 = pure Eq.10-11 ; >0 re-anchors e1 toward the field tangent
% requires: nim.evec/FA/mask ; nim.peaks/npeaks for field='csd'.

o=@(f,d) getdef(options,f,d);
options.field=char(string(o('field','dti')));
options.step_size=o('step_size',0.2); options.step_min=o('step_min',0.02); options.step_max=o('step_max',0.5);
options.max_steps=o('max_steps',3000);
options.angle_thresh=o('angle_thresh',60); options.min_length=o('min_length',15);
options.fa_threshold=o('fa_threshold',0.1); options.termination_fa=o('termination_fa',0.05);
options.bishop_eps=o('bishop_eps',1e-6);
options.mmf_anchor=o('mmf_anchor',0);   % 0 = pure Eq.10-11 (faithful); >0 re-anchors e1 to the field
% integrator = the NUMERICAL STEPPING SCHEME for the connection-form ODE (orthogonal to the
% direction, which comes from the connection form): 'rk4' (fixed step) | 'rkf45' (adaptive
% Dormand-Prince, smaller steps in curved regions). rkf_tol is the RKF45 error tolerance.
options.integrator=char(string(o('integrator','rk4')));
options.rkf_tol=o('rkf_tol',0.02);
% interp_method = HOW the pre-built connection form (kappa/tau/e1) is sampled at
% continuous positions: 'cubic' for smoother sampling, else 'linear' (trilinear).
% Orthogonal to how the connection is BUILT. Peak DIRECTIONS always
% use 'nearest' regardless (sign-ambiguity would break under linear/cubic averaging).
options.interp_method=char(string(o('interp_method','trilinear')));
if strcmpi(options.interp_method,'cubic'), gm='cubic'; else, gm='linear'; end
options.seed_density=o('seed_density',1);
options.wm_mask=o('wm_mask',[]); options.gm_mask=o('gm_mask',[]); options.csf_mask=o('csf_mask',[]);
options.seed_mask=o('seed_mask',[]); options.propagation_mask=o('propagation_mask',[]);

dims=size(nim.FA); gv={1:dims(1),1:dims(2),1:dims(3)};
is_csd=strcmpi(options.field,'csd');
if is_csd
  assert(isfield(nim,'peaks'),'csd field needs nim.peaks/npeaks (run nim_csd)');
  nim.PK=nim.peaks; nim.NP=nim.npeaks;
else
  nim.E=squeeze(nim.evec(:,:,:,:,1));
end
nim.FA_i=griddedInterpolant(gv,nim.FA,gm,'none');

% ===== the moving-frame geometry (Eq 6-9) MUST arrive already built =================
assert(isfield(nim,'mmf_frames'), ...
  ['nim_tractography_mmf_connframe: the connection geometry is missing. It is built ' ...
   'per run by runTractography step 3 (nim = nim_mmf_geometry(nim, options)), not by ' ...
   'this tracker and not by main.m.']);
% ANCHORING IS SPATIALLY VARYING, driven by the tensor's own shape.
%
% A constant anchor forces one compromise everywhere, and the trade-off is
% monotone: on Cingulum_right, anchor 0 gives the best crossing survival (39%
% of streamlines complete the turn) and the worst bundle recall (30%), while
% anchor 0.25 gives the best recall measured anywhere (51%) and a turn rate
% below plain DTI (23%). Neither setting is right in both places, because the
% two places are different: where one fibre population dominates the data
% should be trusted, and where two cross it cannot be.
%
% The tensor says which is which. Linearity CL=(l1-l2)/l1 and planarity
% CP=(l2-l3)/l1 separate the derailment hotspots cleanly - CP is 0.183 there
% against 0.038 over the rest of the bundle, at equal CL. So the anchor is
% scaled by the linear fraction
%
%     a_eff(x) = mmf.anchor * CL / (CL + CP)
%
% which is ~a where the fibre is unambiguous and decays toward 0 where the
% tensor goes planar, letting the connection form carry the curve through the
% crossing. NO NEW PARAMETER: mmf.anchor is now the anchoring applied where the
% data is trustworthy, and the modulation comes from the data itself.
L = sort(reshape(nim.eval, [], 3), 2, 'descend');
CLv = (L(:,1)-L(:,2)) ./ max(L(:,1), eps);
CPv = (L(:,2)-L(:,3)) ./ max(L(:,1), eps);
amod = reshape(CLv ./ max(CLv + CPv, eps), size(nim.FA));
nim.MAmod = griddedInterpolant(gv, amod, gm, 'none');   % anchor modulation in [0,1]

nim.ME1x=griddedInterpolant(gv,nim.mmf_frames(:,:,:,1,1),gm,'none');   % field e1 (anchor)
nim.ME1y=griddedInterpolant(gv,nim.mmf_frames(:,:,:,2,1),gm,'none');
nim.ME1z=griddedInterpolant(gv,nim.mmf_frames(:,:,:,3,1),gm,'none');
nim.MKx =griddedInterpolant(gv,nim.mmf_kappa(:,:,:,1),gm,'none');       % curvature vector
nim.MKy =griddedInterpolant(gv,nim.mmf_kappa(:,:,:,2),gm,'none');
nim.MKz =griddedInterpolant(gv,nim.mmf_kappa(:,:,:,3),gm,'none');
if isfield(nim,'mmf_tau'), nim.MTau=griddedInterpolant(gv,nim.mmf_tau,gm,'none');  % torsion
else, nim.MTau=griddedInterpolant(gv,zeros(dims),gm,'none'); end
% multi-frame (CSD): per-peak direction + curvature interpolants for multiple pathways
nim.mmf_multi = isfield(nim,'mmf_multi') && nim.mmf_multi;
if nim.mmf_multi
  mp=size(nim.mmf_peakdirs,4); nim.MPk=cell(mp,1); nim.MKp=cell(mp,1);
  for p=1:mp
    % peak DIRECTION: nearest-neighbour (sign-ambiguous +/-; linear interp across a
    % sign-flip seam cancels and would spuriously drop the peak).
    nim.MPk{p}={griddedInterpolant(gv,nim.mmf_peakdirs(:,:,:,p,1),'nearest','none'),...
                griddedInterpolant(gv,nim.mmf_peakdirs(:,:,:,p,2),'nearest','none'),...
                griddedInterpolant(gv,nim.mmf_peakdirs(:,:,:,p,3),'nearest','none')};
    % curvature vector is sign-INVARIANT -> linear is safe/smooth.
    nim.MKp{p}={griddedInterpolant(gv,nim.mmf_kappa_p(:,:,:,p,1),gm,'none'),...
                griddedInterpolant(gv,nim.mmf_kappa_p(:,:,:,p,2),gm,'none'),...
                griddedInterpolant(gv,nim.mmf_kappa_p(:,:,:,p,3),gm,'none')};
  end
  fprintf('MMF: %d-peak multi-frame connection; per-approach pathway selection\n', mp);
end
fprintf('MMF: connection-form geometry built this run (field=%s); connection-driven tracing (Eq.10-11)\n', options.field);

if isempty(options.seed_mask)
  sm=nim.FA>options.fa_threshold; if isfield(nim,'mask')&&~isempty(nim.mask), sm=sm&(nim.mask>0.5); end; options.seed_mask=sm;
end
options.seed_mask=logical(options.seed_mask>0);
if ~isempty(options.propagation_mask), prop=logical(options.propagation_mask>0); else, prop=nim.FA>options.termination_fa; if isfield(nim,'mask'), prop=prop&(nim.mask>0.5); end; end
nim.prop_mask=imdilate(prop,ones(3,3,3));

seeds=build_seeds(nim,options,dims,is_csd);
fprintf('MMF [%s]: %d seeds\n', options.field, size(seeds,1));
% Degrees of turning per VOXEL OF ARC (nim_angle_limit), the same units the
% other two trackers use. MMF advances dx/ds=e1 at a fixed step, so the arc a
% step covers is the nominal step size.
[~,cos_thresh]=nim_angle_limit(options.angle_thresh, options.step_size);
ns=size(seeds,1); allt=cell(ns,1); valid=false(ns,1);

use_par=~isempty(ver('parallel'));
if use_par
  pool=gcp('nocreate'); if isempty(pool), w=getenv('HINEC_MAX_WORKERS'); if ~isempty(w)&&~isnan(str2double(w)), nw=max(1,round(str2double(w))); else, nw=8; end; pool=parpool('local',nw); end
  fprintf('Using %d workers\n', pool.NumWorkers);
end
% WHICH SEEDS ARE TRACED. Sampled evenly across the seed list, not a prefix:
% seed order follows the voxel scan, so a prefix is one corner of the volume.
do_tr=false(ns,1);
if getdef(options,'trace',false)
  tmax=getdef(options,'trace_max',200);
  if tmax<=0||tmax>=ns, do_tr(:)=true; else, do_tr(round(linspace(1,ns,tmax)))=true; end
  fprintf('MMF: per-step tracing enabled for %d of %d seeds\n', nnz(do_tr), ns);
end
alltr=cell(ns,1);
% Why every streamline stopped, for EVERY seed -- not just traced ones. This is
% two strings per seed; without it a run reports nothing about termination and
% the low-overlap question ("do streamlines stop before joining the two endpoint
% regions?") cannot be answered from a finished run at all.
allterm=cell(ns,1);
t0=tic;
if use_par
  parfor i=1:ns
    [allt{i},alltr{i},allterm{i}]=track_bi(nim,seeds(i,1:3),seeds(i,4:6),options,cos_thresh,dims,do_tr(i));
    valid(i)=size(allt{i},1)>1;
  end
else
  for i=1:ns
    [allt{i},alltr{i},allterm{i}]=track_bi(nim,seeds(i,1:3),seeds(i,4:6),options,cos_thresh,dims,do_tr(i));
    valid(i)=size(allt{i},1)>1;
  end
end
fprintf('Tracking loop: %.1f s\n', toc(t0));
raw=allt(valid); tracks=cell(numel(raw),1); k=0;
for i=1:numel(raw), t=raw{i}; if sum(sqrt(sum(diff(t,1,1).^2,2)))>=options.min_length, k=k+1; tracks{k}=t; end, end
tracks=tracks(1:k);
info=struct('n_tracks',k,'n_seeds',ns,'field',options.field,'method','mmf-connection-form');
tt=[allterm{:}]; tt=tt(~cellfun(@isempty,tt));
if ~isempty(tt)
    [u,~,ii]=unique(tt); c=accumarray(ii,1); [c,o]=sort(c,'descend');
    info.termination = struct('reason',{u(o)},'count',{c});
    fprintf('MMF termination (%d half-streamlines):', numel(tt));
    for j=1:numel(c), fprintf(' %s=%.1f%%', u{o(j)}, 100*c(j)/numel(tt)); end
    fprintf('\n');
end
% The trace is indexed by SEED and is deliberately NOT subset to surviving
% tracks: the streamlines that died are exactly the ones worth diagnosing.
if any(do_tr)
    ti=find(do_tr(:))';
    info.trace = struct('seed_index', num2cell(ti), ...
                        'seed_point', num2cell(seeds(ti,1:3),2)', ...
                        'forward',  cellfun(@(c) c.forward,  alltr(ti), 'uni', 0)', ...
                        'backward', cellfun(@(c) c.backward, alltr(ti), 'uni', 0)');
end
fprintf('MMF [%s]: %d tracks\n', options.field, k);
end

% ---------------------------------------------------------------------------
function seeds=build_seeds(nim,options,dims,is_csd)
idx=find(options.seed_mask); [x,y,z]=ind2sub(dims,idx);
% exactly seed_density offsets (see nim_seed_offsets); the old inline lattice
% rounded up to the next perfect cube.
off=nim_seed_offsets(options.seed_density); no=size(off,1);
if is_csd, maxK=size(nim.PK,4); else, maxK=1; end
seeds=zeros(numel(idx)*no*maxK,6); r=0;
for i=1:numel(idx)
  if is_csd
    np=nim.NP(x(i),y(i),z(i));
    for p=1:np
      d=[nim.PK(x(i),y(i),z(i),p,1),nim.PK(x(i),y(i),z(i),p,2),nim.PK(x(i),y(i),z(i),p,3)]; nd=norm(d); if nd<1e-6, continue; end
      for s=1:no, r=r+1; seeds(r,:)=[[x(i),y(i),z(i)]+off(s,:), d/nd]; end
    end
  else
    d=squeeze(nim.E(x(i),y(i),z(i),:))'; if norm(d)<1e-6, continue; end; d=d/norm(d);
    for s=1:no, r=r+1; seeds(r,:)=[[x(i),y(i),z(i)]+off(s,:), d]; end
  end
end
seeds=seeds(1:r,:);
end

function [combined,trc,terms]=track_bi(nim,pos,dir0,options,cos_thresh,dims,do_trace)
if nargin<7||isempty(do_trace), do_trace=false; end
[tf,trf,ttf]=track_one(nim,pos, dir0,options,cos_thresh,dims,do_trace);
[tb,trb,ttb]=track_one(nim,pos,-dir0,options,cos_thresh,dims,do_trace);
trc=struct('forward',trf,'backward',trb);
terms={ttf,ttb};
if size(tb,1)>1, tb=flipud(tb(2:end,:)); else, tb=zeros(0,3); end
if size(tf,1)>1, tf=tf(2:end,:); else, tf=zeros(0,3); end
combined=[tb;pos;tf]; if size(combined,1)<=1, combined=zeros(0,3); end
end

function [track,trace,term]=track_one(nim,pos,v,options,cos_thresh,dims,do_trace)
track=zeros(options.max_steps+1,3); track(1,:)=pos; n=1; x=pos; v=v/norm(v);
h=options.step_size;
if nargin<7||isempty(do_trace), do_trace=false; end
trace=[];
if do_trace
  % PER-STEP RECORD FOR THE CONNECTION-FORM TRACKER.
  %
  % The quantity that matters here and has no counterpart in the hinec trace is
  % DRIFT: the angle between the frame the streamline is carrying and the e1 the
  % data actually has at that point. With mmf.anchor = 0 nothing re-syncs them,
  % and that drift is what costs the tracker its bundle - it was inferred from
  % saved polylines all through this investigation and never observed directly.
  % Also recorded: the curvature and torsion read from the connection field, and
  % the anchor modulation actually applied, which is spatially varying.
  tr.pos   = nan(options.max_steps+1,3);   % position BEFORE the step
  tr.e1    = nan(options.max_steps+1,3);   % carried tangent
  tr.drift = nan(options.max_steps+1,1);   % deg between carried e1 and the FIELD e1
  tr.kappa = nan(options.max_steps+1,1);   % |curvature| from the connection field
  tr.tau   = nan(options.max_steps+1,1);   % torsion
  tr.amod  = nan(options.max_steps+1,1);   % anchor modulation CL/(CL+CP) applied
  tr.aeff  = nan(options.max_steps+1,1);   % mmf.anchor * amod
  tr.h     = nan(options.max_steps+1,1);
  tr.turn  = nan(options.max_steps+1,1);   % deg turned this step
  tr.tissue= cell(options.max_steps+1,1);
  tr.psel  = zeros(options.max_steps+1,1);  % which FOD peak's curvature was used
  tr.npeak = zeros(options.max_steps+1,1);  % how many peaks were usable there
  tr.palign= nan(options.max_steps+1,1);    % alignment of the winner with the carried tangent
  ns=0;
end
% carried moving frame (Chun-Peng): e2 by reference-axis projection at the seed (Eq.6),
% then evolved along the streamline by the connection structure equation (Eq.10).
e1f=v; [e2f,e3f]=mmf_reference_axis_frame(v);
term='max_steps';
for step=1:options.max_steps
  if do_trace
    ns=ns+1; tr.pos(ns,:)=x; tr.e1(ns,:)=e1f; tr.h(ns)=h;
    [ef,of]=mmf_field(nim.ME1x,nim.ME1y,nim.ME1z,x,dims,e1f);
    if of, tr.drift(ns)=acosd(min(abs(dot(ef,e1f)),1)); end
    [~,de1t,~,~,okd,ps,nv,pa]=mmf_deriv(nim,x,e1f,e2f,e3f,dims);
    if okd, tr.kappa(ns)=norm(de1t); end
    tr.psel(ns)=ps; tr.npeak(ns)=nv; tr.palign(ns)=pa;
    if isfield(nim,'MTau'), tt=nim.MTau(x(1),x(2),x(3)); if ~isnan(tt), tr.tau(ns)=tt; end, end
    md=1; if isfield(nim,'MAmod'), md=nim.MAmod(x(1),x(2),x(3)); end
    if isnan(md), md=1; end
    tr.amod(ns)=md; tr.aeff(ns)=getdef(options,'mmf_anchor',0)*md;
  end
  [x_new,e1n,e2n,e3n,h,ok]=mmf_step(nim,x,e1f,e2f,e3f,h,dims,options); v_new=e1n;   % h adapts for rkf45
  if ~ok, term='no_direction'; break; end
  c=dot(v_new,v);
  if do_trace, tr.turn(ns)=acosd(min(abs(c),1)); end
  if c<cos_thresh, term='angle'; break; end
  [vok,vwhy]=valid_point(nim,x_new,options,dims);
  if ~vok, term=vwhy; break; end
  tissue=act_tissue(x_new,options,dims);
  if do_trace, tr.tissue{ns}=tissue; end
  if strcmp(tissue,'CSF')||strcmp(tissue,'OUTSIDE'), term=lower(tissue); break; end
  e1f=e1n; e2f=e2n; e3f=e3n; v=v_new; x=x_new; n=n+1; track(n,:)=x;   % carry the moving frame
  if strcmp(tissue,'GM'), term='gm'; break; end
end
track=track(1:n,:);
if do_trace
  k=max(ns,0);
  trace=struct('pos',tr.pos(1:k,:),'e1',tr.e1(1:k,:),'drift',tr.drift(1:k), ...
               'kappa',tr.kappa(1:k),'tau',tr.tau(1:k),'amod',tr.amod(1:k), ...
               'aeff',tr.aeff(1:k),'h',tr.h(1:k),'turn',tr.turn(1:k), ...
               'tissue',{tr.tissue(1:k)},'psel',tr.psel(1:k),'npeak',tr.npeak(1:k), ...
               'palign',tr.palign(1:k),'n_steps',k,'termination',term);
end
end

% ---- MMF connection-driven step (Chun-Peng Eq.10-11) ------------------------
% dx/ds = e1 (Eq.11); d/ds[e1;e2;e3] = [[0 w12 w13];[-w12 0 w23];[-w13 -w23 0]][e1;e2;e3]
% (Eq.10) with curvature vector kappa = grad_{e1}e1 (w12=kappa.e2, w13=kappa.e3) and
% torsion w23=tau read from the precomputed, interpolated connection field. RK4 over the
% coupled (x,e1,e2,e3) system, re-orthonormalized each step.
% The PATH is de1/ds = kappa(x) (a curvature-vector-field streamline); the full frame and
% torsion are evolved (available downstream) but -- as in Frenet, dT/ds = kappa N -- do
% not feed back into dx/ds. options.mmf_anchor in [0,1] optionally blends e1 toward the field.
% Dispatch the numerical stepping scheme (integrator). The DIRECTION always comes from the
% connection form (mmf_deriv); the integrator only decides HOW and HOW FAR to advance.
function [x_new,e1n,e2n,e3n,h_out,ok]=mmf_step(nim,x,e1,e2,e3,h,dims,options)
if strcmpi(options.integrator,'rkf45')
  [x_new,e1n,e2n,e3n,h_out,ok]=mmf_rkf45_step(nim,x,e1,e2,e3,h,dims,options);
else
  [x_new,e1n,e2n,e3n,ok]=mmf_rk4_step(nim,x,e1,e2,e3,h,dims,options); h_out=h;
end
if ok, [e1n,e2n,e3n]=mmf_anchor_blend(nim,x_new,e1n,e2n,e3n,dims,options); end
end

% ---- rk4: classic 4th-order fixed step over the coupled (x,e1,e2,e3) system ----
function [x_new,e1n,e2n,e3n,ok]=mmf_rk4_step(nim,x,e1,e2,e3,h,dims,options) %#ok<INUSD>
x_new=x; e1n=e1; e2n=e2; e3n=e3; ok=false;
[dx1,a1,b1,c1,o1]=mmf_deriv(nim,x,           e1,               e2,               e3,               dims); if ~o1, return; end
[dx2,a2,b2,c2,o2]=mmf_deriv(nim,x+0.5*h*dx1, unitv(e1+0.5*h*a1),unitv(e2+0.5*h*b1),unitv(e3+0.5*h*c1),dims); if ~o2,dx2=dx1;a2=a1;b2=b1;c2=c1;end
[dx3,a3,b3,c3,o3]=mmf_deriv(nim,x+0.5*h*dx2, unitv(e1+0.5*h*a2),unitv(e2+0.5*h*b2),unitv(e3+0.5*h*c2),dims); if ~o3,dx3=dx2;a3=a2;b3=b2;c3=c2;end
[dx4,a4,b4,c4,o4]=mmf_deriv(nim,x+h*dx3,     unitv(e1+h*a3),    unitv(e2+h*b3),    unitv(e3+h*c3),    dims); if ~o4,dx4=dx3;a4=a3;b4=b3;c4=c3;end
x_new = x  + (h/6)*(dx1+2*dx2+2*dx3+dx4);
e1n   = e1 + (h/6)*(a1+2*a2+2*a3+a4);
e2n   = e2 + (h/6)*(b1+2*b2+2*b3+b4);
e3n   = e3 + (h/6)*(c1+2*c2+2*c3+c4);
[e1n,e2n,e3n]=mmf_gram_schmidt(e1n,e2n,e3n);          % keep the moving frame orthonormal
ok=true;
end

% ---- rkf45: adaptive Dormand-Prince (embedded 5/4) over the coupled system ----
% Smaller steps in high-curvature regions via embedded error control; returns the adapted h.
function [x_new,e1n,e2n,e3n,h_out,ok]=mmf_rkf45_step(nim,x,e1,e2,e3,h,dims,options)
persistent A B5 B4
if isempty(A)
  A={ [], [1/5], [3/40 9/40], [44/45 -56/15 32/9], [19372/6561 -25360/2187 64448/6561 -212/729], ...
      [9017/3168 -355/33 46732/5247 49/176 -5103/18656], [35/384 0 500/1113 125/192 -2187/6784 11/84] };
  B5=[35/384 0 500/1113 125/192 -2187/6784 11/84 0];
  B4=[5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40];
end
x_new=x; e1n=e1; e2n=e2; e3n=e3; h_out=h; ok=false;
for attempt=1:6
  h=min(max(h,options.step_min),options.step_max);
  Kx=zeros(7,3); K1=zeros(7,3); K2=zeros(7,3); K3=zeros(7,3); good=true;
  for i=1:7
    xi=x; e1i=e1; e2i=e2; e3i=e3;
    for j=1:i-1, aij=A{i}(j); xi=xi+h*aij*Kx(j,:); e1i=e1i+h*aij*K1(j,:); e2i=e2i+h*aij*K2(j,:); e3i=e3i+h*aij*K3(j,:); end
    [d0,d1,d2,d3,oi]=mmf_deriv(nim,xi,unitv(e1i),unitv(e2i),unitv(e3i),dims);
    if ~oi, if i==1, return; else, good=false; break; end, end
    Kx(i,:)=d0; K1(i,:)=d1; K2(i,:)=d2; K3(i,:)=d3;
  end
  if ~good, h=max(options.step_min, h*0.5); continue; end
  x5=x+h*(B5*Kx); x4=x+h*(B4*Kx); err=norm(x5-x4);
  if err<=options.rkf_tol || h<=options.step_min*1.01
    x_new=x5;
    e1n=unitv(e1+h*(B5*K1)); e2n=unitv(e2+h*(B5*K2)); e3n=unitv(e3+h*(B5*K3));
    [e1n,e2n,e3n]=mmf_gram_schmidt(e1n,e2n,e3n);
    h_out=min(options.step_max, h*min(2.0, 0.9*(options.rkf_tol/max(err,1e-9))^0.2));
    ok=true; return;
  else
    h=max(options.step_min, h*max(0.2, 0.9*(options.rkf_tol/max(err,1e-9))^0.2));
  end
end
end

% optional re-anchor of e1 toward the field tangent (mmf_anchor), shared by both integrators
function [e1n,e2n,e3n]=mmf_anchor_blend(nim,x_new,e1n,e2n,e3n,dims,options)
a=options.mmf_anchor;
% Scale by the local linear fraction: full anchoring where one fibre dominates,
% none where the tensor is planar and the data cannot say which fibre this is.
if a>0 && isfield(nim,'MAmod')
  md = nim.MAmod(x_new(1),x_new(2),x_new(3));
  if ~isnan(md), a = a * md; end
end
if a>0
  [ef,of]=mmf_field(nim.ME1x,nim.ME1y,nim.ME1z,x_new,dims,e1n);
  if of, e1n=unitv((1-a)*e1n + a*ef); [e1n,e2n,e3n]=mmf_gram_schmidt(e1n,e2n,e3n); end
end
end

% structure-equation derivative (Eq.10): curvature vector kappa and torsion tau from the
% interpolated connection field; components taken in the CARRIED frame keep it orthonormal.
function [dx,de1,de2,de3,ok,psel,nvalid,palign]=mmf_deriv(nim,x,e1,e2,e3,dims)
% psel/nvalid/palign are DIAGNOSTIC outputs for the multi-peak path: which peak
% was selected, how many were usable at this point, and how well the winner
% aligned with the carried tangent. Without them there is no way to tell a live
% pathway selection from one that silently always returns peak 1 - the setup
% banner prints either way.
dx=e1; de1=[0 0 0]; de2=[0 0 0]; de3=[0 0 0]; ok=false;
psel=0; nvalid=0; palign=NaN;
if any(x<1)||x(1)>dims(1)||x(2)>dims(2)||x(3)>dims(3), return; end
tau=0;
if nim.mmf_multi
  % MULTIPLE PATHWAYS: select the peak whose direction at x best aligns with the incoming
  % tangent e1, and use THAT peak's curvature -> different approach => different continuation.
  bestal=-inf; kap=[0 0 0]; found=false;
  for p=1:numel(nim.MPk)
    Ip=nim.MPk{p}; d=[Ip{1}(x(1),x(2),x(3)),Ip{2}(x(1),x(2),x(3)),Ip{3}(x(1),x(2),x(3))];
    if any(isnan(d))||norm(d)<0.5, continue; end
    nvalid=nvalid+1;
    d=d/norm(d); al=abs(dot(d,e1));
    if al>bestal, bestal=al; Kp=nim.MKp{p}; psel=p; palign=al;
      kap=[Kp{1}(x(1),x(2),x(3)),Kp{2}(x(1),x(2),x(3)),Kp{3}(x(1),x(2),x(3))]; found=true; end
  end
  if ~found || any(isnan(kap)), ok=true; return; end        % no peak here -> straight advance
else
  kx=nim.MKx(x(1),x(2),x(3)); ky=nim.MKy(x(1),x(2),x(3)); kz=nim.MKz(x(1),x(2),x(3)); tau=nim.MTau(x(1),x(2),x(3));
  if isnan(kx)||isnan(ky)||isnan(kz), ok=true; return; end   % outside field -> straight advance
  if isnan(tau), tau=0; end
  kap=[kx,ky,kz];
end
% ANATOMICAL CURVATURE BOUND. A white-matter fascicle does not bend with a
% radius of a few millimetres; curvature above that is a fit failure, not a
% measurement. Rasterising the ISMRM-2015 ground-truth geometry into voxels puts
% its |kappa| at median 0.042 /mm, p99 0.165, p99.9 0.209 -- only 0.01% of voxels
% exceed 0.25 /mm (R = 4 mm). The estimated fields exceed it constantly: the DTI
% route 4.33% of voxels, the direct-DWI route 7.88% -- up to 788x too often, with
% a p99.9 of 1.19 /mm (R = 0.84 mm). The old bound of 2.0/voxel (= 1.0 /mm, a 1 mm
% radius) never bound anything. Tracking is sequential, so one absurd step ruins a
% streamline: the tail governs the result, not the median.
KAPPA_MAX = 0.5;   % per voxel = 0.25 /mm at 2 mm iso (R = 4 mm), above GT p99.9
kap=kap-dot(kap,e1)*e1; km=norm(kap); if km>KAPPA_MAX, kap=kap*(KAPPA_MAX/km); end
w12=dot(kap,e2); w13=dot(kap,e3);          % curvature components in the carried frame
de1 =  w12*e2 + w13*e3;                     % Eq.10:
de2 = -w12*e1 + tau*e3;
de3 = -w13*e1 - tau*e2;
ok=true;
end
function u=unitv(v), n=norm(v); if n>1e-9, u=v/n; else, u=v; end, end

% interpolate a unit-vector field at x, oriented to agree with ref (sign flip)
function [u,ok]=mmf_field(Ix,Iy,Iz,x,dims,ref)
ok=false; u=ref;
if any(x<1)||x(1)>dims(1)||x(2)>dims(2)||x(3)>dims(3), return; end
a=Ix(x(1),x(2),x(3)); b=Iy(x(1),x(2),x(3)); c=Iz(x(1),x(2),x(3));
if isnan(a)||isnan(b)||isnan(c), return; end
u=[a,b,c]; nu=norm(u); if nu<1e-6, u=ref; return; end
u=u/nu; if dot(u,ref)<0, u=-u; end
ok=true;
end

function [ok,why]=valid_point(nim,x,options,dims)
% WHY the point is unusable, not merely THAT it is.
%
% This predicate collapses four unrelated outcomes into one false: leaving the
% volume, leaving the propagation mask, landing outside the interpolation
% domain, and FA falling under the termination floor. The caller labelled every
% one of them 'outside', so a traced run reported 100% of streamlines ending
% 'outside' and the FA floor looked as though it were never applied. It was -
% the label simply hid it. Same defect, and same fix, as the hinec tracker's
% 'no_direction'.
ok=false; why='outside';
if any(x<1)||any(x>dims), return; end
xv=round(x);
if all(xv>=1)&&all(xv<=dims)&&~nim.prop_mask(xv(1),xv(2),xv(3)), why='outside_mask'; return; end
fa=nim.FA_i(x(1),x(2),x(3));
if isnan(fa), why='outside_domain'; return; end
if fa<options.termination_fa, why='fa'; return; end
ok=true; why='ok';
end

function tissue=act_tissue(pos,options,dims)
tissue='UNKNOWN';
if isempty(options.wm_mask)||isempty(options.gm_mask)||isempty(options.csf_mask), return; end
vp=round(pos); if any(vp<1)||vp(1)>dims(1)||vp(2)>dims(2)||vp(3)>dims(3), tissue='OUTSIDE'; return; end
li=sub2ind(dims,vp(1),vp(2),vp(3));
if options.csf_mask(li)>0.5, tissue='CSF'; elseif options.gm_mask(li)>0.5, tissue='GM'; elseif options.wm_mask(li)>0.5, tissue='WM'; end
end

function v=getdef(s,f,d), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end
