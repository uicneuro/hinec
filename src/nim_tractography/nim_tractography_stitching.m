function [tracks, meta] = nim_tractography_stitching(nim, options)
% NIM_TRACTOGRAPHY_STITCHING Dense short tractlets + compatible endpoint graph.
% See docs/STITCHING.md. All lengths are voxel-space arc lengths. No GT or
% bundle gates are read. DTI/CSD supply directions; stitching.geometry=mmf
% evolves a local moving frame using the supplied connection field.
s=options.stitching; dims=size(nim.FA); h=options.step_size;
assert(s.neighbors==round(s.neighbors)&&s.max_fragments==round(s.max_fragments),'stitching:integer','neighbors and max_fragments must be integers.');
assert(s.min_fragment_arc<=s.fragment_arc&&s.fragment_arc<=options.max_arc,'stitching:arc','Require min_fragment_arc <= fragment_arc <= max_arc.');
if ~ismember(lower(options.field),{'dti','csd'}),error('stitching:field','Use field dti or csd.');end
if ~strcmpi(options.integrator,'rk2')||~strcmpi(options.interp_method,'trilinear')||options.upsample~=1
 error('stitching:settings','Stitching currently requires rk2, trilinear, and upsample=1.');
end
if ~strcmpi(options.seed_strategy,'uniform'),error('stitching:seeding','Stitching requires uniform seeding.');end
if options.act_enabled,error('stitching:act','ACT is not supported by stitching; use act: false.');end
iscsd=strcmpi(options.field,'csd');ismmf=strcmpi(s.geometry,'mmf');
if iscsd,assert(isfield(nim,'peaks')&&isfield(nim,'npeaks'),'stitching:csd','CSD peaks required.');end
if ismmf,assert(isfield(nim,'mmf_kappa'),'stitching:mmf','Pipeline must build MMF geometry first.');end
if ismmf&&iscsd,assert(isfield(nim,'mmf_kappa_p'),'stitching:mmfCsd','Per-peak MMF curvature required.');end
gv={1:dims(1),1:dims(2),1:dims(3)};FA=griddedInterpolant(gv,double(nim.FA),'linear','none');
v=double(reshape(nim.evec(:,:,:,:,1),[dims 3]));pairs=[1 1;2 2;3 3;1 2;1 3;2 3];D=cell(1,6);
for j=1:6,D{j}=griddedInterpolant(gv,v(:,:,:,pairs(j,1)).*v(:,:,:,pairs(j,2)),'linear','none');end
KI={};TI=[];
if ismmf&&~iscsd
 for j=1:3,KI{j}=griddedInterpolant(gv,double(nim.mmf_kappa(:,:,:,j)),'linear','none');end
 TI=griddedInterpolant(gv,double(nim.mmf_tau),'linear','none');
end
[offsets,~]=nim_seed_offsets(options.seed_density);
[ix,iy,iz]=ind2sub(dims,find(options.seed_mask));centres=[ix iy iz];ns=size(centres,1)*size(offsets,1);
np=1;if iscsd,np=size(nim.peaks,4);end
if ns*np>s.max_fragments
 error('stitching:budget','Potential %d fragments exceeds stitching.max_fragments=%d. Reduce ROI/density or raise the explicit budget.',ns*np,s.max_fragments);
end
fragments=cell(ns*np,1);seedids=zeros(ns*np,1);seedpos=zeros(ns*np,3);geo=zeros(2*ns*np,2);nf=0;si=0;
for offset=1:size(offsets,1)
 for c=1:size(centres,1)
  si=si+1;seed=centres(c,:)+offsets(offset,:);
  if ~alive(seed),continue;end
  if iscsd
   ij=round(seed);nr=double(nim.npeaks(ij(1),ij(2),ij(3)));
   starts=reshape(nim.peaks(ij(1),ij(2),ij(3),1:nr,:),nr,3);
  else,[d,~,~]=field(seed,[]);starts=d;end
  for peak=1:size(starts,1)
   d=starts(peak,:);if norm(d)<.5||any(~isfinite(d)),continue;end;d=unit(d);
   a=half(seed,-d);b=half(seed,d);p=[flipud(a(2:end,:));b];
   if size(p,1)<3||sum(vecnorm(diff(p),2,2))<s.min_fragment_arc,continue;end
   nf=nf+1;fragments{nf}=p;seedids(nf)=si;seedpos(nf,:)=seed;
   if ismmf
    for endpoint=1:2
     if endpoint==1,x=p(1,:);ref=unit(p(2,:)-p(1,:));else,x=p(end,:);ref=unit(p(end,:)-p(end-1,:));end
     [~,k,tau]=field(x,ref);geo(2*nf-2+endpoint,:)=[norm(k) tau];
    end
   end
  end
 end
 fprintf('stitching: seed offset %d/%d; %d fragments\n',offset,size(offsets,1),nf);
end
fragments=fragments(1:nf);seedids=seedids(1:nf);seedpos=seedpos(1:nf,:);geo=geo(1:2*nf,:);
if ~ismmf,geo=[];end
[tracks,graph]=nim_stitch_fragments(fragments,options,@bridge_ok,geo);
nt=numel(tracks);first=zeros(nt,1);source_seeds=cell(nt,1);
for i=1:nt,ids=graph.fragment_ids{i};first(i)=ids(1);source_seeds{i}=seedids(ids)';end
meta=struct('n_seeds',ns,'n_fragments',nf,'seed_index',seedids(first)', ...
 'seed_points',seedpos(first,:),'fragment_seed_indices',{source_seeds},'graph',graph, ...
 'geometry',s.geometry,'field',options.field,'torsion_available',ismmf&&~iscsd);
if options.trace,meta.fragments=fragments;meta.fragment_geometry=geo;meta.fragment_seeds=seedpos;end
fprintf('stitching: %d fragments, %d candidate joins, %d accepted joins, %d retained chains\n',nf,graph.candidates,size(graph.joins,1),nt);
 function ok=alive(p)
  ok=all(isfinite(p))&&all(p>=1)&&all(p<=dims);
  if ~ok,return;end
  ij=round(p);
  if isfield(nim,'mask')&&~isempty(nim.mask),ok=nim.mask(ij(1),ij(2),ij(3))>.5;end
  if ok&&isfield(options,'propagation_mask')&&~isempty(options.propagation_mask),ok=options.propagation_mask(ij(1),ij(2),ij(3))>0;end
  if ok,ok=FA(p(1),p(2),p(3))>=options.termination_fa;end
 end
 function [d,k,tau]=field(p,ref)
  d=[];k=[0 0 0];tau=NaN;if ~alive(p),return;end
  if ~iscsd
   z=zeros(1,6);for t=1:6,z(t)=D{t}(p(1),p(2),p(3));end
   d=nim_principal_dir(z(1),z(2),z(3),z(4),z(5),z(6));d=d(:)';
   if ismmf
    for t=1:3,k(t)=KI{t}(p(1),p(2),p(3));end
    tau=TI(p(1),p(2),p(3));
   end
  else
   if isempty(ref),ij=round(p);ref=reshape(nim.peaks(ij(1),ij(2),ij(3),1,:),1,3);end
   base=min(floor(p),dims-1);u=p-base;A=zeros(3);weight=0;
   for a=0:1,for b0=0:1,for c0=0:1
    ij=base+[a b0 c0];w=prod(([a b0 c0].*u)+(1-[a b0 c0]).*(1-u));if w<eps,continue;end
    nr=double(nim.npeaks(ij(1),ij(2),ij(3)));if nr<1,continue;end
    ps=double(reshape(nim.peaks(ij(1),ij(2),ij(3),1:nr,:),nr,3));norms=vecnorm(ps,2,2);ps=ps./max(norms,eps);
    score=abs(ps*ref(:));score(norms<.5)=-inf;[best,j]=max(score);if ~isfinite(best),continue;end
    v0=ps(j,:);A=A+w*(v0'*v0);weight=weight+w;
    if ismmf,k=k+w*double(reshape(nim.mmf_kappa_p(ij(1),ij(2),ij(3),j,:),1,3));end
   end,end,end
   if weight<.5,return;end
   d=nim_principal_dir(A(1,1),A(2,2),A(3,3),A(1,2),A(1,3),A(2,3));d=d(:)';k=k/weight;
  end
  if isempty(d)||any(~isfinite(d))||any(~isfinite(k)),d=[];return;end
  if ~isempty(ref)&&dot(d,ref)<0,d=-d;end
 end
 function p=half(seed,d)
  [n,b0]=mmf_reference_axis_frame(d);p=seed;x=seed;arc=0;
  while arc<s.fragment_arc/2-1e-10
   dt=min(h,s.fragment_arc/2-arc);
   if ismmf
    [~,k,tau]=field(x,d);[a,c1,c2]=deriv(d,n,b0,k,tau);
    [dm,nm,bm]=mmf_gram_schmidt(d+.5*dt*a,n+.5*dt*c1,b0+.5*dt*c2);
    mid=x+.5*dt*d;[ff,km,tm]=field(mid,dm);if isempty(ff),break;end
    [a,c1,c2]=deriv(dm,nm,bm,km,tm);
    next=x+dt*dm;[dn,nn,bn]=mmf_gram_schmidt(d+dt*a,n+dt*c1,b0+dt*c2);
   else
    [f,~,~]=field(x,d);if isempty(f),break;end
    [fm,~,~]=field(x+.5*dt*f,f);if isempty(fm),break;end
    next=x+dt*fm;[dn,~,~]=field(next,fm);if isempty(dn),break;end
   end
   if ~alive(next),break;end
   if options.angle_thresh>0&&acosd(min(1,max(-1,dot(d,dn))))>options.angle_thresh*dt,break;end
   p(end+1,:)=next;x=next;d=dn;arc=arc+dt; %#ok<AGROW>
   if ismmf,n=nn;b0=bn;end
  end
 end
 function [a,b,c]=deriv(t,n,bn,k,tau)
  k=k-dot(k,t)*t;km=norm(k);if km>s.max_curvature,k=k*s.max_curvature/km;end
  if ~isfinite(tau),tau=0;end % CSD has no per-peak torsion: no binormal twist.
  w1=dot(k,n);w2=dot(k,bn);a=w1*n+w2*bn;b=-w1*t+tau*bn;c=-w2*t-tau*n;
 end
 function ok=bridge_ok(p,t)
  ok=true;
  for z=1:size(p,1)
   if ~alive(p(z,:)),ok=false;return;end
   [f,~,~]=field(p(z,:),t(z,:));
   if isempty(f)||abs(dot(f,t(z,:)))<cosd(s.bridge_angle),ok=false;return;end
  end
 end
end
function v=unit(v),v=v/max(norm(v),eps);end
