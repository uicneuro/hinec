function edges = nim_stitch_forward_candidates(P,V,geometry,s,tree)
% Search the whole forward cone before limiting compatible endpoints.
% Eight strata (four azimuth quadrants x two gap shells) prevent one dense
% cluster from taking the entire budget. Within each stratum prefer low cost;
% round-robin across strata, with cost and endpoint ID breaking ties.
n=size(P,1); blocks=cell(ceil(n/256),1); ca=cosd(s.join_angle);r=s.join_radius;
if s.join_angle<=45,shift=r/(2*ca);radius=shift;
else,shift=r*ca;radius=r*sind(s.join_angle);end
for start=1:256:n
 ii=start:min(start+255,n);center=P(ii,:)+shift*V(ii,:);
 ids=rangesearch(tree,center,radius+1e-10,'SortIndices',false);
 block=zeros(numel(ii)*s.neighbors,3);count=0;
 for row=1:numel(ii)
  a=ii(row);b=ids{row}(:);b(ceil(b/2)==ceil(a/2))=[];
  if isempty(b),continue;end
  delta=P(b,:)-P(a,:);gap=vecnorm(delta,2,2);u=delta./max(gap,eps);
  al=-V(b,:)*V(a,:)';
  good=gap<=r & al>=ca & (gap<=1e-8 | (u*V(a,:)'>=ca & sum(-V(b,:).*u,2)>=ca));
  b=b(good);gap=gap(good);u=u(good,:);al=al(good);penalty=zeros(size(gap));
  if isempty(b),continue;end
  if ~isempty(geometry)
   dk=abs(geometry(b,1)-geometry(a,1));good=dk<=s.curvature_tolerance;penalty=dk/s.curvature_tolerance;
   dt=abs(geometry(b,2)-geometry(a,2));use=min(geometry(b,1),geometry(a,1))>.01 & isfinite(dt);
   good=good & (~use | dt<=s.torsion_tolerance);penalty(use)=penalty(use)+dt(use)/s.torsion_tolerance;
   b=b(good);gap=gap(good);u=u(good,:);al=al(good);penalty=penalty(good);
  end
  if isempty(b),continue;end
  cost=gap/r+1-al+.25*penalty;
  [~,axis]=min(abs(V(a,:)));ref=zeros(1,3);ref(axis)=1;
  e1=cross(V(a,:),ref);e1=e1/norm(e1);e2=cross(V(a,:),e1);
  sector=1+(u*e1'>=0)+2*(u*e2'>=0)+4*(gap>r/2);
  ranked=sortrows([sector cost b],[1 2 3]);rank=zeros(size(b));
  for bin=1:8,loc=find(ranked(:,1)==bin);rank(loc)=(1:numel(loc))';end
  ranked=sortrows([rank ranked(:,2:3)],[1 2 3]);ranked=ranked(1:min(s.neighbors,size(ranked,1)),:);
  nb=size(ranked,1);dest=ranked(:,3);block(count+(1:nb),:)=[min(a,dest),max(a,dest),ranked(:,2)];count=count+nb;
 end
 blocks{ceil(start/256)}=block(1:count,:);
 if mod(start-1,256*100)==0,fprintf('stitching forward search: %d/%d endpoints\n',start,n);end
end
edges=vertcat(blocks{:});
end
