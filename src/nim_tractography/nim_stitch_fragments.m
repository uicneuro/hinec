function [tracks, info] = nim_stitch_fragments(fragments, options, bridge_ok, geometry)
% NIM_STITCH_FRAGMENTS Assemble an endpoint graph into disjoint, acyclic chains.
% Each fragment/end is used at most once. Candidate edges are local KNN queries;
% proximity alone never creates a join. bridge_ok(points,tangents) must validate
% every bridge sample against propagation support and the direction field.
if nargin<4, geometry=[]; end
s=options.stitching; n=numel(fragments); tracks=cell(0,1);
info=struct('fragment_ids',{{}},'joins',zeros(0,3),'candidates',0,'bridge_rejections',0,'short_chains',0);
if n==0, return; end
P=zeros(2*n,3); V=P; L=zeros(n,1);
for i=1:n
 p=fragments{i}; P(2*i-1:2*i,:)=[p(1,:);p(end,:)];
 V(2*i-1:2*i,:)=[unit(p(1,:)-p(2,:));unit(p(end,:)-p(end-1,:))];
 L(i)=sum(vecnorm(diff(p),2,2));
end
edges=zeros(0,3);
if s.join && n>1
 % Bounding the neighborhood makes memory O(fragments * neighbors), not O(n^2).
 K=min(2*n,s.neighbors+1);
 tree=KDTreeSearcher(P);
 if strcmp(s.search,'forward')
  edges=nim_stitch_forward_candidates(P,V,geometry,s,tree);
 else
 for begin=1:4096:2*n
  ii=begin:min(begin+4095,2*n); [idx,dist]=knnsearch(tree,P(ii,:),'K',K);
  block=zeros(numel(ii)*K,3); count=0;
  for row=1:numel(ii)
   a=ii(row);fa=ceil(a/2);
   for k=1:K
    b=idx(row,k);gap=dist(row,k);fb=ceil(b/2);
    if b==a||fa==fb||gap>s.join_radius, continue; end
    al=dot(V(a,:),-V(b,:)); if al<cosd(s.join_angle),continue;end
    if gap>1e-8
     u=(P(b,:)-P(a,:))/gap;
     if dot(V(a,:),u)<cosd(s.join_angle)||dot(-V(b,:),u)<cosd(s.join_angle),continue;end
    end
    penalty=0;
    if ~isempty(geometry)
     % kappa and tau are reversal-invariant scalar Frenet quantities. Torsion
     % is intentionally ignored for near-straight fragments (undefined normal).
     ka=geometry(a,1);kb=geometry(b,1);
     if abs(ka-kb)>s.curvature_tolerance,continue;end
     penalty=abs(ka-kb)/s.curvature_tolerance;
     if min(ka,kb)>.01 && all(isfinite(geometry([a b],2)))
      dt=abs(geometry(a,2)-geometry(b,2));
      if dt>s.torsion_tolerance,continue;end
      penalty=penalty+dt/s.torsion_tolerance;
     end
    end
    count=count+1;block(count,:)=[min(a,b) max(a,b) gap/s.join_radius+(1-al)+.25*penalty];
   end
  end
  edges=[edges;block(1:count,:)]; %#ok<AGROW>
 end
 end
end
if ~isempty(edges),[~,keep]=unique(edges(:,1:2),'rows','stable');edges=edges(keep,:);end
info.candidates=size(edges,1);
if ~isempty(edges),edges=sortrows(edges,[3 1 2]);end
mate=zeros(2*n,1);parent=(1:n)';component_length=L;bridges=cell(2*n,1);
accepted=zeros(min(n-1,size(edges,1)),3);na=0;
for e=1:size(edges,1)
 a=edges(e,1);b=edges(e,2);
 if mate(a)||mate(b),continue;end
 ra=root(ceil(a/2));rb=root(ceil(b/2));if ra==rb,continue;end
 p=hermite_bridge(P(a,:),P(b,:),V(a,:),-V(b,:),options.step_size);
 tang=unit_rows(diff(p)); len=sum(vecnorm(diff(p),2,2));
 if len<1e-8
  if ~bridge_ok(p(1,:),V(a,:)),info.bridge_rejections=info.bridge_rejections+1;continue;end
  tang=zeros(0,3);
 end
 if component_length(ra)+component_length(rb)+len>options.max_arc,continue;end
 if ~isempty(tang)
  turn=[acosd(min(1,max(-1,dot(V(a,:),tang(1,:))))); ...
        acosd(min(1,max(-1,sum(tang(1:end-1,:).*tang(2:end,:),2)))); ...
        acosd(min(1,max(-1,dot(tang(end,:),-V(b,:)))))];
  if any(turn>s.join_angle)||~bridge_ok(p,[tang;tang(end,:)])
   info.bridge_rejections=info.bridge_rejections+1;continue;
  end
 end
 mate(a)=b;mate(b)=a;bridges{a}=p;bridges{b}=flipud(p);
 parent(rb)=ra;component_length(ra)=component_length(ra)+component_length(rb)+len;
 na=na+1;accepted(na,:)=edges(e,:);
end
info.joins=accepted(1:na,:);visited=false(n,1);ids=cell(0,1);
for endpoint=find(mate==0)'
 if visited(ceil(endpoint/2)),continue;end
 p=zeros(0,3);chain=[];entry=endpoint;
 while true
  f=ceil(entry/2);visited(f)=true;chain(end+1)=f; %#ok<AGROW>
  q=fragments{f};if mod(entry,2)==0,q=flipud(q);end
  if isempty(p),p=q;else,p=[p;q(2:end,:)];end %#ok<AGROW>
  other=2*f-mod(entry+1,2); % opposite endpoint: odd -> even, even -> odd
  next=mate(other);if next==0,break;end
  q=bridges{other};p=[p;q(2:end,:)];entry=next; %#ok<AGROW>
 end
 p=p([true;vecnorm(diff(p),2,2)>1e-9],:);
 if size(p,1)>=2 && sum(vecnorm(diff(p),2,2))>=options.min_length
  tracks{end+1,1}=p;ids{end+1,1}=chain; %#ok<AGROW>
 else,info.short_chains=info.short_chains+1;
 end
end
info.fragment_ids=ids;
 function r=root(i)
  r=i;while parent(r)~=r,r=parent(r);end
  while parent(i)~=i,j=parent(i);parent(i)=r;i=j;end
 end
end
function v=unit(v),v=v/max(norm(v),eps);end
function v=unit_rows(v),v=v./max(vecnorm(v,2,2),eps);end
function p=hermite_bridge(a,b,ta,tb,h)
g=norm(b-a);if g<1e-8,p=[a;b];return;end
u=linspace(0,1,max(3,ceil(2*g/h)+1))';
p=(2*u.^3-3*u.^2+1)*a+(u.^3-2*u.^2+u)*(g*ta)+(-2*u.^3+3*u.^2)*b+(u.^3-u.^2)*(g*tb);
end
