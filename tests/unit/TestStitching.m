classdef TestStitching < matlab.unittest.TestCase
 methods(TestClassSetup)
  function setup(tc)
   addpath(genpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'src')));
  end
 end
 methods(Test)
  function forwardSearchFindsConnectionBehindIncompatibleClutter(tc)
   o=opts();o.min_length=0;o.stitching.neighbors=1;
   f={line([0 0 0],[2 0 0]);line([2.8 0 0],[4.8 0 0]); ...
      line([2.05 0 0],[2.05 2 0]);line([2.75 0 0],[2.75 2 0])};
   [~,m]=nim_stitch_fragments(f,o,@(~,~)true);
   tc.verifyFalse(any(all(m.joins(:,1:2)==[2 3],2)));
   o.stitching.search='forward';[~,m]=nim_stitch_fragments(f,o,@(~,~)true);
   tc.verifyTrue(any(all(m.joins(:,1:2)==[2 3],2)));
  end
  function forwardSearchKeepsMaskAndGeometryGates(tc)
   o=opts();o.min_length=0;o.stitching.search='forward';
   f={line([2 5 5],[4 5 5]);line([4.5 5 5],[6.5 5 5])};
   [~,m]=nim_stitch_fragments(f,o,@(~,~)false);tc.verifyEmpty(m.joins);
   [~,m]=nim_stitch_fragments(f,o,@(~,~)true,[0 0;0 0;1 0;1 0]);tc.verifyEmpty(m.joins);
   [~,m]=nim_stitch_fragments(f,o,@(~,~)true);tc.verifyEqual(size(m.joins,1),1);
  end
  function joinsCollinearFragmentsAndPrunes(tc)
   o=opts();o.min_length=6;o.stitching.join_radius=.6;
   f={line([2 5 5],[4 5 5]);line([4.5 5 5],[6.5 5 5]);line([7 5 5],[9 5 5]);line([2 10 5],[3 10 5])};
   [t,m]=nim_stitch_fragments(f,o,@(~,~)true);
   tc.verifyNumElements(t,1);tc.verifyEqual(t{1}(1,:),[2 5 5]);tc.verifyEqual(t{1}(end,:),[9 5 5]);
   tc.verifyEqual(numel(unique(m.fragment_ids{1})),3);tc.verifyEqual(size(m.joins,1),2);tc.verifyEqual(m.short_chains,1);
  end
  function joinsCoincidentEnds(tc)
   o=opts();o.min_length=0;[t,m]=nim_stitch_fragments({line([2 5 5],[4 5 5]);line([4 5 5],[6 5 5])},o,@(~,~)true);
   tc.verifyNumElements(t,1);tc.verifyEqual(size(m.joins,1),1);tc.verifyTrue(all(vecnorm(diff(t{1}),2,2)>0));
  end
  function doesNotConnectCrossingOrReverseGap(tc)
   o=opts();o.min_length=0;
   f={line([2 5 5],[4 5 5]);line([4.1 5 5],[4.1 7 5]);line([2.1 5 5],[3.9 5 5])};
   [t,m]=nim_stitch_fragments(f,o,@(~,~)true);tc.verifyNumElements(t,3);tc.verifyEmpty(m.joins);
  end
  function refusesUnsupportedBridge(tc)
   o=opts();o.min_length=0;
   [t,m]=nim_stitch_fragments({line([2 5 5],[4 5 5]);line([4.5 5 5],[6.5 5 5])},o,@(~,~)false);
   tc.verifyNumElements(t,2);tc.verifyEmpty(m.joins);tc.verifyGreaterThan(m.bridge_rejections,0);
  end
  function geometryAndLengthConstrainGraph(tc)
   o=opts();o.min_length=0;f={line([2 5 5],[4 5 5]);line([4.5 5 5],[6.5 5 5])};
   [~,m]=nim_stitch_fragments(f,o,@(~,~)true,[0 0;0 0;1 0;1 0]);tc.verifyEmpty(m.joins);
   o.max_arc=4;[~,m]=nim_stitch_fragments(f,o,@(~,~)true);tc.verifyEmpty(m.joins);
  end
  function denseDtiMakesLongStraightChains(tc)
   [nim,o]=phantom();[t,m]=nim_tractography_stitching(nim,o);nim_check_tracker_output(t,m);
   tc.verifyNotEmpty(t);tc.verifyGreaterThan(size(m.graph.joins,1),0);
   for i=1:numel(t),tc.verifyLessThan(max(abs(t{i}(:,2)-5)),1e-8);end
   o.stitching.join=false;[t,m]=nim_tractography_stitching(nim,o);tc.verifyEmpty(t);tc.verifyEqual(size(m.graph.joins,1),0);
  end
  function csdRetainsBothCrossingDirections(tc)
   [nim,o]=phantom();o.field='csd';o.min_length=0;o.stitching.join=false;
   o.seed_mask(:)=false;o.seed_mask(8,5,5)=true;
   nim.peaks=zeros(18,10,10,2,3);nim.peaks(:,:,:,1,1)=1;nim.peaks(:,:,:,2,2)=1;nim.npeaks=2*ones(18,10,10);
   [t,~]=nim_tractography_stitching(nim,o);tc.verifyNumElements(t,2);
   directions=cellfun(@(p)(p(end,:)-p(1,:))/norm(p(end,:)-p(1,:)),t,'UniformOutput',false);
   tc.verifyLessThan(abs(dot(directions{1},directions{2})),1e-8);
  end
  function zeroConnectionMmfPreservesStraightLines(tc)
   [nim,o]=phantom();o.stitching.geometry='mmf';nim.mmf_kappa=zeros(18,10,10,3);nim.mmf_tau=zeros(18,10,10);
   [t,m]=nim_tractography_stitching(nim,o);tc.verifyNotEmpty(t);tc.verifyTrue(m.torsion_available);
   for i=1:numel(t),tc.verifyLessThan(max(abs(t{i}(:,2)-5)),1e-8);end
  end
  function masksBlockGapBridges(tc)
   [nim,o]=phantom();nim.mask(9,:,:)=false;o.min_length=0;o.stitching.join_radius=3;
   [t,~]=nim_tractography_stitching(nim,o);
   for i=1:numel(t),tc.verifyFalse(min(t{i}(:,1))<8.5&&max(t{i}(:,1))>9.5);end
  end
  function curvedDtiAndMmfFollowLocalCircle(tc)
   [x,y,~]=ndgrid(1:18,1:18,1:5);dx=x-9;dy=y-9;r=max(hypot(dx,dy),1);
   nim.FA=ones(18,18,5);nim.mask=r>2;nim.evec=zeros(18,18,5,3,3);
   nim.evec(:,:,:,1,1)=-dy./r;nim.evec(:,:,:,2,1)=dx./r;
   nim.mmf_kappa=cat(4,-dx./r.^2,-dy./r.^2,zeros(size(r)));nim.mmf_tau=zeros(size(r));
   o=opts();o.seed_mask=false(size(r));o.seed_mask(14,9,3)=true;o.min_length=0;o.stitching.join=false;
   for mode={'direction','mmf'}
    o.stitching.geometry=mode{1};[t,~]=nim_tractography_stitching(nim,o);tc.verifyNumElements(t,1);
    rr=hypot(t{1}(:,1)-9,t{1}(:,2)-9);tc.verifyLessThan(max(abs(rr-5)),.03);
   end
  end
  function csdMmfUsesPerPeakCurvature(tc)
   [nim,o]=phantom();o.field='csd';o.stitching.geometry='mmf';o.min_length=0;o.stitching.join=false;
   o.seed_mask(:)=false;o.seed_mask(8,5,5)=true;
   nim.peaks=zeros(18,10,10,2,3);nim.peaks(:,:,:,1,1)=1;nim.peaks(:,:,:,2,2)=1;nim.npeaks=2*ones(18,10,10);
   nim.mmf_kappa=zeros(18,10,10,3);nim.mmf_kappa_p=zeros(18,10,10,2,3);
   [t,m]=nim_tractography_stitching(nim,o);tc.verifyNumElements(t,2);tc.verifyFalse(m.torsion_available);
  end
  function graphCannotCloseCycle(tc)
   o=opts();o.min_length=0;f=cell(3,1);
   for i=1:3,a=linspace((i-1)*2*pi/3,i*2*pi/3,21)';f{i}=[10+5*cos(a),10+5*sin(a),ones(size(a))*5];end
   [t,m]=nim_stitch_fragments(f,o,@(~,~)true);tc.verifyNumElements(t,1);tc.verifyEqual(size(m.joins,1),2);
   tc.verifyEqual(numel(unique(m.fragment_ids{1})),3);
  end
  function rejectsOversizedSeedPlan(tc)
   [nim,o]=phantom();o.stitching.max_fragments=1;
   tc.verifyError(@()nim_tractography_stitching(nim,o),'stitching:budget');
  end
 end
end
function o=opts()
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));c=load_config_yaml(fullfile(root,'config','stitching_dti.yml'));o=nim_config_to_options(c);o.seed_density=1;o.max_arc=100;o.min_length=6;o.trace=true;
end
function [nim,o]=phantom()
o=opts();nim.FA=ones(18,10,10);nim.mask=true(18,10,10);nim.evec=zeros(18,10,10,3,3);nim.evec(:,:,:,1,1)=1;
o.seed_mask=false(18,10,10);o.seed_mask(3:15,5,5)=true;
end
function p=line(a,b)
u=linspace(0,1,9)';p=(1-u)*a+u*b;
end
