classdef TestTrackerInterface < matlab.unittest.TestCase
% The tracker hand-off: what runTractography gives a tracker, what it must get
% back, and that the shipped template honours both.
%
% Three things are pinned here, because each has silently drifted before:
%   * the input description covers every field a run actually hands over - an
%     "(undocumented)" line in tracker_input.txt is a boundary defect;
%   * the output check rejects the shapes a new tracker is most likely to get
%     wrong (N x 2 points, meta indexed by seed instead of by kept track);
%   * nim_tractography_template, the worked example people will copy, tracks a
%     known field correctly and returns exactly the contract.

    properties
        Nim
        Opts
    end

    methods (TestMethodSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fileparts(fileparts(here))));
            % Rigid rotation about (C,C): circles of every radius, FA 1 in an
            % annulus. Same phantom as TestIntegratorOrder.
            C = 30.5; n = 60; nz = 9;
            [X, Y, ~] = ndgrid(1:n, 1:n, 1:nz);
            dx = X - C; dy = Y - C; r = sqrt(dx.^2 + dy.^2); r(r < 1e-9) = 1e-9;
            evec = zeros(n, n, nz, 3, 3);
            evec(:,:,:,1,1) = -dy./r; evec(:,:,:,2,1) = dx./r;
            evec(:,:,:,1,2) =  dx./r; evec(:,:,:,2,2) = dy./r; evec(:,:,:,3,3) = 1;
            FA = zeros(n, n, nz); FA(r > 4 & r < 27) = 1;
            tc.Nim = struct('FA', FA, 'evec', evec, ...
                'eval', repmat(reshape([3 1 1], 1, 1, 1, 3), n, n, nz), 'mask', FA > 0);
            seed = false(size(FA)); seed(round(C + 15), round(C), 5) = true;

            % The full option set a config produces, plus what runTractography
            % adds in steps 4 and 4b - i.e. what a tracker really receives.
            cfg = load_config_yaml(fullfile(fileparts(fileparts(here)), 'config', 'hinec_dti.yml'));
            tc.Opts = nim_config_to_options(cfg);
            tc.Opts.seed_mask = seed;
            tc.Opts.wm_mask = []; tc.Opts.gm_mask = []; tc.Opts.csf_mask = [];
            tc.Opts.step_size = 0.25; tc.Opts.max_arc = 60; tc.Opts.max_steps = 240;
            tc.Opts.min_length = 0; tc.Opts.seed_density = 1;
        end
    end

    methods (Test)
        function everyHandedOverFieldIsDocumented(tc)
            out = fullfile(tempdir, 'tracker_input_test.txt');
            txt = evalc('nim_describe_tracker_input(tc.Nim, tc.Opts, ''hinec'', out);');
            tc.verifyTrue(isfile(out), 'tracker_input.txt was not written');
            tc.verifyEmpty(strfind(txt, '(undocumented)'), ...
                sprintf('undocumented field in the hand-off:\n%s', txt));
            % The layouts a new author needs must be stated, not implied.
            tc.verifySubstring(txt, '[X Y Z 3 3]');
            tc.verifySubstring(txt, 'seed_mask');
            tc.verifySubstring(txt, 'angle_thresh');
            tc.verifySubstring(txt, 'output contract');
            delete(out);
        end

        function outputCheckRejectsWrongShapes(tc)
            good = {[1 2 3; 2 3 4; 3 4 5]};
            evalc('nim_check_tracker_output(good, struct());');   % must not throw
            tc.verifyError(@() nim_check_tracker_output({[1 2; 2 3]}, struct()), 'tracker:outputContract');
            tc.verifyError(@() nim_check_tracker_output({[1 2 3]}, struct()), 'tracker:outputContract');
            tc.verifyError(@() nim_check_tracker_output(good, struct('seed_index', [1 2])), 'tracker:outputContract');
            tc.verifyError(@() nim_check_tracker_output(good, struct('seed_points', [1 2 3; 4 5 6])), 'tracker:outputContract');
            tc.verifyError(@() nim_check_tracker_output([1 2 3], struct()), 'tracker:outputContract');
        end

        function templateFollowsTheCircleAndHonoursTheContract(tc)
            evalc('[T, meta] = nim_tractography_template(tc.Nim, tc.Opts);');
            evalc('nim_check_tracker_output(T, meta);');
            tc.assertNumElements(T, 1);
            tc.verifyEqual(meta.n_seeds, 1);
            tc.verifyEqual(meta.seed_index, 1);
            tc.verifyEqual(size(meta.seed_points), [1 3]);

            tr = T{1};
            % The seed is on the track, once, and the halves are joined through it.
            seedrow = find(all(abs(tr - meta.seed_points) < 1e-12, 2));
            tc.verifyNumElements(seedrow, 1);
            tc.verifyGreaterThan(seedrow, 1);  tc.verifyLessThan(seedrow, size(tr, 1));
            % Consecutive points are one step apart: it really is the polyline.
            steps = sqrt(sum(diff(tr, 1, 1).^2, 2));
            tc.verifyEqual(steps, 0.25 * ones(size(steps)), 'AbsTol', 1e-9);
            % Euler on a circle of radius 15 with h = 0.25 drifts outward by
            % ~h/(2R) per step; over 60 voxels of arc that is well under a voxel.
            r = sqrt(sum((tr(:, 1:2) - 30.5).^2, 2));
            tc.verifyLessThan(max(abs(r - 15)), 1.0);
            tc.verifyGreaterThan(size(tr, 1), 400);   % both halves ran to max_arc
        end
    end
end
