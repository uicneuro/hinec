classdef TestDyadicInterpolation < matlab.unittest.TestCase
% v1 is a LINE field: the sign the eigensolver assigns each voxel is arbitrary and
% carries no anatomical meaning. Interpolating the signed components across two
% voxels that disagree on sign computes a difference rather than an average, so
% the interpolated vector collapses toward zero and its direction is noise.
%
% These tests use a field whose true orientation is CONSTANT, with the stored
% sign flipped on alternating slabs. Any correct direction interpolator must
% return the same orientation everywhere and produce a perfectly straight
% streamline; an interpolator that averages signed components cannot.

    properties
        Root
    end

    methods (TestMethodSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            tc.Root = fileparts(fileparts(here));
            addpath(genpath(tc.Root));
        end
    end

    methods
        function nim = stripedField(~, flip_signs)
            % True orientation is +x everywhere. When flip_signs is true the
            % STORED sign alternates every 4 voxels along x - a faithful model of
            % what an eigensolver actually writes out.
            n = 40; nz = 9;
            [X, ~, ~] = ndgrid(1:n, 1:n, 1:nz);
            sgn = ones(size(X));
            if flip_signs
                sgn(mod(floor((X-1)/4), 2) == 1) = -1;
            end
            evec = zeros(n, n, nz, 3, 3);
            evec(:,:,:,1,1) = sgn;          % v1 = sgn * [1 0 0]
            evec(:,:,:,2,2) = 1;
            evec(:,:,:,3,3) = 1;
            FA = ones(n, n, nz);
            nim = struct('FA', FA, 'evec', evec, ...
                         'eval', repmat(reshape([3 1 1],1,1,1,3), n, n, nz), ...
                         'mask', FA > 0);
        end

        function tr = trackOne(~, nim, step)
            seed = false(size(nim.FA));
            seed(20, 20, 5) = true;
            opts = struct('seed_mask', seed, 'seed_density', 1, 'seed_strategy', "uniform", ...
                'step_size', step, 'angle_thresh', 0, 'integration_order', 4, ...
                'interp_method', "trilinear", 'termination_fa', 0.5, ...
                'seed_fa_threshold', 0.5, 'min_length', 0, 'max_steps', ceil(80/step), ...
                'enable_diagnostics', false, 'act_enabled', false, ...
                'wm_mask', [], 'gm_mask', [], 'csf_mask', []);
            evalc('tk = nim_tractography_hinec(nim, opts);');
            tr = [];
            for i = 1:numel(tk)
                if size(tk{i},1) > size(tr,1), tr = tk{i}; end
            end
        end
    end

    methods (Test)

        function closedFormMatchesEig(tc)
            % The analytic principal eigenvector must agree with eig() to machine
            % precision on the matrices this actually sees: weighted sums of
            % dyadics, which is what an interpolated dyadic field is.
            rng(11); worst = 0;
            for t = 1:2000
                M = zeros(3);
                for j = 1:randi([1 4])
                    u = randn(1,3); u = u/norm(u); M = M + rand*(u'*u);
                end
                M = (M+M')/2;
                v = nim_principal_dir(M(1,1),M(2,2),M(3,3),M(1,2),M(1,3),M(2,3));
                tc.assertNotEmpty(v);
                [V,D] = eig(M); [~,i] = max(diag(D)); ref = V(:,i)';
                if dot(v,ref) < 0, ref = -ref; end
                worst = max(worst, norm(v - ref));
            end
            tc.verifyLessThan(worst, 1e-9, ...
                'Closed-form principal eigenvector disagrees with eig().');
        end

        function dyadicIsInvariantToStoredSign(tc)
            % The headline property. Flipping the stored sign of half the voxels
            % changes nothing about the anatomy, so it must change nothing about
            % the interpolated orientation.
            for k = 1:200
                p = [1 1 1] + rand(1,3).*[30 30 6];
                a = tc.dirAt(tc.stripedField(false), p);
                b = tc.dirAt(tc.stripedField(true),  p);
                tc.assertNotEmpty(a); tc.assertNotEmpty(b);
                tc.verifyGreaterThan(abs(dot(a,b)), 1 - 1e-9, ...
                    'Interpolated orientation changed when the stored sign flipped.');
            end
        end

        function streamlineStaysStraightThroughSignFlips(tc)
            % End to end. The true field is constant, so the streamline is a
            % straight line along x and must cross the whole volume. Under signed
            % interpolation the vector vanishes at every sign boundary and
            % tracking stops there instead.
            tr = tc.trackOne(tc.stripedField(true), 0.25);
            tc.assertGreaterThan(size(tr,1), 2, 'No streamline was produced.');
            span = max(tr(:,1)) - min(tr(:,1));
            tc.verifyGreaterThan(span, 30, ...
                'Streamline stopped early - sign boundaries are still blocking it.');
            tc.verifyLessThan(max(abs(tr(:,2) - tr(1,2))), 1e-6, ...
                'Streamline drifted off the y axis; the field is constant in +x.');
            tc.verifyLessThan(max(abs(tr(:,3) - tr(1,3))), 1e-6, ...
                'Streamline drifted off the z axis; the field is constant in +x.');
        end

        function signFlipsDoNotShortenTheStep(tc)
            % The mechanism behind the RK4 chord collapse: disagreeing stage
            % vectors produce a step shorter than h. On a constant field every
            % step must be exactly h regardless of the stored signs.
            tr = tc.trackOne(tc.stripedField(true), 0.25);
            L = sqrt(sum(diff(tr,1,1).^2, 2));
            tc.verifyLessThan(max(abs(L - 0.25)), 1e-9, ...
                'Step length deviates from the nominal step on a constant field.');
        end
    end

    methods
        function d = dirAt(~, nim, pos)
            % Build the dyadic interpolants exactly as the tracker does and read
            % the orientation at one point.
            gv = {1:size(nim.FA,1), 1:size(nim.FA,2), 1:size(nim.FA,3)};
            vx = squeeze(nim.evec(:,:,:,1,1));
            vy = squeeze(nim.evec(:,:,:,2,1));
            vz = squeeze(nim.evec(:,:,:,3,1));
            comps = {vx.*vx, vy.*vy, vz.*vz, vx.*vy, vx.*vz, vy.*vz};
            c = zeros(1,6);
            for i = 1:6
                gi = griddedInterpolant(gv, comps{i}, 'linear', 'none');
                c(i) = gi(pos(1), pos(2), pos(3));
            end
            d = nim_principal_dir(c(1), c(2), c(3), c(4), c(5), c(6));
        end
    end
end
