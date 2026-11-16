classdef TestDtRefitZeroed < matlab.unittest.TestCase
% TestDtRefitZeroed: nim_dt_refit_zeroed, the least-squares tensor of the voxels of
% the mask that the tensor fit left at zero (preprocessing.refit_zeroed_tensors).
%
% A 3 x 3 x 2 nim with 12 gradient directions at b = 1000 s/mm^2 and noise-free
% signals S = S0 exp(-b g' D g), so that a refitted tensor is the generating one to
% round-off. The voxels, with their mask and tensor as nim_dt_spd, nim_eig and nim_fa
% leave them (a hole: DT, eval and FA zero, evec the identity):
%   (1,1,1) a hole whose signals come from a positive definite tensor: refitted
%   (2,1,1) a hole whose signals come from a tensor with a negative eigenvalue (as
%           the voxels nim_dt_spd sends to BFGS): refitted, that eigenvalue floored
%   (3,1,1) a hole with only five positive samples: left at zero
%   (1,2,1) a hole with b0 = 0: left at zero
%   (2,2,1) a hole whose tensor exceeds the fit's bound of 0.01 mm^2/s: left at zero
%   (3,2,1) a voxel with its tensor: not touched
%   (1,3,1) outside the mask, without a tensor, with signals (the voxel an improved
%           mask may add): not touched until it is in the mask
%   the other voxels: outside the mask, without signals.

    properties
        G       % 12 x 3 unit gradient directions
        R       % rotation of the generating tensors
    end

    methods (TestClassSetup)
        function addPaths(tc) %#ok<MANU>
            projRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
            addpath(genpath(fullfile(projRoot, 'src')));
        end
    end

    methods (TestMethodSetup)
        function geometry(tc)
            g = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 -1 0; 1 0 -1; 0 1 -1; 1 1 1; 1 -1 1; -1 1 1];
            tc.G = g ./ sqrt(sum(g .^ 2, 2));
            a = deg2rad(30); c = deg2rad(20);
            Rz = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
            Rx = [1 0 0; 0 cos(c) -sin(c); 0 sin(c) cos(c)];
            tc.R = Rz * Rx;
        end
    end

    methods
        function D = tensor(tc, l)
            D = tc.R * diag(l) * tc.R.';
        end
        function s = signal(tc, D, s0)
            s = s0 * exp(-1000 * sum((tc.G * D) .* tc.G, 2));
        end
        function nim = phantom(tc)
            dims = [3 3 2];
            nim = struct();
            nim.bval = [0; 1000 * ones(12, 1)];
            nim.bvec = [0 0 0; tc.G];
            nim.thrsh_b0 = 5;
            nim.img_b0 = 1000 * ones(dims);
            nim.img_bi = zeros([dims 12]);
            nim.mask = zeros(dims, 'int16');
            nim.DT = zeros([dims 6]);
            nim.evec = zeros([dims 3 3]);
            nim.eval = zeros([dims 3]);
            nim.FA = zeros(dims);
            good = tc.tensor([1.7e-3 0.3e-3 0.3e-3]);
            bad  = tc.tensor([1.7e-3 0.3e-3 -0.1e-3]);
            huge = tc.tensor([2e-2 1e-3 1e-3]);
            five = tc.signal(good, 1000); five(6:end) = 0;
            sig = {1, 1, tc.signal(good, 1000); 2, 1, tc.signal(bad, 1000); 3, 1, five; ...
                   1, 2, tc.signal(good, 1000); 2, 2, tc.signal(huge, 1000); 3, 2, tc.signal(good, 1000); ...
                   1, 3, tc.signal(good, 1000)};                               % (1,3,1): outside the mask
            for k = 1:size(sig, 1)
                nim.img_bi(sig{k, 1}, sig{k, 2}, 1, :) = reshape(sig{k, 3}, 1, 1, 1, 12);
            end
            nim.img_b0(1, 2, 1) = 0;
            nim.mask(1:3, 1:2, 1) = 1;
            for x = 1:3                                                      % nim_eig on a zero DT: evec I
                for y = 1:2
                    nim.evec(x, y, 1, :, :) = eye(3);
                end
            end
            % (3,2,1) carries its tensor, as nim_dt_spd, nim_eig and nim_fa leave it
            nim.DT(3, 2, 1, :) = [good(1,1) good(2,2) good(3,3) good(1,2) good(2,3) good(1,3)];
            [Q, L] = eig(good); [l, o] = sort(diag(L).', 'descend');
            nim.evec(3, 2, 1, :, :) = Q(:, o); nim.eval(3, 2, 1, :) = l; nim.FA(3, 2, 1) = 0.8;
        end
    end

    methods (Test)
        function holesAreRefitted(tc)
            nim = tc.phantom();
            evalc('[out, n_refit, n_zero] = nim_dt_refit_zeroed(nim);');
            tc.verifyEqual(n_zero, 5);                                 % the five holes of the mask
            tc.verifyEqual(n_refit, 2);
            good = tc.tensor([1.7e-3 0.3e-3 0.3e-3]);
            dt = squeeze(out.DT(1, 1, 1, :)).';
            tc.verifyEqual(dt, [good(1,1) good(2,2) good(3,3) good(1,2) good(2,3) good(1,3)], 'AbsTol', 1e-14);
            tc.verifyEqual(squeeze(out.eval(1, 1, 1, :)).', [1.7e-3 0.3e-3 0.3e-3], 'AbsTol', 1e-14);
            e1 = squeeze(out.evec(1, 1, 1, :, 1));
            tc.verifyEqual(abs(e1.' * tc.R(:, 1)), 1, 'AbsTol', 1e-10);
            l = [1.7e-3 0.3e-3 0.3e-3];
            fa = sqrt(1/2) * sqrt((l(1) - l(2))^2 + (l(2) - l(3))^2 + (l(3) - l(1))^2) / sqrt(sum(l .^ 2));
            tc.verifyEqual(out.FA(1, 1, 1), fa, 'AbsTol', 1e-12);
        end
        function aNegativeEigenvalueIsFloored(tc)
            nim = tc.phantom();
            evalc('out = nim_dt_refit_zeroed(nim);');
            tc.verifyEqual(squeeze(out.eval(2, 1, 1, :)).', [1.7e-3 0.3e-3 1e-7], 'AbsTol', 1e-14);
            D = nim_reshape_d(squeeze(out.DT(2, 1, 1, :)));
            tc.verifyGreaterThan(min(eig(D)), 0);                     % positive definite
            tc.verifyEqual(abs(squeeze(out.evec(2, 1, 1, :, 1)).' * tc.R(:, 1)), 1, 'AbsTol', 1e-10);
            evalc('out = nim_dt_refit_zeroed(nim, 1e-6);');
            tc.verifyEqual(out.eval(2, 1, 1, 3), 1e-6, 'AbsTol', 1e-18);
        end
        function whatCannotBeFittedStaysAtZero(tc)
            nim = tc.phantom();
            evalc('out = nim_dt_refit_zeroed(nim);');
            for v = [3 1; 1 2; 2 2].'                                 % five samples, b0 = 0, above 0.01
                tc.verifyEqual(squeeze(out.DT(v(1), v(2), 1, :)), zeros(6, 1));
                tc.verifyEqual(squeeze(out.eval(v(1), v(2), 1, :)), zeros(3, 1));
                tc.verifyEqual(out.FA(v(1), v(2), 1), 0);
            end
        end
        function theRestIsNotTouched(tc)
            nim = tc.phantom();
            evalc('out = nim_dt_refit_zeroed(nim);');
            changed = any(reshape(out.DT ~= nim.DT, 3, 3, 2, []), 4);
            tc.verifyEqual(find(changed).', [1 2]);                   % (1,1,1) and (2,1,1) only
            tc.verifyEqual(out.DT(3, 2, 1, :), nim.DT(3, 2, 1, :));
            tc.verifyEqual(out.FA(3, 2, 1), 0.8);
            tc.verifyEqual(out.DT(1, 3, 1, :), nim.DT(1, 3, 1, :));    % outside the mask
            tc.verifyEqual(out.mask, nim.mask);
            tc.verifyEqual(out.img_bi, nim.img_bi);
            % a voxel whose eigenvalues are zero but whose DT is not (nim_eig skipped it,
            % e.g. for a mask value other than 1) is no hole and keeps its tensor
            nim.eval(3, 2, 1, :) = 0; nim.mask(3, 2, 1) = 2;
            evalc('[out, ~, n_zero] = nim_dt_refit_zeroed(nim);');
            tc.verifyEqual(n_zero, 5);
            tc.verifyEqual(out.DT(3, 2, 1, :), nim.DT(3, 2, 1, :));
        end
        function theVoxelsTheImprovedMaskAddsAreRefitted(tc)
            % main calls the refit again before the save: a voxel the improved mask adds
            % was outside the mask of the fit, so it has no tensor, and it is refitted then.
            nim = tc.phantom();
            nim.mask(1, 3, 1) = 1;
            evalc('[out, n_refit, n_zero] = nim_dt_refit_zeroed(nim);');
            tc.verifyEqual([n_refit n_zero], [3 6]);
            tc.verifyEqual(squeeze(out.eval(1, 3, 1, :)).', [1.7e-3 0.3e-3 0.3e-3], 'AbsTol', 1e-14);
        end
        function aNimWithoutSignalsOrFillableHolesIsReturnedAsItIs(tc)
            nim = rmfield(tc.phantom(), 'img_bi');
            out = '';
            txt = evalc('[out, n_refit, n_zero] = nim_dt_refit_zeroed(nim);');
            tc.verifyEqual(out, nim);
            tc.verifyEqual([n_refit n_zero], [0 5]);
            tc.verifySubstring(txt, 'img_bi');
            nim = tc.phantom(); nim.img_bi = nim.img_bi(:, :, :, 1:11);  % signals that do not match the b-values
            txt = evalc('[out, n_refit] = nim_dt_refit_zeroed(nim);');
            tc.verifyEqual(out, nim);
            tc.verifyEqual(n_refit, 0);
            tc.verifySubstring(txt, 'do not match');
            nim = struct('FA', ones(3, 4, 5));                        % the stub of TestExternalAtlas
            txt = evalc('out = nim_dt_refit_zeroed(nim);');
            tc.verifyEqual(out, nim);
            tc.verifyEmpty(txt);
            evalc('nim = nim_dt_refit_zeroed(tc.phantom());');        % refitted: what is left cannot be
            txt = evalc('[out, n_refit, n_zero] = nim_dt_refit_zeroed(nim);');
            tc.verifyEqual(out, nim);
            tc.verifyEqual([n_refit n_zero], [0 3]);
            tc.verifySubstring(txt, '0 of the 3');
        end
        function theKeyIsInTheSchema(tc)
            p = [tempname '.yml'];
            c = onCleanup(@() delete(p));
            fid = fopen(p, 'w'); fprintf(fid, 'tractography:\n  algorithm: mmf\n'); fclose(fid);
            evalc('cfg = load_config_yaml(p);');
            tc.verifyTrue(cfg.preprocessing.refit_zeroed_tensors);
            fid = fopen(p, 'w'); fprintf(fid, 'preprocessing:\n  refit_zeroed_tensors: false\ntractography:\n  algorithm: mmf\n'); fclose(fid);
            evalc('cfg = load_config_yaml(p);');
            tc.verifyFalse(cfg.preprocessing.refit_zeroed_tensors);
        end
    end
end
