classdef TestParcellationFromMasks < matlab.unittest.TestCase
% The mask-derived parcellation and the affine reader it depends on.
%
% Every case here is a bug that would be invisible at runtime. A mask placed at
% the origin because its geometry was in the qform rather than the sform still
% produces a plausible-looking mask; a parcellation that drops overlaps still
% produces a plausible-looking label volume. Nothing throws. The numbers are just
% quietly wrong, which is why these are asserted rather than eyeballed.

    properties
        Root
        Tmp
        Ref
    end

    methods (TestMethodSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            tc.Root = fileparts(fileparts(here));
            addpath(genpath(tc.Root));
            tc.Tmp = tempname; mkdir(tc.Tmp);
            tc.Ref = tempname; mkdir(tc.Ref);
        end
    end
    methods (TestMethodTeardown)
        function teardown(tc)
            if isfolder(tc.Tmp), rmdir(tc.Tmp, 's'); end
            if isfolder(tc.Ref), rmdir(tc.Ref, 's'); end
        end
    end

    methods
        function f = writeMask(tc, name, vol, affine, use_qform, dest)
            if nargin < 6 || isempty(dest), dest = tc.Tmp; end
            % Write a NIfTI whose geometry is carried by the sform or, when
            % use_qform is true, ONLY by the qform - the shape a third of the
            % ISMRM ground-truth masks actually ship in.
            f = fullfile(dest, [name '.nii']);
            info = struct();
            niftiwrite(uint8(vol), f);
            % Patch the header in place: niftiwrite does not expose qform/sform
            % codes directly, and the point of the test is the raw bytes.
            fid = fopen(f, 'r+', 'l');
            if use_qform
                fwrite_at(fid, 252, int16([1 0]), 'int16');            % qform=1, sform=0
                fwrite_at(fid,  76, single(1), 'single');              % qfac
                fwrite_at(fid,  80, single(abs(diag(affine(1:3,1:3))')), 'single');
                q = mat2quat(affine(1:3,1:3));
                fwrite_at(fid, 256, single(q(2:4)), 'single');         % b, c, d
                fwrite_at(fid, 268, single(affine(1:3,4)'), 'single'); % qoffset
                fwrite_at(fid, 280, single(zeros(1,12)), 'single');    % blank srow
            else
                fwrite_at(fid, 252, int16([0 1]), 'int16');            % qform=0, sform=1
                srow = affine(1:3,:)';
                fwrite_at(fid, 280, single(srow(:)'), 'single');
            end
            fclose(fid);
        end
    end

    methods (Test)

        function sformAndQformAgree(tc)
            % The same geometry written two ways must read back identically. If
            % the qform branch is wrong, masks silently land somewhere else.
            A = [-1 0 0 0.5; 0 -1 0 0.5; 0 0 1 -0.5; 0 0 0 1];
            v = zeros(6,6,6); v(2:4,2:4,2:4) = 1;
            fs = tc.writeMask('viaS', v, A, false, tc.Ref);
            fq = tc.writeMask('viaQ', v, A, true,  tc.Ref);
            As = nim_nifti_affine(fs);
            Aq = nim_nifti_affine(fq);
            tc.verifyEqual(As, A, 'AbsTol', 1e-5, 'sform was not read back correctly.');
            tc.verifyEqual(Aq, A, 'AbsTol', 1e-5, ...
                'qform fallback disagrees with the geometry it was written from.');
        end

        function realIsmrmMasksAgreeAcrossBothConventions(tc)
            % The case that motivated this: 29 of the 75 ISMRM ROI masks carry
            % only a qform. They must land in the same frame as the ones that
            % carry an sform, or every overlap measured against them is wrong.
            R = fullfile(tc.Root,'data','ismrm2015','scoring_data_Renauld2023','ROI');
            tc.assumeTrue(isfolder(R), 'ISMRM scoring data not present.');
            A_s = nim_nifti_affine(fullfile(R,'all_masks','UF_right.nii.gz'));     % sform
            A_q = nim_nifti_affine(fullfile(R,'all_masks','CC_temporal.nii.gz'));  % qform only
            tc.verifyEqual(A_q, A_s, 'AbsTol', 1e-6, ...
                'qform-only ISMRM masks resolve to a different frame than sform ones.');
        end

        function overlapsSurviveInMasksButNotInLabels(tc)
            % Two regions sharing voxels. The label volume can keep only one
            % owner per voxel - that is what a label volume IS - but .masks must
            % return each region whole, and the overlap must be reported rather
            % than silently dropped.
            A = eye(4);
            big   = false(8,8,8); big(2:7, 2:7, 2:7) = true;   % 216 voxels
            small = false(8,8,8); small(3:5, 3:5, 3:5) = true; %  27, entirely inside big
            tc.writeMask('big',   big,   A, false);
            tc.writeMask('small', small, A, false);
            ref = tc.writeMask('ref', zeros(8,8,8), A, false, tc.Ref);

            P = nim_parcellation_from_masks({tc.Tmp}, ref, [8 8 8]);
            mb = P.masks('big'); ms = P.masks('small');
            tc.verifyEqual(sum(mb(:)), 216, 'big region lost voxels');
            tc.verifyEqual(sum(ms(:)),  27, 'small region lost voxels');
            tc.verifyEqual(P.overlap.voxels_multi, 27, ...
                'the overlap between the two regions was not reported');
            % smallest wins, so the small region is not swallowed by the big one
            small_id = find(strcmp(values(P.map), 'small'));
            tc.verifyEqual(sum(P.labels(:) == small_id), 27, ...
                'the smaller region was overwritten by the larger one in the label volume');
        end

        function attachBundleRoisKeepsThePreviousParcellation(tc)
            % The bundle parcellation must be reproducible from the repo. It was
            % originally built by an ad-hoc script and baked into a gitignored
            % .mat, so nothing in version control could recreate it - the utility
            % existed but nothing called it.
            A = eye(4);
            big   = false(8,8,8); big(2:7,2:7,2:7) = true;
            small = false(8,8,8); small(3:5,3:5,3:5) = true;
            am = fullfile(tc.Tmp, 'all_masks'); mkdir(am);
            tc.writeMask('big',   big,   A, false, am);
            tc.writeMask('small', small, A, false, am);
            ref = tc.writeMask('ref', zeros(8,8,8), A, false, tc.Ref);

            nim = struct('FA', zeros(8,8,8), 'parcellation_mask', ones(8,8,8), ...
                         'atlas_type', 'jhu', ...
                         'atlas_labels', struct('map', containers.Map(1,'old')));
            evalc('out = nim_attach_bundle_rois(nim, tc.Tmp, ref);');

            tc.verifyEqual(double(out.atlas_labels.map.Count), 2, ...
                'both bundle masks should become labels');
            tc.verifyTrue(isfield(out, 'parcellation_mask_jhu'), ...
                'the outgoing atlas parcellation must be preserved, not discarded');
            tc.verifyTrue(isa(out.roi_masks, 'containers.Map'), ...
                'the true masks must be attached for ROI selection');
            m = evalc_mask(@() nim_roi_mask(out, {'big'}, 0));
            tc.verifyEqual(sum(m(:)), 216, ...
                'a bundle must resolve to all of itself, not the label-volume remnant');
        end

        function roiMaskPrefersTheTrueMaskOverTheLabelVolume(tc)
            % The consequence of the above for ROI selection: asking for the
            % region that lost the label-volume tie must still return all of it.
            A = eye(4);
            big   = false(8,8,8); big(2:7,2:7,2:7) = true;
            small = false(8,8,8); small(3:5,3:5,3:5) = true;
            tc.writeMask('big', big, A, false);
            tc.writeMask('small', small, A, false);
            ref = tc.writeMask('ref', zeros(8,8,8), A, false, tc.Ref);
            P = nim_parcellation_from_masks({tc.Tmp}, ref, [8 8 8]);

            nim = struct('FA', zeros(8,8,8), 'parcellation_mask', P.labels, ...
                         'atlas_labels', struct('map', P.map), 'roi_masks', P.masks);
            m = evalc_mask(@() nim_roi_mask(nim, {'big'}, 0));
            tc.verifyEqual(sum(m(:)), 216, ...
                ['nim_roi_mask returned the label-volume remnant of "big" (its ' ...
                 'non-overlapping part) instead of the region as defined.']);
        end
    end
end

% =========================================================================
function fwrite_at(fid, offset, data, type)
    fseek(fid, offset, 'bof');
    fwrite(fid, data, type);
end

function q = mat2quat(R)
% Rotation matrix -> unit quaternion [a b c d]. The scale is divided out first;
% NIfTI stores orientation and voxel size separately.
    R = R ./ vecnorm(R);
    a = sqrt(max(0, 1 + R(1,1) + R(2,2) + R(3,3))) / 2;
    if a > 1e-6
        b = (R(3,2)-R(2,3)) / (4*a);
        c = (R(1,3)-R(3,1)) / (4*a);
        d = (R(2,1)-R(1,2)) / (4*a);
    else
        b = sqrt(max(0, 1 + R(1,1) - R(2,2) - R(3,3))) / 2;
        c = sqrt(max(0, 1 - R(1,1) + R(2,2) - R(3,3))) / 2;
        d = sqrt(max(0, 1 - R(1,1) - R(2,2) + R(3,3))) / 2;
    end
    q = [a b c d];
end

function m = evalc_mask(fn)
    evalc('m = fn();');
end
