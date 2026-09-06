function A = nim_nifti_affine(f)
% nim_nifti_affine: The voxel-to-world affine of a NIfTI-1 file, as a 4x4
% mapping 0-BASED voxel indices to world RAS millimetres.
%
%   A = nim_nifti_affine('bundle.nii.gz')
%
% Prefers the sform, falls back to the qform, errors if neither is set.
%
% WHY BOTH. Roughly a third of the ISMRM 2015 ground-truth ROI masks carry
% sform_code = 0 and describe their geometry only through the qform quaternion.
% Reading srow blindly gets an all-zero matrix from those files, and a pipeline
% that does not notice places them at the origin - the masks silently land in the
% wrong part of the brain and every overlap number computed against them is
% wrong, with nothing obviously broken to point at.
%
% Read from the raw header rather than through a NIfTI toolbox: the toolboxes
% disagree about whether their affine is 0-based or 1-based, and getting that
% wrong shifts everything by one voxel - here 1-2 mm, easy to mistake for a
% registration error.

    path = f;
    if endsWith(lower(f), '.gz')
        d = tempname; mkdir(d);
        g = gunzip(f, d); path = g{1};
        c = onCleanup(@() rmdir(d, 's'));
    end

    fid = fopen(path, 'r', 'l');
    if fid < 0, error('nim_nifti_affine:open', 'Cannot open %s', f); end
    cl = onCleanup(@() fclose(fid));

    fseek(fid,  76, 'bof'); pixdim = fread(fid, 8, 'single');
    fseek(fid, 252, 'bof'); codes  = fread(fid, 2, 'int16');   % qform_code, sform_code
    fseek(fid, 256, 'bof'); quat   = fread(fid, 3, 'single');  % b, c, d
    fseek(fid, 268, 'bof'); qoff   = fread(fid, 3, 'single');
    fseek(fid, 280, 'bof'); srow   = fread(fid, 12, 'single');

    if codes(2) > 0
        A = [reshape(srow, 4, 3)'; 0 0 0 1];
        if any(any(A(1:3,1:3) ~= 0)), return; end
    end

    if codes(1) > 0
        b = quat(1); c = quat(2); d = quat(3);
        a2 = 1 - (b*b + c*c + d*d);
        if a2 < 0, a2 = 0; end          % rounding in the stored quaternion
        a = sqrt(a2);
        R = [ a*a+b*b-c*c-d*d, 2*(b*c-a*d),     2*(b*d+a*c); ...
              2*(b*c+a*d),     a*a+c*c-b*b-d*d, 2*(c*d-a*b); ...
              2*(b*d-a*c),     2*(c*d+a*b),     a*a+d*d-b*b-c*c ];
        qfac = pixdim(1); if qfac == 0, qfac = 1; end   % pixdim(1) is qfac, not a size
        S = diag([pixdim(2), pixdim(3), pixdim(4)*qfac]);
        A = [R*S, qoff(:); 0 0 0 1];
        return;
    end

    error('nim_nifti_affine:noTransform', ...
        '%s has neither an sform nor a qform; its position in the world is undefined.', f);
end
