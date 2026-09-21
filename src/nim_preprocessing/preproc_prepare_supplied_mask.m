function preproc_prepare_supplied_mask(source, reference, destination)
% Validate a supplied mask on the DWI grid and stage it as compressed NIfTI.
source=char(source);reference=char(reference);destination=char(destination);
if ~isfile(source)
    error('preproc:missingMask','Configured preprocessing.mask_file not found: %s',source);
end
raw=niftiinfo(reference);mask=niftiinfo(source);
if numel(mask.ImageSize)~=3 || ~isequal(raw.ImageSize(1:3),mask.ImageSize)
    error('preproc:maskGrid','Supplied mask must be 3-D and match the DWI dimensions.');
end
if ~strcmp(raw.SpaceUnits,mask.SpaceUnits) || ...
        any(abs(raw.PixelDimensions(1:3)-mask.PixelDimensions)>1e-5) || ...
        any(abs(raw.Transform.T-mask.Transform.T)>1e-5,'all')
    error('preproc:maskGrid','Supplied mask must have the same spatial transform and voxel spacing as the DWI.');
end
if strcmp(source,destination),return;end
if endsWith(source,'.nii.gz','IgnoreCase',true)
    copyfile(source,destination);
else
    % Copying uncompressed bytes to a .gz filename does not compress them.
    scratch=tempname;mkdir(scratch);
    cleanup=onCleanup(@()rmdir(scratch,'s')); %#ok<NASGU>
    files=gzip(source,scratch);
    copyfile(files{1},destination);
end
end
