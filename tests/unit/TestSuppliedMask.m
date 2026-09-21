classdef TestSuppliedMask < matlab.unittest.TestCase
    properties
        Folder
    end
    methods(TestMethodSetup)
        function setup(tc)
            addpath(genpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'src')));
            tc.Folder=tempname;mkdir(tc.Folder);
            tc.addTeardown(@()rmdir(tc.Folder,'s'));
            niftiwrite(ones(3,4,5,2,'single'),fullfile(tc.Folder,'dwi.nii'));
        end
    end
    methods(Test)
        function acceptsCompressedAndUncompressed(tc)
            data=ones(3,4,5,'uint8');data(1,:,:)=0;
            for compressed=[false true]
                source=fullfile(tc.Folder,'mask.nii');
                niftiwrite(data,source,'Compressed',compressed);
                if compressed,source=[source '.gz'];end
                out=fullfile(tc.Folder,'staged.nii.gz');
                preproc_prepare_supplied_mask(source,fullfile(tc.Folder,'dwi.nii'),out);
                tc.verifyEqual(niftiread(out),data);
                fid=fopen(out,'r');magic=fread(fid,2,'uint8');fclose(fid);
                tc.verifyEqual(magic,[31;139]);
            end
        end
        function rejectsDifferentDimensions(tc)
            source=fullfile(tc.Folder,'mask.nii');niftiwrite(ones(4,4,5,'uint8'),source);
            tc.verifyError(@()preproc_prepare_supplied_mask(source,fullfile(tc.Folder,'dwi.nii'),fullfile(tc.Folder,'out.nii.gz')),'preproc:maskGrid');
        end
        function rejectsShiftedGrid(tc)
            source=fullfile(tc.Folder,'mask.nii');data=ones(3,4,5,'uint8');niftiwrite(data,source);
            info=niftiinfo(source);T=info.Transform.T;T(4,1)=T(4,1)+3;info.Transform=affine3d(T);
            info.TransformName='Sform';
            niftiwrite(data,source,info);
            tc.verifyError(@()preproc_prepare_supplied_mask(source,fullfile(tc.Folder,'dwi.nii'),fullfile(tc.Folder,'out.nii.gz')),'preproc:maskGrid');
        end
    end
end
