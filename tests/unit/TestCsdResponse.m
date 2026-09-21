classdef TestCsdResponse < matlab.unittest.TestCase
    methods(TestClassSetup)
        function setup(tc)
            addpath(genpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'src')));
        end
    end
    methods(Test)
        function lowFaDataHasNonzeroResponse(tc)
            nim=lowFaPhantom();
            nim=nim_csd(nim,struct('lmax',2,'n_iter',2));
            tc.verifyTrue(all(isfinite(nim.response)));
            tc.verifyGreaterThan(norm(nim.response),0);
            tc.verifyGreaterThan(nnz(nim.npeaks),0);
        end
        function emptyMaskRejectsResponseEstimation(tc)
            nim=lowFaPhantom();nim.mask(:)=false;
            tc.verifyError(@()nim_csd(nim,struct('lmax',2,'n_iter',2)),'nim_csd:noResponseVoxels');
        end
    end
end
function nim=lowFaPhantom()
n=32;i=(0:n-1)';z=1-2*(i+.5)/n;phi=i*pi*(3-sqrt(5));
g=[sqrt(1-z.^2).*cos(phi),sqrt(1-z.^2).*sin(phi),z];
nim.bval=[0;1500*ones(n,1)];nim.bvec=[0 0 0;g];
signal=[1;exp(-1500*(.0009+.0002*g(:,3).^2))];
nim.img=repmat(reshape(signal,1,1,1,[]),3,3,3,1);
nim.FA=.1*ones(3,3,3);nim.mask=true(3,3,3);
nim.eval=repmat(reshape([.0011 .0009 .0009],1,1,1,3),3,3,3,1);
nim.evec=zeros(3,3,3,3,3);nim.evec(:,:,:,3,1)=1;
end
