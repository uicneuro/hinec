function wijk = nim_connection_form(frames, mask, options)
% nim_connection_form: Compute the connection 1-form field omega from a frame field.
%
% The connection 1-form [omega] = dA * A' records how the orthonormal frame
% A = [e1 e2 e3] rotates through space. It is antisymmetric, so it has only 3
% independent scalar 1-forms (omega12, omega13, omega23). Evaluated in a frame
% direction e_k it is a scalar:
%
%     omega_ij(e_k) = sum_c  e_j^c * ( grad e_i^c . e_k )        (Book B-3.96)
%
% i.e. differentiate each component field of e_i in space, contract with e_k to
% get the directional derivative of e_i along e_k, then dot with e_j.
% See MMF_Tractography_ImplementationPlan §2.4 / §6.3.
%
% Arguments:
%   frames  - [D1 D2 D3 3 3], frames(:,:,:,c,i) = component c of e_i
%             (produced by nim_build_frames).
%   mask    - [D1 D2 D3] logical/numeric brain mask (>0.5 = inside). Optional;
%             pass [] to compute everywhere.
%   options - struct (currently unused knobs reserved; smoothing is applied to
%             the frames upstream in nim_build_frames, not here).
%
% Returns:
%   wijk    - [D1 D2 D3 3 3 3], wijk(:,:,:,i,j,k) = omega_ij(e_k).
%             Antisymmetry is enforced (wijk(:,:,:,i,i,k) = 0,
%             wijk(:,:,:,i,j,k) = -wijk(:,:,:,j,i,k)).

if nargin < 3
    options = struct(); %#ok<NASGU>
end

dims = size(frames);
dims = dims(1:3);

if nargin < 2 || isempty(mask)
    mask = true(dims);
else
    mask = mask > 0.5;
end

% --- Spatial Jacobian of the frame: MFJacobian(:,:,:,i,c,d) = d e_i^c / d x_d --
% MATLAB gradient returns [Fx, Fy, Fz] = [d/dim2, d/dim1, d/dim3]. We reorder to
% (dim1, dim2, dim3) = the HINEC (x, y, z) voxel-coordinate axes.
MFJacobian = zeros([dims, 3, 3, 3]);
for i = 1:3
    for c = 1:3
        F = frames(:, :, :, c, i);
        [d2, d1, d3] = gradient(F);   % d2 = d/dim2, d1 = d/dim1, d3 = d/dim3
        MFJacobian(:, :, :, i, c, 1) = d1;
        MFJacobian(:, :, :, i, c, 2) = d2;
        MFJacobian(:, :, :, i, c, 3) = d3;
    end
end

% --- omega_ij(e_k) = sum_c e_j^c * ( grad e_i^c . e_k ) -----------------------
wijk = zeros([dims, 3, 3, 3]);
for i = 1:3
    for k = 1:3
        % Directional derivative of e_i along e_k: a 3-vector field (index c).
        % dEi_dk_c = sum_d MFJacobian(:,:,:,i,c,d) * e_k^d
        dEi_dk = zeros([dims, 3]);
        for c = 1:3
            acc = zeros(dims);
            for d = 1:3
                acc = acc + MFJacobian(:, :, :, i, c, d) .* frames(:, :, :, d, k);
            end
            dEi_dk(:, :, :, c) = acc;
        end
        for j = 1:3
            % Contract with e_j.
            s = zeros(dims);
            for c = 1:3
                s = s + frames(:, :, :, c, j) .* dEi_dk(:, :, :, c);
            end
            wijk(:, :, :, i, j, k) = s;
        end
    end
end

% --- Enforce antisymmetry (numerical hygiene) --------------------------------
for k = 1:3
    for i = 1:3
        wijk(:, :, :, i, i, k) = 0;
        for j = (i+1):3
            a = wijk(:, :, :, i, j, k);
            b = wijk(:, :, :, j, i, k);
            wijk(:, :, :, i, j, k) = (a - b) / 2;
            wijk(:, :, :, j, i, k) = -(a - b) / 2;
        end
    end
end

% --- Zero outside the mask ----------------------------------------------------
if ~all(mask(:))
    outside = ~mask;
    for i = 1:3
        for j = 1:3
            for k = 1:3
                slc = wijk(:, :, :, i, j, k);
                slc(outside) = 0;
                wijk(:, :, :, i, j, k) = slc;
            end
        end
    end
end

end
