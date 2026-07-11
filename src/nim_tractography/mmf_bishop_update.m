function [e2n, e3n] = mmf_bishop_update(t0, t1, e2, eps_tol)
% mmf_bishop_update: Rotation-minimizing (Bishop) transport of a frame vector.
%
% Given the previous tangent t0 and new tangent t1 (both unit), transport e2 by
% the minimal rotation carrying t0 to t1 (Rodrigues), then re-orthogonalize
% against t1 to remove drift. This is the parallel-transport / Bishop frame update
% used by Phase 1. See MMF_Tractography_ImplementationPlan §2.3.
%
% Arguments:
%   t0, t1  - previous / new unit tangents (1x3 or 3x1).
%   e2      - current second frame vector to transport.
%   eps_tol - tolerance below which the tangents are treated as parallel
%             (no rotation). Default 1e-6.
% Returns:
%   e2n, e3n - transported, re-orthonormalized frame vectors (row vectors).

if nargin < 4 || isempty(eps_tol)
    eps_tol = 1e-6;
end
t0 = t0(:).'; t1 = t1(:).'; e2 = e2(:).';

a = cross(t0, t1);
na = norm(a);
if na < eps_tol
    e2n = e2;                                    % tangents parallel: no rotation
else
    ah = a / na;
    phi = atan2(na, dot(t0, t1));
    e2n = e2 * cos(phi) + cross(ah, e2) * sin(phi) + ah * (dot(ah, e2)) * (1 - cos(phi));
end

% Re-orthogonalize against the NEW tangent, then normalize.
e2n = e2n - dot(e2n, t1) * t1;
e2n = e2n / norm(e2n);
e3n = cross(t1, e2n);
end
