function [e2, e3] = mmf_reference_axis_frame(e1)
% mmf_reference_axis_frame: Build e2, e3 from e1 alone (reference-axis projection).
%
% Robust construction of the orthonormal complement of a unit tangent e1, used to
% initialize/define the moving frame. Avoids the DTI e2/e3 degeneracy when
% lambda2 ~ lambda3. See MMF_Tractography_ImplementationPlan §2.2 (D-6).
%
%   r  = Cartesian unit axis least aligned with e1 (smallest |e1.axis|)
%   e2 = normalize( r - (r.e1) e1 )
%   e3 = e1 x e2
%
% Arguments:
%   e1 - 1x3 or 3x1 unit vector (fiber tangent).
% Returns:
%   e2, e3 - orthonormal vectors completing the frame (row vectors).

e1 = e1(:).';
e1 = e1 / norm(e1);

[~, ax] = min(abs(e1));
r = [0 0 0];
r(ax) = 1;

e2 = r - dot(r, e1) * e1;
e2 = e2 / norm(e2);
e3 = cross(e1, e2);
end
