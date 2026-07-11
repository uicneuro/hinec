function [e1, e2, e3] = mmf_gram_schmidt(e1, e2, e3)
% mmf_gram_schmidt: Re-orthonormalize a nearly-orthonormal right-handed frame.
%
% Guards against numerical drift of the moving frame during integration:
% normalize e1, remove its component from e2 and normalize, then e3 = e1 x e2.
%
% Arguments/Returns: e1, e2, e3 as 1x3 or 3x1 vectors (returned as rows).

e1 = e1(:).'; e2 = e2(:).'; e3 = e3(:).'; %#ok<NASGU>

e1 = e1 / norm(e1);
e2 = e2 - dot(e2, e1) * e1;
e2 = e2 / norm(e2);
e3 = cross(e1, e2);
end
