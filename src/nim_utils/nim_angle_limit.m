function [allowed_deg, cos_allowed] = nim_angle_limit(angle_rate, arc)
% nim_angle_limit: Turn budget for one integration step of a streamline.
%
%   [allowed_deg, cos_allowed] = nim_angle_limit(angle_rate, arc)
%
%   angle_rate - degrees of turning permitted per VOXEL OF ARC. Expressing the
%                limit as a RATE is what makes it step-invariant: halving the
%                step halves the budget, so the same physical curvature bound is
%                enforced at every step size. It is a minimum radius of
%                curvature in disguise: R = (180/pi)/angle_rate voxels.
%
%   arc        - arc length this step advances, in voxels.
%
% PASS THE NOMINAL STEP SIZE, NOT THE REALIZED CHORD.
%   For dx/ds = v(x) with |v| = 1 a step advances `arc` of arc BY CONSTRUCTION.
%   The chord between the two endpoints is shorter, for two very different
%   reasons: genuine curvature shortens it by O(h^3) (negligible), while
%   disagreement among the RK stage vectors shortens it a great deal. The second
%   is a discretisation artefact, not a statement that less arc was covered.
%   Scaling the budget by the chord therefore tightens the criterion exactly
%   where the direction field is most ambiguous, and makes it METHOD-DEPENDENT:
%   Euler and RK2 advance by h*(unit vector) so their chord is identically h,
%   while RK4 averages four stage vectors and its chord runs as low as 0.25*h on
%   this data. Same config, same step size, a 4x tighter angle limit for RK4
%   alone - which is fatal in an experiment whose purpose is to compare methods.
%
% THE 90 DEGREE CLAMP IS A CEILING, NOT A SAFETY RAIL.
%   v1 is a LINE field: defined up to sign. Every tangent is sign-aligned to its
%   predecessor before the turn is measured, which confines the measured turn to
%   [0, 90] degrees. A budget above 90 degrees can never be exceeded, so the
%   criterion is INERT there. With angle_rate = 225 deg/voxel that happens for
%   every step above 0.4 voxels - the criterion silently switches itself off at
%   the coarse end of a step-size ladder and switches on again as the ladder
%   refines. Keep angle_rate * step <= 90 across a ladder if the criterion is
%   meant to be active throughout.

    if nargin < 2 || isempty(arc), arc = 1; end

    % angle_rate <= 0 DISABLES the criterion outright. This exists so that a
    % control run - "what does this tracker do with no angle termination at all?"
    % - is expressible exactly, rather than approximated by setting the rate so
    % high it happens to stop binding. Those are not the same experiment, and the
    % 90 degree ceiling above is precisely why the approximation is treacherous:
    % a "large" rate is not a loose criterion, it is a silently absent one.
    if angle_rate <= 0
        allowed_deg = Inf;
        cos_allowed = -1;   % dot(u,v) >= -1 always, so the test can never fire
        return;
    end

    allowed_deg = min(90, angle_rate * arc);
    cos_allowed = cos(deg2rad(allowed_deg));
end
