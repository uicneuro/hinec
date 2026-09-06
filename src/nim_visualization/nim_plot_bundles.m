function fig = nim_plot_bundles(items, opts)
% nim_plot_bundles: Render one or more streamline sets in a shared voxel frame.
%
%   fig = nim_plot_bundles(items, opts)
%
% items is a struct array with fields
%   .tracks  - cell array of Nx3 voxel-coordinate polylines
%   .name    - label for the legend
%   .color   - 1x3 RGB, or [] to colour each segment by its direction
%   .alpha   - line alpha (default 0.35)
%   .width   - line width (default 1.2)
%   .clip    - optional logical volume; keep only the longest contiguous run of
%              each streamline that lies inside it, and report how much length
%              was dropped in items(i).frac_outside
%
% opts
%   .views   - cell of {[az el], title} rows (default: sagittal/coronal/axial/oblique)
%   .title   - figure title
%   .save    - path stem; writes <stem>.png
%   .dpi     - export resolution (default 150)
%   .maxper  - streamlines drawn per item (default 900)
%
% DIRECTION COLOURING (.color empty) is the convention the field reads by
% default: |dx| -> red, |dy| -> green, |dz| -> blue, so left-right tracts are red,
% anterior-posterior green and superior-inferior blue. It is applied per SEGMENT,
% which is why a whole-brain tractogram renders as anatomy rather than a monotone
% hairball.

    if nargin < 2, opts = struct(); end
    views  = getf(opts, 'views', { [90 0],'sagittal'; [180 0],'coronal'; ...
                                   [0 90],'axial';    [-37 22],'oblique' });
    maxper = getf(opts, 'maxper', 900);
    dpi    = getf(opts, 'dpi', 150);
    nv     = size(views, 1);

    % shared limits so the panels are comparable
    allP = [];
    for i = 1:numel(items)
        T = items(i).tracks;
        if isfield(items, 'clip') && ~isempty(items(i).clip)
            [T, items(i).frac_outside] = clip_to(T, items(i).clip);  %#ok<AGROW>
        end
        T = subsample(T, maxper);
        items(i).drawn = T;                                     %#ok<AGROW>
        allP = [allP; vertcat(T{:})];                           %#ok<AGROW>
    end
    lo = min(allP,[],1) - 3; hi = max(allP,[],1) + 3;

    nc = min(nv, 2); nr = ceil(nv/nc);
    fig = figure('Visible','off','Color','w','Position',[0 0 720*nc 620*nr]);
    tl = tiledlayout(fig, nr, nc, 'Padding','compact', 'TileSpacing','compact');

    for v = 1:nv
        ax = nexttile(tl); hold(ax,'on');
        for i = 1:numel(items)
            a = getf_s(items(i), 'alpha', 0.35);
            w = getf_s(items(i), 'width', 1.2);
            if isempty(items(i).color)
                draw_dec(ax, items(i).drawn, a, w);
            else
                draw_flat(ax, items(i).drawn, items(i).color, a, w);
            end
        end
        xlim(ax,[lo(1) hi(1)]); ylim(ax,[lo(2) hi(2)]); zlim(ax,[lo(3) hi(3)]);
        daspect(ax,[1 1 1]); grid(ax,'on'); box(ax,'on');
        ax.GridAlpha = 0.10; ax.Color = [0.995 0.995 0.995]; ax.FontSize = 9;
        view(ax, views{v,1});
        title(ax, views{v,2}, 'FontWeight','normal', 'FontSize', 11);
        xlabel(ax,'x (vox)'); ylabel(ax,'y (vox)'); zlabel(ax,'z (vox)');
    end

    if isfield(opts,'title') && ~isempty(opts.title)
        title(tl, opts.title, 'FontSize', 14, 'FontWeight','bold');
    end
    if isfield(opts,'subtitle') && ~isempty(opts.subtitle)
        subtitle(tl, opts.subtitle, 'FontSize', 11.5, 'Interpreter','tex');
    end
    if isfield(opts,'save') && ~isempty(opts.save)
        exportgraphics(fig, [opts.save '.png'], 'Resolution', dpi);
        fprintf('wrote %s.png\n', opts.save);
    end
end

% =========================================================================
function [out, frac_outside] = clip_to(T, mask)
% Keep, per streamline, the longest contiguous run of points inside the mask.
%
% WHY CLIP RATHER THAN REJECT. The ISMRM scorer defines a bundle by containment:
% a streamline with any point outside the corridor is not that bundle, and the
% whole streamline goes. That is the right rule for scoring and the wrong one for
% a picture - a streamline that follows a bundle correctly and then overshoots
% past its end is deleted entirely, taking the correct portion with it, so the
% figure shows a bundle that stops short. That reads as the tracker failing to
% reach when what happened is that it reached too far. Clipping shows the reach
% that is genuinely inside the bundle and reports the length that was not, rather
% than hiding the overshoot by deleting its streamline.
%
% Rejection still belongs in nim_filter_tracks_roi, which is what scoring uses.
    dims = size(mask);
    out = {}; len_in = 0; len_tot = 0;
    for t = 1:numel(T)
        p = T{t};
        if size(p,1) < 2, continue; end
        len_tot = len_tot + sum(sqrt(sum(diff(p,1,1).^2, 2)));
        v = round(p);
        ok = all(v >= 1, 2) & v(:,1) <= dims(1) & v(:,2) <= dims(2) & v(:,3) <= dims(3);
        inside = false(size(p,1),1);
        inside(ok) = mask(sub2ind(dims, v(ok,1), v(ok,2), v(ok,3)));
        d = diff([false; inside; false]);
        st = find(d == 1); en = find(d == -1) - 1;
        if isempty(st), continue; end
        [rl, k] = max(en - st + 1);
        if rl < 2, continue; end
        q = p(st(k):en(k), :);
        out{end+1} = q;                                          %#ok<AGROW>
        len_in = len_in + sum(sqrt(sum(diff(q,1,1).^2, 2)));
    end
    out = out(:);
    frac_outside = 1 - len_in / max(len_tot, eps);
end

function draw_flat(ax, T, c, a, w)
% One line object for the whole set, streamlines separated by NaN. A plot3 per
% streamline creates thousands of graphics objects and takes minutes.
    seg = cell(numel(T),1);
    for i = 1:numel(T)
        p = decimate_pts(T{i});
        if isempty(p), continue; end
        seg{i} = [p; nan(1,3)];
    end
    P = vertcat(seg{:});
    if isempty(P), return; end
    plot3(ax, P(:,1), P(:,2), P(:,3), '-', 'Color', [c a], 'LineWidth', w);
end

function draw_dec(ax, T, a, w)
% Direction-encoded colour. Segments are bucketed by orientation and each bucket
% drawn as one line object: per-segment colouring with a surface or patch is
% exact but produces a figure too heavy to export at print resolution.
    nb = 12;
    buckets = repmat({{}}, nb, nb);
    for i = 1:numel(T)
        p = decimate_pts(T{i});
        if size(p,1) < 2, continue; end
        d = abs(diff(p,1,1));
        n = sqrt(sum(d.^2,2)); n(n<1e-9) = 1;
        d = d ./ n;
        % index by the two largest-variation axes; blue falls out as the rest
        ia = min(nb, max(1, ceil(d(:,1)*nb)));
        ib = min(nb, max(1, ceil(d(:,2)*nb)));
        for k = 1:size(d,1)
            buckets{ia(k),ib(k)}{end+1} = [p(k:k+1,:); nan(1,3)];
        end
    end
    for ia = 1:nb
        for ib = 1:nb
            if isempty(buckets{ia,ib}), continue; end
            P = vertcat(buckets{ia,ib}{:});
            r = (ia-0.5)/nb; g = (ib-0.5)/nb;
            b = sqrt(max(0, 1 - r^2 - g^2));
            m = max([r g b]); if m > 0, c = [r g b]/m; else, c = [0.5 0.5 0.5]; end
            plot3(ax, P(:,1), P(:,2), P(:,3), '-', 'Color', [c a], 'LineWidth', w);
        end
    end
end

function p = decimate_pts(p)
    if size(p,1) < 2, p = []; return; end
    if size(p,1) > 220, p = p(round(linspace(1,size(p,1),220)), :); end
end

function S = subsample(T, n)
    T = T(cellfun(@(p) size(p,1) > 1, T));
    if numel(T) <= n, S = T; return; end
    S = T(round(linspace(1, numel(T), n)));
end

function v = getf(s,f,d),   if isfield(s,f) && ~isempty(s.(f)), v=s.(f); else, v=d; end, end
function v = getf_s(s,f,d), if isfield(s,f) && ~isempty(s.(f)), v=s.(f); else, v=d; end, end
