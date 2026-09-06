function fig = nim_plot_vs_groundtruth(tracks_src, gt_trk, ref_nii, opts)
% nim_plot_vs_groundtruth: Our streamlines against an ISMRM 2015 ground-truth
% bundle, in one shared voxel frame.
%
%   fig = nim_plot_vs_groundtruth(run_dir_or_tracks_mat, gt_trk, ref_nii, opts)
%
%   opts.max_gt      - ground-truth streamlines to draw   (default 700)
%   opts.max_ours    - our streamlines to draw            (default 700)
%   opts.title       - figure title
%   opts.save        - path stem; writes <stem>.png and <stem>.fig
%   opts.roi         - JHU label(s) to outline as the seed region
%   opts.nim         - nim struct, needed only when opts.roi is used
%
% Ground truth is grey, ours is orange. The two are drawn in the SAME frame:
% the ground truth arrives on a 1mm 180x216x180 grid and is mapped through world
% RAS into the DWI's 2mm 90x108x90 voxel space by nim_read_trk. Skipping that
% step draws the ground truth at half scale in a corner, which looks like a
% tracking failure and is really a units bug.

    if nargin < 4, opts = struct(); end
    max_gt   = getf(opts, 'max_gt',   700);
    max_ours = getf(opts, 'max_ours', 700);

    ours = load_tracks(tracks_src);
    gt   = nim_read_trk(gt_trk, ref_nii);
    fprintf('drawing %d ground-truth and %d of our streamlines\n', numel(gt), numel(ours));

    gt_s = subsample(gt, max_gt);
    if numel(ours) <= max_ours
        idx = 1:numel(ours);
    else
        idx = round(linspace(1, numel(ours), max_ours));
    end
    ours_s = ours(idx);
    if isfield(opts,'inside') && ~isempty(opts.inside)
        opts.inside_s = opts.inside(idx);
    end

    C_GT   = [0.62 0.63 0.66];
    C_OURS = [0.91 0.38 0.09];
    C_OK   = [0.10 0.69 0.48];    % inside the ground-truth bundle
    C_BAD  = [0.82 0.20 0.31];    % outside it

    % Shared limits across all four panels. Without this each panel autoscales to
    % its own projection and the bundles look like different sizes.
    allP = [vertcat(gt_s{:}); vertcat(ours_s{:})];
    lo = min(allP,[],1) - 3; hi = max(allP,[],1) + 3;

    fig = figure('Visible','off','Color','w','Position',[100 100 1460 1320]);
    tl = tiledlayout(fig, 2, 2, 'Padding','compact', 'TileSpacing','compact');
    views = { [90 0],   'sagittal (down +x)'; ...
              [180 0],  'coronal (down +y)'; ...
              [0 90],   'axial (down +z)'; ...
              [-37 22], 'oblique' };

    for k = 1:4
        ax = nexttile(tl); hold(ax,'on');
        draw(ax, gt_s, C_GT, 1.0, 0.16);
        if isfield(opts,'inside') && ~isempty(opts.inside)
            % Per-SEGMENT verdict: split each streamline into runs that agree, so
            % the drawing stays continuous instead of dashing. A streamline is not
            % right or wrong as a unit - it is right until the point where it
            % leaves the bundle, and that crossover is the thing worth seeing.
            [ok_seg, bad_seg] = split_by_verdict(ours_s, opts.inside_s);
            draw(ax, bad_seg, C_BAD, 1.4, 0.42);
            draw(ax, ok_seg,  C_OK,  1.4, 0.42);
        else
            draw(ax, ours_s, C_OURS, 1.3, 0.38);
        end
        if isfield(opts,'roi') && ~isempty(opts.roi) && isfield(opts,'nim')
            m = nim_roi_mask(opts.nim, opts.roi, 0);
            [i,j,l] = ind2sub(size(m), find(m));
            plot3(ax, i, j, l, '.', 'Color', [0.06 0.30 0.72], 'MarkerSize', 11);
        end
        xlim(ax,[lo(1) hi(1)]); ylim(ax,[lo(2) hi(2)]); zlim(ax,[lo(3) hi(3)]);
        daspect(ax,[1 1 1]); grid(ax,'on'); box(ax,'on');
        ax.GridAlpha = 0.12; ax.Color = [0.995 0.995 0.995]; ax.FontSize = 9;
        view(ax, views{k,1});
        title(ax, views{k,2}, 'FontWeight','normal', 'FontSize', 11);
        xlabel(ax,'x (vox)'); ylabel(ax,'y (vox)'); zlabel(ax,'z (vox)');
    end

    ttl = getf(opts, 'title', 'tractography vs ground truth');
    if isfield(opts,'inside') && ~isempty(opts.inside)
        sub = sprintf(['\\color[rgb]{0.62,0.63,0.66}ground truth (%d)' ...
                       '    \\color[rgb]{0.10,0.69,0.48}ours: inside bundle' ...
                       '    \\color[rgb]{0.82,0.20,0.31}ours: strayed' ...
                       '    \\color[rgb]{0.06,0.30,0.72}seed ROI'], numel(gt));
    else
        sub = sprintf(['\\color[rgb]{0.62,0.63,0.66}ground truth (%d)' ...
                       '    \\color[rgb]{0.91,0.38,0.09}ours (%d)' ...
                       '    \\color[rgb]{0.06,0.30,0.72}seed ROI'], numel(gt), numel(ours));
    end
    if isfield(opts,'note') && ~isempty(opts.note)
        sub = [sub sprintf('\n\\color[rgb]{0.2,0.2,0.2}\\fontsize{11}%s', opts.note)];
    end
    title(tl, ttl, 'FontSize', 14, 'FontWeight', 'bold');
    subtitle(tl, sub, 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');

    if isfield(opts,'save') && ~isempty(opts.save)
        print(fig, [opts.save '.png'], '-dpng', '-r140');
        savefig(fig, [opts.save '.fig']);
        fprintf('wrote %s.png and %s.fig\n', opts.save, opts.save);
    end
end

% =========================================================================
function draw(ax, T, c, lw, a)
% One line object for the whole bundle, streamlines separated by NaN. Drawing a
% plot3 per streamline creates thousands of graphics objects and takes minutes;
% this takes a moment and renders identically.
    seg = cell(numel(T), 1);
    for i = 1:numel(T)
        p = T{i};
        if size(p,1) < 2, continue; end
        if size(p,1) > 300, p = p(round(linspace(1, size(p,1), 300)), :); end
        seg{i} = [p; nan(1,3)];
    end
    P = vertcat(seg{:});
    if isempty(P), return; end
    plot3(ax, P(:,1), P(:,2), P(:,3), '-', 'Color', [c a], 'LineWidth', lw);
end

function S = subsample(T, n)
    if numel(T) <= n, S = T; return; end
    S = T(round(linspace(1, numel(T), n)));
end

function tracks = load_tracks(src)
    if iscell(src), tracks = src; return; end
    f = src;
    if isfolder(f)
        d = dir(fullfile(f, '**', 'tracks*.mat'));
        if isempty(d), error('nim_plot_vs_groundtruth:notFound','No tracks*.mat under %s', f); end
        [~,i] = max([d.datenum]); f = fullfile(d(i).folder, d(i).name);
    end
    S = load(f); tracks = S.tracks;
end

function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% =========================================================================
function [ok_seg, bad_seg] = split_by_verdict(T, inside)
% Break each streamline at every verdict change, keeping each run as its own
% polyline. Colouring segment-by-segment would draw isolated 2-point stubs;
% this keeps the geometry readable and puts a visible seam exactly where the
% streamline leaves the bundle.
    ok_seg = {}; bad_seg = {};
    for i = 1:numel(T)
        p = T{i}; v = inside{i};
        if size(p,1) < 2 || numel(v) ~= size(p,1)-1, continue; end
        b = [true; diff(v(:)) ~= 0];
        starts = find(b); stops = [starts(2:end)-1; numel(v)];
        for r = 1:numel(starts)
            run = p(starts(r):stops(r)+1, :);
            if v(starts(r)), ok_seg{end+1} = run; else, bad_seg{end+1} = run; end %#ok<AGROW>
        end
    end
    ok_seg = ok_seg(:); bad_seg = bad_seg(:);
end
