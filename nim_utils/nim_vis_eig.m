function fig = nim_vis_eig(nim, plane, slice, opts)
  arguments
    nim
    plane string {mustBeMember(plane, ["xy", "xz", "yz"])}
    slice
    opts.Mask {mustBeMember(opts.Mask, ["on", "off"])} = "on"
    opts.VoxColor {mustBeMember(opts.VoxColor, ["FA", "default"])} = "default"
    opts.VoxColorMap string = "sky"
    opts.LineSpec string = "-k"
    opts.LineColor {mustBeMember(opts.LineColor, ["FA", "default"])} = "default"
    opts.LineWidth double = 3
    opts.Position (4,1) double = [100 100 750 750]
  end
  enable_mask = opts.Mask == "on";
  if ~enable_mask
    disp("Brain mask disabled. Visaulizng whole brain.");
  end

  % Determine plane
  xrange = 0:0;
  yrange = 0:0;
  zrange = 0:0;
  hi = 0;
  vi = 0;
  size_2d = 0;
  axlabels = { "L <- x -> R", "P <- y -> A", "I <- z -> S" };

  if plane == "xy"
    xrange = 1:nim.xdim;
    yrange = 1:nim.ydim;
    zrange = slice:slice;
    size_2d = nim.xdim * nim.ydim;
    hi = 1;
    vi = 2;
  elseif plane == "xz"
    xrange = 1:nim.xdim;
    yrange = slice:slice;
    zrange = 1:nim.zdim;
    size_2d = nim.xdim * nim.zdim;
    hi = 1;
    vi = 3;
  elseif plane == "yz"
    xrange = slice:slice;
    yrange = 1:nim.ydim;
    zrange = 1:nim.zdim;
    size_2d = nim.ydim * nim.zdim;
    hi = 2;
    vi = 3;
  end

  size_v = nim.hdr.ImageSize(vi);
  size_h = nim.hdr.ImageSize(hi);

  e = zeros(size_2d, 6);
  FA_2d = zeros(size_2d, 1);

  i = 1;
  for x=xrange
    for y=yrange
      for z=zrange
        if enable_mask && nim.mask(x, y, z) == 0
          e(i, :) = [ x y z 0 0 0 ];
        else
          % Visualize principle eigenvector
          q_i = reshape(nim.evec(x, y, z, 1, :), [3 1]);

          % Normalize
          e_i = q_i ./ vecnorm(q_i);
          if any(isnan(e_i))
            % Prevent divide by zero
            e_i = double(zeros(3, 1));
          end

          e(i, :) = [ x y z e_i' ];

          % Color voxels or line with given values
          FA_2d(i) = nim.FA(x,y,z);
        end
        i = i + 1;
      end
    end
  end


  fig = figure("Name", "Principle Eigenvectors", "Position", opts.Position);

  % Color the voxels
  if opts.VoxColor ~= "default"
    % p = pcolor(reshape(vox_color, [nim.hdr.ImageSize(vi) nim.hdr.ImageSize(hi)]));
    p = pcolor(reshape(FA_2d, [size_v size_h]));
    set(gca, "Color", "#e2f0f9");
    p.EdgeColor = "none";
    colormap(opts.VoxColorMap);
    clim([0 1]);
    cbar = colorbar('westoutside');
    cbar.Label.String = opts.VoxColor;
    cbar.Label.FontSize = 26;
    hold on;
  end

  % Translate vectors to place them at the center of each voxel
  qX = e(:, hi) + 0.5 - e(:, hi+3) .* 0.5;
  qY = e(:, vi) + 0.5 - e(:, vi+3) .* 0.5;

  q = quiver(qX, qY, e(:, hi+3), e(:, vi+3), ...
    opts.LineSpec, 'LineWidth', opts.LineWidth, 'ShowArrowHead', 'off', 'AutoScale', 'off');
  daspect([1 1 1]);
  ax = gca;
  ax.FontSize = 22;
  xlabel(axlabels(hi), 'FontSize', 26);
  ylabel(axlabels(vi), 'FontSize', 26);
  title("Principle Eigenvectors", 'FontSize', 26);
  subtitle(plane + " plane, slice=" + slice + ", size=" + size_h + "x" + size_v, "FontSize", 24);

  % Highlight region to use in interp
  hxr = [61 80];
  hyr = [51 70];

  [HX,HY] = meshgrid(hxr,hyr);
  hold on; plot(HX,HY,'-y',HX',HY','-y', "LineWidth", 4); hold off;
  
  % Color the vectors
  if opts.LineColor ~= "default"
    % colorquiver(q, lin_color, opts.LineColor);
    colorquiver(q, reshape(FA_2d, [], 1), opts.LineColor);
  end
  hold off;
  axis equal;

  % Zoom in to brain
  if isfield(nim, "mask")
    mask_2d = reshape(nim.mask(xrange, yrange, zrange), [size_h size_v])';
    mask_h = max(mask_2d);
    mask_v = max(mask_2d, [], 2);

    hr = [ find(mask_h, 1) find(mask_h, 1, "last") ];
    vr = [ find(mask_v, 1) find(mask_v, 1, "last") ];
    hsz = hr(2) - hr(1);
    vsz = vr(2) - vr(1);
    szdiff = abs(hsz - vsz);

    if hsz > vsz
      xlim([ hr(1)-3 hr(2)+3 ]);
      ylim([ vr(1)-(szdiff/2)-3 vr(2)+(szdiff/2)+3]);
    else
      ylim([ vr(1)-3 vr(2)+3 ]);
      xlim([ hr(1)-(szdiff/2)-3 hr(2)+(szdiff/2)+3 ]);
    end
    daspect([1 1 1]);
  end

end


function colorquiver(q, value, label)
  arguments
    q
    value
    label string = ""
  end
  cvalue = reshape(value, [], 1);
  ccmap = colormap(gca, graycmap(0.3));
  [~, ~, ind] = histcounts(cvalue, size(ccmap, 1));
  cmap = uint8(ind2rgb(ind(:), ccmap) * 255);
  cmap(:, :, 4) = 255;
  cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
  set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(cmap(1:3, :, :), [], 4).');
  set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(cmap(1:2, :, :), [], 4).');
  cbar = colorbar(gca, 'eastoutside');

  if strlength(label) > 0
    cbar.Label.String = label;
  end
end

function alphaquiver(q, value, color)
  arguments
    q
    value
    color (1, 3) = [0 0 0]
  end
  cvalue = reshape(value, [], 1);
  ccmap = colormap(repmat(color, 256, 1));
  [~, ~, ind] = histcounts(cvalue, size(ccmap, 1));
  cmap = uint8(ind2rgb(ind(:), ccmap) * 255);
  cmap(:, :, 4) = ind - 1;
  cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
  set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(cmap(1:3, :, :), [], 4).');
  set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(cmap(1:2, :, :), [], 4).');
end

function cmap = graycmap(thrsh, lim)
  arguments
    thrsh double
    lim double = 1.0
  end

  thrsh256 = floor(256 * thrsh / lim);
  cmap = repmat(flip([zeros(1, thrsh256) linspace(0, 1, 256 - thrsh256)])', 1, 3);
end

