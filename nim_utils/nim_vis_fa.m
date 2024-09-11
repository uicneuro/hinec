function fig = nim_vis_fa(nim, plane, slice, opts)
  arguments
    nim
    plane
    slice
    opts.Mask {mustBeMember(opts.Mask, ["on", "off"])} = "on"
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
  size2 = 0;
  axlabels = { "L <- x -> R", "P <- y -> A", "I <- z -> S" };

  if plane == "xy"
    xrange = 1:nim.xdim;
    yrange = 1:nim.ydim;
    zrange = slice:slice;
    size2 = nim.xdim * nim.ydim;
    hi = 1;
    vi = 2;
  elseif plane == "xz"
    xrange = 1:nim.xdim;
    yrange = slice:slice;
    zrange = 1:nim.zdim;
    size2 = nim.xdim * nim.zdim;
    hi = 1;
    vi = 3;
  elseif plane == "yz"
    xrange = slice:slice;
    yrange = 1:nim.ydim;
    zrange = 1:nim.zdim;
    size2 = nim.ydim * nim.zdim;
    hi = 2;
    vi = 3;
  end

  FA = zeros(size2, 4);
  i = 1;
  for x=xrange
    for y=yrange
      for z=zrange
        if enable_mask && nim.mask(x, y, z) == 0
          FA(i, :) = [ x y z 0 ];
        else
          FA(i, :) = [ x y z nim.FA(x, y, z) ];
        end
        i = i + 1;
      end
    end
  end

  fig = figure("Name", "Fractional Anisotropy (FA)");
  p = pcolor(reshape(FA(:, 4), [nim.hdr.ImageSize(vi) nim.hdr.ImageSize(hi)]));
  p.EdgeColor = "none";
  title("Fractional Anisotropy (FA)");
  subtitle(plane + " plane, z=" + slice + ", size=" + nim.hdr.ImageSize(hi) + "*" + nim.hdr.ImageSize(vi));
  colormap('gray');
  colorbar('eastoutside');
  axis equal;

end
