function plot_nim_interp(grid, indx, indy, indz)


loadfile = 'nim_interp_grid' + string(grid);
load(loadfile)

XXp = XX(:,:,:,indx,indy,indz);
YYp = YY(:,:,:,indx,indy,indz);
ZZp = ZZ(:,:,:,indx,indy,indz);

Vxp = Vx(:,:,:,indx,indy,indz);
Vyp = Vy(:,:,:,indx,indy,indz);
Vzp = Vz(:,:,:,indx,indy,indz);

figure(1)
quiver3(XXp, YYp, ZZp, Vxp, Vyp, Vzp, '-k', 'LineWidth', 2);
axis equal

end