% colormap_for_rr

function colormap_nico(ha,flag,lim)

switch nargin
    case 1; flag = 1; lim = get(ha,'clim');
    case 2; lim = get(ha,'clim');
end

switch flag
    case 1
    % colormap for rr or any 0+ field
    axes(ha);
    caxis(lim);
    cmap = colormap;
    cmap(1,:) = 1;
    colormap(cmap);
    colorbar;
end

return
