function Map_new=FillSmoothResize3D(Map,mask,mtx_new)
% properly fill, smooth and resize 3D map
Map(mask==0)=NaN;
Map(Map==0)=NaN;
Map1 = inpaintn(Map);
Map2=fiex3d(Map1,2,2);
[x, y, z] = meshgrid(1:size(Map2,1), 1:size(Map2,2), 1:size(Map2,3));
[xq, yq, zq] = meshgrid(linspace(1, size(Map2,1), mtx_new(1)), ...
                      linspace(1, size(Map2,2), mtx_new(2)), ...
                      linspace(1, size(Map2,3), mtx_new(3)));
Map_new = interp3(x, y, z, Map2, xq, yq, zq, 'spline');
end