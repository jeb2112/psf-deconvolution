function [band_limit] = band_limit_filter(sizes, steps, starts, bandlimit);
%    """compute band limit filter 
%    sizes  :  (Nx, Nphi)   dimensions of image
%    steps  :  (sx, sphi)   voxel size in x and l dimensions
%    starts :  (cx, cphi)   centre of 0,0 voxel

% alloc array
band_limit = zeros(sizes(1),sizes(2));
% find the right bandlimit (so the code says!!)
bw_i = bandlimit*(sizes(1))/(4.5*10);

for j=1:sizes(2)  % for each angular frequency
  for i=1:sizes(1)  % for each detector spatial frequency
%            phi = starts(2) + j*steps(2)
    phi = starts(2) - sizes(2)/2*steps(2) + j*steps(2);
%          rx = starts(1) + i*steps(1)
    rx  = (-sizes(1)/2+i)/(steps(1)*sizes(1));
    if ( abs(rx) > bw_i)
      band_limit(i,j) = 0.0;
    elseif ( abs(rx) < 0.9*bw_i)
      band_limit(i,j) = 1.0;
    elseif ( abs(rx) > 0.9*bw_i & abs(rx) < bw_i)
      band_limit(i,j) = cos(pi/2.0*(rx-0.9*bw_i)/(0.1*bw_i))^2;
    end
  end
end       
