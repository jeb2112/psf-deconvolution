function [roll_off] = roll_off_filter(sizes, steps, starts, lmax, lmin);
%sizes  :  (Nx, Nphi)   dimensions of image
%steps  :  (sx, sphi)   voxel size in x and l dimensions
%starts :  (cx, cphi)   centre of 0,0 voxel

roll_off = zeros(sizes(1),sizes(2));
% evaluate each voxel
for j =1:sizes(2) % for each angular frequency
  for i=1:sizes(1)   % for each detector spatial frequency
    Phi(i,j) = ( - sizes(2)/2+j)/(steps(2)*sizes(2));
    Rx(i,j)  = (-sizes(1)/2+i+0.5)/(steps(1)*sizes(1));
    
      l(i,j) = -Phi(i,j)/Rx(i,j); 
      if (l(i,j) > lmax+lmin)
        roll_off(i,j) = 0.0;
      elseif (l(i,j) <=lmax+lmin & l(i,j) > lmax-lmin)
        roll_off(i,j) = cos(pi/2.0*abs(lmax-lmin-l(i,j))/(2*lmin));
      elseif (l(i,j) <= lmax & l(i,j) > lmin)
        roll_off(i,j) = 1.0;
      elseif ( l(i,j) > -lmin & l(i,j) < lmin)
        roll_off(i,j) = cos(pi/2.0*abs(lmin-l(i,j))/(2*lmin));
      else
        roll_off(i,j) = 0.0;
      end
  end
end
