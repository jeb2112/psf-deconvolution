function [limit_recovery] = limit_recovery_filter(sizes, steps, starts, fdr_inv, fdrlimit=5);
%    sizes  :  (Nx, Nphi)   dimensions of image
%    steps  :  (sx, sphi)   voxel size in x and l dimensions
%    starts :  (cx, cphi)   centre of 0,0 voxel
%    fdr_inv    :  inverse fdr filter

% alloc array
limit_recovery = zeros([sizes(1),sizes(2)], complex_);
%% defined in jwcode
Ct=fdrlimit;
Cr=fdrlimit/2;
%%defined in paper
%Ct=1e-3
%Cr=1e-4
for j=1:sizes(2)        
% iterate over detector elements
  for i=1:sizes(1)
    val1 = abs(fdr_inv[i,j]);
    if(val1 > Ct)
      val2 = Ct+Cr - Cr*exp(-(val1-Ct)/Cr);
      limit_recovery[i,j] *= val2/val1;
    else
      limit_recovery[i,j] = fdr_inv[i,j];
    end
  end
end

return
