function [img]= simulate_2d_gaussian_psf(sizes, steps, starts, centre, width, depth_dependent_width);
%    """simulate the image of point source at location centre in
%    translated along dimension l and having a gaussian profile in x
%    with a minimum width specified by width and an optional l (or
%    depth dependence)

%    sizes :  (Nx, Nl)   dimensions of simulated image
%    steps :  (sx, sl)   voxel size in x and l dimensions
%    starts : (cx, cl)   centre of 0,0 voxel
%    width  :  width in x of gaussian
%    depth_dependent_width : additional width term which is proportional to l-cl and combines in quadrature with base width

% alloc array
img = struct('sizes',[1 sizes(1) sizes(2)],'dimension_names',{{'zspace', ...
                    'xspace','lspace'}},'nc_data_type','NC_SHORT', ...
             'typecode','float_','starts',[0 starts(1) starts(2)],'steps',[1  steps(1) steps(2)],'array',zeros(1,sizes(1),sizes(2)));
%img = newimage((1, sizes(1), sizes(2)), dimension_names = ('zspace', 'xspace', 'lspace'), nc_data_type = py_minc.NC_SHORT, typecode = float_);
%img.set_starts((0, starts(1), starts(2)));
%img.set_separations((1, steps(1), steps(2)));

% iterate over voxels
for j=1:sizes(2)
  l = starts(2) + j*steps(2);
  dl = l - centre(2); % coordinate relative to focal plane

  combined_width_sq = width^2 + (dl*depth_dependent_width)^2;
  c = 1.0/sqrt(2*pi*combined_width_sq);

  for i=1:sizes(1)
    x = starts(1) + i*steps(1);
    dx = x - centre(1);   % coordinate relative to centre
    
    img.array(1,i,j) = c*exp(-.5*dx^2/combined_width_sq); % evaluate gaussian
  end
end

