function [h,hinv] = psf_fdr(sizes, steps, starts, psf, psf_starts,psf_steps,psf_sizes);
%    compute the psf in FDR space where the latter has coordinates of spatial frequency and angular frequency

%   sizes  :  (Nx, Nphi)   dimensions of image
%   steps  :  (sx, sphi)   voxel size in x and l dimensions
%   starts :  (cx, cphi)   centre of 0,0 voxel
%   psf    :  ArrayVolume specifying psf

% alloc array
h = zeros(sizes(1), sizes(2));
zero_offset = size(psf,2)/2;

% scale the zoomed psf measurement according to eq 2 in the patent
% guessing values for the various parameters
NA=0.0155;
lamda=473e-9;
zoom=12.8/2.45;
nbath=1.33;
e=6.45e-6;
DOFmax=nbath*1.305*lamda/NA^2;
% stole this table of NA/zoom values from the sphagetti code. mw is saying 
% that eq 2 is incorrect, but also denies this NA data exists. however, 
% if this table is true then the scaling can be made to make sense. 
% note that this table is for the microscope either
% from the prototype or the bioptonics and we are using the values 
% on the custom opt. so while not exactly accurate, it is at least plausible. 
zoomedNA = [ 0.80  0.0125; 1.00  0.0155;
                    1.30  0.0180;
                    1.60  0.0225;
                    2.00  0.0270;
                    2.50  0.0320;
                    3.20  0.0390;
                    4.00  0.0465;
                    5.00  0.0505;
                    6.30  0.0620;
                    8.00  0.0625;
                    10.0  0.0625  ];
scaledNA=interp1(zoomedNA(:,1),zoomedNA(:,2),zoom);
lscale=DOFmax/(nbath*(lamda/scaledNA^2+e/zoom/scaledNA));
% since nobody is providing the position of the focal plane either, guess 
% a value as 1/2 the DOF of the experiment, converting from metres to mm
pfp=DOFmax/2 * 1000;

% compute fft of psf along detector element dimension (x)
Rx_psf = fftshift(fft(fftshift((psf(:,:)),2),[],2) ,2);
% calculate R (Fourier space) coord
Rx_values = ([-zero_offset+1:psf_sizes(2)-zero_offset]+0.5)/(psf_steps(2)*psf_sizes(2)); 
% calculate the l coord, scale it to the brain dataset, and offset by the pfp
l_values = (psf_starts(3) + [1:psf_sizes(3)]*psf_steps(3))*lscale + pfp;

figure(1),subplot(1,3,2),imagesc(Rx_values', l_values',abs(Rx_psf),[0 2]);
xlabel('R (1/mm)');
ylabel('distance from focal plane (mm)');

% generate matrices from the coordinate arrays for bilinear interp
[l_values,Rx_values]=meshgrid(l_values,Rx_values);

for j=1:sizes(2)  % for each angular frequency
  for i=1:sizes(1)  % for each detector spatial frequency

% these coords are calculated in Fourier space, not sinogram space
% it has to be in Fourier space in order to get units of 
% mm for l. 
    Phi(i,j) = ( - sizes(2)/2+j)/(steps(2)*sizes(2));
    Rx(i,j)  = (-sizes(1)/2+i+0.5)/(steps(1)*sizes(1));
    l(i,j) = -Phi(i,j)/Rx(i,j);
  end
end
lmax=DOFmax*1000*1.5;
lmin=-DOFmax*1000/1.5;

h = interp2(l_values, Rx_values,Rx_psf.',l, Rx,'linear');
hinv = 1./h;

% a filter to smooth off the sharp edges of the nullspace at the lmax/lmin 
% cutoff points, while also preserving a small part of the nullspace 
% (ie the "vee") very close to the origin where there looks to be more 
% than just noise, since we also know the FDR is only an approximation 
% to begin with.

% some hardcoded cutoff values for defining the region of origin preservation
% need to parameterize these later.
Pmin=-1; % 1/rad
Rmin=2; % 1/mm
% generate the cutoff line P(Rx) as a piecewise constant line in Phi,R space,
% and smooth the vertices a little bit
clear P;
for i=1:sizes(1) 
  P(i)=-Rx(i)*lmax ./(1+exp(-(Rx(i)-Rmin))) + Pmin / (1+exp(-(Rmin-Rx(i)))) + Pmin / (1+exp(Rx(i)+Rmin))   + -Rx(i)*lmin ./ (1+exp((Rx(i)+Rmin)));
end
% map that line to a Phi-based sigmoidal roll off function
P=repmat(P,[sizes(2) 1]);
sg=1./(1+exp(-(Phi-P')));
% generate correct values for the other half of Phi,R space by a 
% reflection symmetry argument.
sg(:,end/2+1:end)=flipud(fliplr(sg(:,1:end/2)));

% original slope-based rolloff filter
%w=2;
%sg2 = zeros(size(l));
%sg2(find(l<0))=1.0;
%sg2(find(l<w&l>0))=cos(pi/2*abs(l(find(l<w&l>0)))/w).^2;
%sg2(find(l>w))=0;

h = h.* sg ;
hinv = hinv .* sg;

figure(1),subplot(1,3,3),imagesc(Rx(:,1),Phi(1,:),(abs(h)'),[0 3]);
xlabel('R (1/mm)');
ylabel('\Phi (1/rad)');
title('psf in \Phi/R space for JoeS data');
