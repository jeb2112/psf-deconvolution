% Load psf 
infile = openimage('psf512_order.mnc');
psf = getimages(infile,1);
psf = psf/max(psf(:));
psf_sizes = getimageinfo(infile,'dimsizes');
psf_sizes = psf_sizes(2:4);
psf=reshape(psf,[psf_sizes(3) psf_sizes(2)]);
psf_steps = getimageinfo(infile,'steps');
psf_starts = getimageinfo(infile,'starts');
% psf data are actually 2.45 um pixels due to zoom factor 12307,
% the scaling for zoom is moved to psf_in_fdr_space.m
psf_steps=[0.00245 0.00245 0.00245];
psf_starts = [0 -512*psf_steps(2) -384.71*psf_steps(3)];
closeimage(infile);

figure(1),subplot(1,3,1),imagesc([psf_starts(2) ...
                    psf_starts(2)+(psf_sizes(2)-1)*psf_steps(2)],[psf_starts(3), psf_starts(3)+(psf_sizes(3)-1)*psf_steps(3)],psf(:,:),[0 1]);
ylabel('distance from focal plane (mm)');
xlabel('detector position (mm)');
title('psf in l,r space');

%load 3d psf fft
load3d=0;
if load3d
  psf3d=zeros(785,1024,1024,'int16');
  infile = openimage('psf_fft_mags.mnc');
  for i=1:785
    tmp = getimages(infile,i);
    psf3d(i,:,:) = reshape(int16(tmp/max(tmp(:))*32767),[1024 1024]);
  end
  psf3d_sizes = getimageinfo(infile,'dimsizes');
  psf3d_sizes = psf3d_sizes(2:4);
  psf3d_steps = getimageinfo(infile,'steps');
  psf3d_starts = getimageinfo(infile,'starts');
% psf data are actually 2.45 um pixels due to zoom factor 12307,
% the scaling for zoom is moved to psf_in_fdr_space.m
%  psf3d_steps=[0.00245 0.00245 0.00245];
%  psf3d_starts = [0 -512*psf_steps(2) -384.71*psf_steps(3)];
  closeimage(infile);
end

% Load data
% Joes brain
%infilename='slice512_order.mnc';
%infile = openimage(infilename);zoomdata=1;zoompsf=5.22;cor=0;
% mw embryo
infilename='embryo2.mnc';
infile = openimage(infilename);zoomdata=2.871;zoompsf=5.22;cor=3
% ga
%infilename='ga_slice512_order.mnc';
%infile = openimage(infilename);zoomdata=2.038;zoompsf=5.22;cor=16;
sino = getimages(infile,1);
sizes = getimageinfo(infile,'dimsizes');
sizes=sizes(3:4);

steps = getimageinfo(infile,'steps');
% if the step values in the minc arent consistent with the zoom values?
%steps = steps * zoompsf/zoomdata;
steps = [steps(2) 2*pi/sizes(2)];
xsteps=steps(1);
starts = getimageinfo(infile,'starts');
starts = [starts(2) 0];
xstarts = starts(1);
closeimage(infile);

% note on the data matrices in these calculations... the sizes dimensions 
% are a bit mixed up in this code so there are a lot of transposes (the .' or '
% operators in matlab) going on here and there. just make sure that r space 
% is the 1024 size and phi space is the 1200 size
sino=reshape(sino,[sizes(2) sizes(1)]);
% kludge for mw embryo and ga data, not used in joes data
sino = fliplr(sino);
figure(2),subplot(2,2,1),imagesc([(starts(1)-sizes(1)/2)*steps(1)  starts(1)+(sizes(1)/2-1)*steps(1)], [starts(2), starts(2)+(sizes(2)-1)*steps(2)],sino,[0 25]);
ylabel('angle (rad)');
xlabel('detector position (mm)');
title('original sinogram');
S=fftshift(fft2(fftshift(sino)));
figure(2),subplot(2,2,2),imagesc([-sizes(1)/2:sizes(1)/2-1]/(steps(1)*sizes(1)),[-sizes(2)/2:sizes(2)/2-1]/(steps(2)*sizes(2)),abs(S),[0 50000]);
xlabel('R (1/mm)');
ylabel(' \Phi (1/rad)');
title('original Fourier space');

%start of the Wiener filter. draw an roi for noise power measurement
if ~exist('nmask','var')
  disp('Draw ROI in noise area of Phi-R space in fig 4');
  figure(4),clf,imagesc(abs(S),[0 10000]);
  nmask=roipoly;
  delete(gcf);
else
  disp('Not updating nmask');
end
% within the noise ROI, calculate the mean square noise power
n=std(S(find(nmask)));
if ~exist('nscale','var')
  nscale=1;
else
  disp(sprintf('nscale = %d',nscale));
end
N=abs(n)^2;
% estimate of the signal power spectrum is obtained by subtracting
% the estimated noise power spectrum
s=abs(S.^2-N);
% some arbitrary smoothing of the signal power spectrum so that
% the it does not contribute noise
[l,m]=meshgrid(-21:21,-21:21);
sig=11;
lofil=exp(-(l.^2+m.^2)/sig^2);
lofil=lofil/sum(lofil(:));
s=conv2(s,lofil,'same');
% plot the signal/noise power ratio
figure(3),subplot(1,4,1),imagesc([-sizes(1)/2:sizes(1)/2-1]/(steps(1)*sizes(1)), ...
                          [-sizes(2)/2:sizes(2)/2-1]/(steps(2)*sizes(2)),log10(s./(N)),[0 3]);
xlabel('R (1/mm)');ylabel('\Phi (1/rad)');title('Wiener filter');



% mapping the psf to the Phi-R space of the imaging data
if ~exist('h','var')
  [h,hinv] = psf_in_fdr_space(sizes, steps, starts, zoomdata, psf, psf_starts, ...
                              psf_steps,psf_sizes,zoompsf);
else
  disp('Not updating h,hinv');
end
figure(1),subplot(1,3,3),imagesc([-sizes(1)/2:sizes(1)/2-1]/(steps(1)*sizes(1)), ...
                          [-sizes(2)/2:sizes(2)/2-1]/(steps(2)*sizes(2)),(abs(h)'),[0 1]);
xlabel('R (1/mm)');
ylabel('\Phi (1/rad)');
title('psf in \Phi/R space');
% take the inverse. this is 1/OTF in king
figure(3),subplot(1,4,2),imagesc([-sizes(1)/2:sizes(1)/2-1]/(steps(1)*sizes(1)), ...
                          [-sizes(2)/2:sizes(2)/2-1]/(steps(2)*sizes(2)),abs(hinv.'),[0 max(abs(hinv(:)))/1000]);
xlabel('R (1/mm)');ylabel('\Phi (1/rad)');title('Inverse PSF');

%
% form three alternative reconstructions to compare the value of using
% the psf to roll off the restoration filter and thereby limit noise 
% amplification b
%

% Wiener filter, with no additional rolloff
LPF0 = 1./(1 + (N*nscale)./(s));
%figure(3),subplot(1,4,3),imagesc([-sizes(1)/2:sizes(1)/2-1]/(steps(1)*sizes(1)), ...   [-sizes(2)/2:sizes(2)/2-1]/(steps(2)*sizes(2)),abs(LPF0),[0 1]);
%xlabel('R (1/mm)');ylabel('\Phi (1/rad)');title('Wiener filter');

% use of 1/h for the rolloff function
LPF1 = abs(h.')./(abs(h.') + (N)./(s));
%figure(3),subplot(1,4,3),imagesc([-sizes(1)/2:sizes(1)/2-1]/(steps(1)*sizes(1)), ...  [-sizes(2)/2:sizes(2)/2-1]/(steps(2)*sizes(2)),abs(LPF1),[0 1]);
%xlabel('R (1/mm)');ylabel('\Phi (1/rad)');title('Low-pass filter, 1/h rolloff');

% use of 1/h^2 for the rolloff function
LPF2 = abs(h.').^2./(abs(h.').^2 + (N)./(s));
figure(3),subplot(1,4,3),imagesc([-sizes(1)/2:sizes(1)/2-1]/(steps(1)*sizes(1)), ...
                          [-sizes(2)/2:sizes(2)/2-1]/(steps(2)*sizes(2)),abs(LPF2),[0 1]);
xlabel('R (1/mm)');ylabel('\Phi (1/rad)');title('Low-pass filter, 1/h^2 rolloff');


% the combined filter
hW0 = hinv.'.*LPF0;
hW1 = hinv.'.*LPF1;
hW2 = hinv.'.*LPF2;

% providing a normalization here as it has been lacking in the previous
% steps
hW0 = hW0/abs(hW0(sizes(2)/2,sizes(1)/2));
hW1 = hW1/abs(hW1(sizes(2)/2,sizes(1)/2));
hW2 = hW2/abs(hW2(sizes(2)/2,sizes(1)/2));

% provide a local interpolation to fill in and/or demphasize Fourier
% space near the origin since the FDR is only an approximation and
% does not govern these data. currently the local region of Fourier
% space is hardcoded to the 60 central lines of Phi space, and +/10
% points on either side of the vee in Rx space. This region can be
% parameterized as a percentage of the space (about 5%) in future.
for j=sizes(2)/2+30:-1:sizes(2)/2-30
  [jset,iset]=find(isnan(hW0(j,:)));
  if (isempty(iset))
    iset=sizes(1)/2;
  end
  iset=[min(iset)-10:min(iset)-1 iset max(iset)+1:max(iset)+10];
  frange=[min(iset)-20:min(iset)-1 max(iset)+1:max(iset)+20];
  hW0(j,iset)=interp1(frange,hW0(j,frange),iset,'spline');
  hW1(j,iset)=interp1(frange,hW1(j,frange),iset,'spline');
  hW2(j,iset)=interp1(frange,hW2(j,frange),iset,'spline');
end

% get rid of NaNs in the nullspace, as they are not part of the problem
h(find(isnan(h)))=0;
hW0(find(isnan(hW0)))=0;
hW1(find(isnan(hW1)))=0;
hW2(find(isnan(hW2)))=0;

figure(3),subplot(1,4,4),imagesc([-sizes(1)/2:sizes(1)/2-1]/(steps(1)*sizes(1)), ...
                          [-sizes(2)/2:sizes(2)/2-1]/(steps(2)*sizes(2)),abs(hW2),[0 5]);
xlabel('R (1/mm)');ylabel('\Phi (1/rad)');title('1/h^2 rolloff');


%applying complete filter
rest_S = S .* hW2;
rest_S(find(isnan(rest_S)))=0;
figure(2),subplot(2,2,4),imagesc([-sizes(1)/2:sizes(1)/2-1]/(steps(1)*sizes(1)),[-sizes(2)/2:sizes(2)/2-1]/(steps(2)*sizes(2)),abs(rest_S),[0 50000]);
xlabel('R (1/mm)');
ylabel(' \Phi (1/rad)');
title('restored Fourier space');

% fft back to sinogram space
rest_sino = fftshift(ifft2(fftshift(rest_S)));
figure(2),subplot(2,2,3),imagesc([(starts(1)-sizes(1)/2)*steps(1)  starts(1)+(sizes(1)/2-1)*steps(1)], [starts(2), starts(2)+(sizes(2)-1)*steps(2)],real(rest_sino),[0 250]);
ylabel('angle (rad)');
xlabel('detector position (mm)');
title('restored sinogram');
% minc sinogram output
delete('restored_sinogram.mnc');
delete('restored_sinogram.raw');
%rest_sino=rest_sino.'; use this on joes brain but not the embryo
%outfile = newimage('rsino.mnc', [0 1 sizes(2) sizes(1)],infilename,'float','transverse');
%putimages(outfile,real(rest_sino(:)),1);
%closeimage(outfile);
fid=fopen('restored_sinogram.raw','wb','l');
count=fwrite(fid,abs(rest_sino(:)),'float32');
fclose(fid);
unix(['rawtominc -input restored_sinogram.raw -float -clobber -xyz -xstart  '  num2str(xstarts) ' -xstep '  num2str(xsteps) ' -ystart ' num2str(starts(1)) ' -ystep ' num2str(steps(1)) ' -zstart ' num2str(starts(2)) ' -zstep ' num2str(steps(2)) ' restored_sinogram.mnc 1 ' num2str(sizes(1)) ' ' num2str(sizes(2))]);

delete('original_sinogram.mnc');
%sino_out=sino'; %for brain but not embryo
sino_out=sino;
%outfile = newimage('original_sinogram.mnc', [0 1 sizes(2) sizes(1)],infilename,'float');
%putimages(outfile,sino_out(:),1);
%closeimage(outfile);
fid=fopen('original_sinogram.raw','wb','l');
count=fwrite(fid,abs(sino_out(:)),'float32');
fclose(fid);
unix(['rawtominc -input original_sinogram.raw -float -clobber -xyz -xstart  '  num2str(xstarts) ' -xstep '  num2str(xsteps) ' -ystart ' num2str(starts(1)) ' -ystep ' num2str(steps(1)) ' -zstart ' num2str(starts(2)) ' -zstep ' num2str(steps(2)) ' original_sinogram.mnc 1 ' num2str(sizes(1)) ' ' num2str(sizes(2))]);

% transform to the image domain

% original image
Io=iradon(circshift(sino.',[cor 0]),[1:sizes(2)]*360/sizes(2),sizes(1)); 
delete('original_image.mnc');
%outfile = newimage('original_image.mnc', [0 1 sizes(1) sizes(1)],infilename,'float');
%putimages(outfile,abs(Io(:)),1);
%closeimage(outfile);
fid=fopen('original_image.raw','wb','l');
count=fwrite(fid,abs(Io(:)),'float32');
fclose(fid);
unix(['rawtominc -input original_image.raw -float -clobber -xyz -xstart  '  num2str(xstarts) ' -xstep '  num2str(xsteps) ' -ystart ' num2str(starts(1)) ' -ystep ' num2str(steps(1)) ' -zstart ' num2str(starts(2)) ' -zstep ' num2str(steps(2)) ' original_image.mnc 1 ' num2str(sizes(1)) ' ' num2str(sizes(1))]);

% 1/h^2 image
I2=iradon(circshift(rest_sino.',[cor 0]),[1:sizes(2)]*360/sizes(2),sizes(1));

% 1/h and wiener images
rest_S = S .* hW1;
rest_sino = fftshift(ifft2(fftshift(rest_S)));
I1=iradon(circshift(rest_sino.',[cor 0]),[1:sizes(2)]*360/sizes(2),sizes(1));
rest_S = S .* hW0;
rest_sino = fftshift(ifft2(fftshift(rest_S)));
I0=iradon(circshift(rest_sino.',[cor 0]),[1:sizes(2)]*360/sizes(2),sizes(1));

%f4,subplot(1,2,1),imagesc([(starts(1)-sizes(1)/2)*steps(1)  starts(1)+(sizes(1)/2-1)*steps(1)],[(starts(1)-sizes(1)/2)*steps(1)  starts(1)+(sizes(1)/2-1)*steps(1)],abs(I0),[0 1]);set(gca,'dataaspect',[1 1 1]);
%title('original');;xlabel('x1 (mm)');ylabel('x2 (mm)');
%f4,subplot(1,2,2),imagesc([(starts(1)-sizes(1)/2)*steps(1)  starts(1)+(sizes(1)/2-1)*steps(1)],[(starts(1)-sizes(1)/2)*steps(1)  starts(1)+(sizes(1)/2-1)*steps(1)],abs(I),[0 1]);set(gca,'dataaspect',[1 1 1]);colormap('gray');
%title('restored');xlabel('x1 (mm)');ylabel('x2 (mm)');

delete(['restored_image_2.mnc']);
delete(['restored_image_2.raw']);
fid=fopen('restored_image_2.raw','wb','l');
count=fwrite(fid,abs(I2(:)),'float32');
fclose(fid);
unix(['rawtominc -input restored_image_2.raw -float -clobber -xyz -xstart  '  num2str(xstarts) ' -xstep '  num2str(xsteps) ' -ystart ' num2str(starts(1)) ' -ystep ' num2str(steps(1)) ' -zstart ' num2str(starts(2)) ' -zstep ' num2str(steps(2)) ' restored_image_2.mnc 1 ' num2str(sizes(1)) ' ' num2str(sizes(1))]);
%outfile = newimage(['restored_image_2.mnc'], [0 1 sizes(1) sizes(1)],infilename,'float');
%putimages(outfile,abs(I2(:)),1);
%closeimage(outfile);

delete(['restored_image_1.mnc']);
fid=fopen('restored_image_1.raw','wb','l');
count=fwrite(fid,abs(I1(:)),'float32');
fclose(fid);
unix(['rawtominc -input restored_image_1.raw -float -clobber -xyz -xstart  '  num2str(xstarts) ' -xstep '  num2str(xsteps) ' -ystart ' num2str(starts(1)) ' -ystep ' num2str(steps(1)) ' -zstart ' num2str(starts(2)) ' -zstep ' num2str(steps(2)) ' restored_image_1.mnc 1 ' num2str(sizes(1)) ' ' num2str(sizes(1))]);
%outfile = newimage(['restored_image_1.mnc'], [0 1 sizes(1) sizes(1)],infilename,'float');
%putimages(outfile,abs(I1(:)),1);
%closeimage(outfile);

delete(['restored_image_0.mnc']);
fid=fopen('restored_image_0.raw','wb','l');
count=fwrite(fid,abs(I0(:)),'float32');
fclose(fid);
unix(['rawtominc -input restored_image_0.raw -float -clobber -xyz -xstart  '  num2str(xstarts) ' -xstep '  num2str(xsteps) ' -ystart ' num2str(starts(1)) ' -ystep ' num2str(steps(1)) ' -zstart ' num2str(starts(2)) ' -zstep ' num2str(steps(2)) ' restored_image_0.mnc 1 ' num2str(sizes(1)) ' ' num2str(sizes(1))]);
%outfile = newimage(['restored_image_0.mnc'], [0 1 sizes(1) sizes(1)],infilename,'float');
%putimages(outfile,abs(I0(:)),1);
%closeimage(outfile);


% comparison of wiener, 1/h, 1/h^2 rolloff filters based on h and SNR
hp=10.^(-2:.01:0);
SNR=10.^(-1:.1:2);
[hp,SNR]=meshgrid(hp,SNR);
h2=(hp.^2./(hp.^2+1./SNR.^2))./hp;
h1=(hp./(hp+1./SNR.^2))./hp;
h0=1./(1+1./SNR.^2) .* 1./hp;
cvec=[1 2 5 10 20 50 100];
figure(5),clf,subplot(1,3,1),[cs,hh]=contour(h0,cvec,'b');
clabel(cs,hh,'labelspacing',1000,'fontsize',15);
set(gca,'xtick',[1 100.5 201],'xticklabel',[.01 .1 1], 'ytick',[1:10:31],...
        'yticklabel',[.1 1 10 100],'fontsize',15);xlabel('H','fontsize',20);ylabel('SNR','fontsize',20); ...
    title('Amplification: Wiener filter','fontsize',15);
subplot(1,3,2),[cs,hh]=contour(h1,cvec,'b');
clabel(cs,hh,'labelspacing',1000,'fontsize',15);
set(gca,'xtick',[1 100.5 201],'xticklabel',[.01 .1 1], 'ytick',[1:10:31],...
        'yticklabel',[.1 1 10 100],'fontsize',15);xlabel('H','fontsize',20);ylabel('SNR','fontsize',20); ...
    title('1/H roll-off','fontsize',15);
subplot(1,3,3),[cs,hh]=contour(h2,cvec,'b');
clabel(cs,hh,'labelspacing',1000,'fontsize',15);
set(gca,'xtick',[1 100.5 201],'xticklabel',[.01 .1 1], 'ytick',[1:10:31],...
        'yticklabel',[.1 1 10 100],'fontsize',15);xlabel('H','fontsize',20);ylabel('SNR','fontsize',20); ...
    title('1/H^2 roll-off','fontsize',15);



