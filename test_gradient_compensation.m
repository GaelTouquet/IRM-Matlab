addpath(genpath('.\CODES\'))
%% define data
datafolder = 'D:\DATA\';
filename = 'ph_11072019_2035184_5_2_wip_phyllo_classicV4';
load([datafolder filename '_coils']);
csm = csm/max(csm(:));% normalize the coil map intensities
coil_rss = (sum(csm.*conj(csm),4));
load([datafolder filename '_kdata']); load([datafolder filename '_channel_images_regrid']); 
load([datafolder filename '_k']); load([datafolder filename '_dcf'])

[nx,ntviews,nz,nc,npc] = size(kdata);
k = reshape(k,[nx*ntviews*nz 3]);
dcf = reshape(dcf,[nx*ntviews*nz 1]);
precision = 1E-2;
siz = [180 180 180];
%% 

kdata_corr = gradient_delay_compensation_HU( kdata );


%% 

rawdata_corr = reshape(kdata_corr,[nx*ntviews*nz nc]);
k2im_corr = zeros(size(siz),nc);
for coil = 1:nc
    k2im_corr(:,:,:,coil) = nufft3_type1(double(k), double(rawdata_corr(:,coil).*dcf), siz, +1,precision);
end
images_corr = sum(k2im .* conj(csm),4)./coil_rss;
%% 


rawdata = reshape(kdata,[nx*ntviews*nz nc]);
k2im = zeros(size(siz),nc);
for coil = 1:nc
    k2im(:,:,:,coil) = nufft3_type1(double(k), double(rawdata(:,coil).*dcf), siz, +1,precision);
end
images = sum(k2im .* conj(csm),4)./coil_rss;