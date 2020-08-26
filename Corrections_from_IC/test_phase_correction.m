addpath(genpath('C:\Users\gaelr\OneDrive\Documents\MATLAB\codes_4_gael\CODES\'))
datafolder = 'C:\Users\gaelr\OneDrive\Documents\MATLAB\codes_4_gael\DATA_PH\';
filename = 'ph_11072019_2035184_5_2_wip_phyllo_classicV4';
load([datafolder filename '_coils']);
csm = csm/max(csm(:));% normalize the coil map intensities
coil_rss = (sum(csm.*conj(csm),4));
load([datafolder filename '_kdata']); load([datafolder filename '_channel_images_regrid']); 
load([datafolder filename '_k']); load([datafolder filename '_dcf'])
kdata_phyllo = kdata; k_phyllo = k; dcf_phyllo = dcf; chims_phyllo = squeeze(channel_images);
[nx,ntviews,nz,nc,npc] = size(kdata_phyllo);

corr = linear_phase_correction(kdata_phyllo,3,0,1)