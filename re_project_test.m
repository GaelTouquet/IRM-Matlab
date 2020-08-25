addpath(genpath('C:\Users\gaelr\OneDrive\Documents\MATLAB\codes_4_gael\CODES\'))


%% Load data
datafolder = 'C:\Users\gaelr\OneDrive\Documents\MATLAB\codes_4_gael\DATA_PH\';
filename = 'ph_11072019_2035184_5_2_wip_phyllo_classicV4';

load([datafolder filename '_channel_images_regrid'])
load([datafolder filename '_kdata'])
load([datafolder filename '_k'])
load([datafolder filename '_dcf']) % obtenue par une autre methode que celle dans le code ci-dessous

[nx, ntviews, nz, nc] = size(kdata); % nc = number of coils
%% %% Load Coil sensitivities
% load([datafolder filename '_Sobject_recRes']) 
% csm = S.Sensitivity;
load([datafolder filename '_coils'])
csm = csm/max(csm(:));% normalize the coil map intensities
coil_rss = (sum(csm.*conj(csm),4));

%% do recon with CPU nufft
kdatacpu = reshape(kdata,[nx*ntviews*nz nc]); % samples
kcpu = reshape(k,[nx*ntviews*nz 3]); % position of the evaluation
dcf = reshape(dcf,[nx*ntviews*nz 1]);
precision = 1E-2;
siz = size(coil_rss); % size of grid
for coil = 1:nc
    k2im(:,:,:,coil) = nufft3_type1(double(kcpu), double(kdatacpu(:,coil).*dcf), siz, +1,precision);
end
images = sum(k2im .* conj(csm),4)./coil_rss;
%% load my own k-space point distribution
%%load([datafolder 'matlab_test'])

%% Project back into k-space with inverse nufft, only one coil
for coil = 1:nc
    im2k(:,coil) = nufft3_type2(double(kcpu), double(k2im(:,:,:,coil)),+1,precision);
    test_im(:,:,:,coil) = nufft3_type1(double(kcpu), double(im2k(:,coil).*dcf), siz, +1,precision);
end
test_images = sum(test_im .* conj(csm),4)./coil_rss;
imagine(test_images)
% imagine(images)