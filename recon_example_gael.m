addpath(genpath('C:\Users\gaelr\OneDrive\Documents\MATLAB\codes_4_gael\CODES\'))

%% Load data
datafolder = 'C:\Users\gaelr\OneDrive\Documents\MATLAB\codes_4_gael\DATA_PH\';
filename = 'ph_11072019_2035184_5_2_wip_phyllo_classicV4';

% load([datafolder filename '_MRobject_with_data'])
% kdata = squeeze(MR.Data);
load([datafolder filename '_channel_images_regrid'])
load([datafolder filename '_kdata'])
load([datafolder filename '_k'])
% load([datafolder filename '_dcf']) % obtenue par une autre methode que celle dans le code ci-dessous

[nx, ntviews, nz, nc] = size(kdata); % nc = number of coils

%% Load Coil sensitivities
% load([datafolder filename '_Sobject_recRes']) 
% csm = S.Sensitivity;
load([datafolder filename '_coils'])
csm = csm/max(csm(:));% normalize the coil map intensities
coil_rss = (sum(csm.*conj(csm),4));



%% compute the DCF density compensation function
kdataloc = reshape(kdata,[nx*ntviews*nz nc]);
k4dcf = reshape(k,[nx*ntviews*nz 3]);
k4dcf = permute(k4dcf,[3 2 1]);
verbose = 1; osf = 1.5; effMtx = 100; numIter = 10;
pre = sdc3_MAT(k4dcf,numIter,effMtx,verbose,osf);
osf = 2.1; numIter = 30;
DCF = sdc3_MAT(k4dcf,numIter,effMtx,verbose,osf,pre);
% DCF = permute(DCF,[2 1]);

%% choose recon resolution
oversampl = 1.25; % MR.Parameter.Encoding.KxOversampling;
acqFOV = [300 300 300]; % MR.Parameter.Scan.FOV; % actualFOV = acqFOV * 2; % because it is 3D radial 
acqRes = [2.2059 2.2059 2.2059]; % MR.Parameter.Scan.AcqVoxelSize;
acqImageSize = acqFOV./acqRes; % of the selected FOV, not the oversampled one
recRes = [2.0833 2.0833 2.0833]; % MR.Parameter.Scan.RecVoxelSize;
desiredRes = recRes; % [2 2 2];
recImageSize = uint32(round(oversampl*(acqFOV./desiredRes)));
recImageSize_nooversampl = uint32(round((acqFOV./desiredRes)));
k_factor_for_interp = desiredRes/acqRes;
k = k*k_factor_for_interp;

%% coil combination 
images = sum(channel_images .* conj(csm),4)./coil_rss; %imagine(images) % images from MRrecon
% images_sos = sum(channel_images .* conj(channel_images),4); imagine(images_sos); % sum of squares combination
S.('original_full_image') = images;
save('Images_n_spokes.mat','-struct','S')
clear S
clear image

%% do recon with CPU nufft
% my_DCF = reshape(DCF,size(kdatacpu(:,:,1)));
kdatacpu = reshape(kdata,[nx*ntviews*nz nc]); % samples
kcpu = reshape(k,[nx*ntviews*nz 3]); % position of the evaluation
DCF = reshape(DCF,[nx*ntviews*nz 1]); % compensate for the density of points in the k-space
precision = 1E-2;
% siz = size(coil_rss); % size of grid
siz = [180 180 180];
for coil = 1:nc
    k2im(:,:,:,coil) = nufft3_type1(double(kcpu), double(kdatacpu(:,coil).*DCF), siz, +1,precision);
end
% images = sum(k2im .* conj(csm),4)./coil_rss; % combinaison des images par antenne with antenna (imagine(csm) donne antenne map; imagine(k2m) pour voir les images des differentes antennes)
%% Project back into k-space with inverse nufft, same k distribution as a temoin
for coil = 1:nc
    im2k(:,coil) = nufft3_type2(double(kcpu), double(k2im(:,:,:,coil)),-1,precision);
    test_im(:,:,:,coil) = nufft3_type1(double(kcpu), double(im2k(:,coil).*DCF), siz, +1,precision);
end
test_images = sum(test_im .* conj(csm),4)./coil_rss;%imagine(test_images)
S.('original_reproj_image') = test_images;
save('Images_n_spokes.mat','-struct','S','-append')
clear S
clear test_images
clear im2k
clear test_im

%% create out file
% save('Images_n_spokes.mat','nc')

%% %% load my own k-space point distribution
% load([datafolder 'matlab_test_full'], 'phyllotaxis_8')
% 
% ReProjectAndSave(kcpu, 'Images_trajectories.mat', 'control_im',nc,k2im,csm,coil_rss)
% 
% ReProjectAndSave(phyllotaxis, 'Images_trajectories.mat', 'phyllotaxis',nc,k2im,csm,coil_rss)
% clear phyllotaxis
% ReProjectAndSave(phyllotaxis_4, 'Images_n_spokes.mat', 'phyllotaxis_4_cardiacphase',nc,k2im,csm,coil_rss)
% clear phyllotaxis_4
% 
% ReProjectAndSave(phyllotaxis_6, 'Images_n_spokes.mat', 'phyllotaxis_6_cardiacphase',nc,k2im,csm,coil_rss)
% clear phyllotaxis_6
% ReProjectAndSave(phyllotaxis_8, 'Images_n_spokes.mat', 'phyllotaxis_8',nc,k2im,csm,coil_rss)
% clear phyllotaxis_8
% load([datafolder 'matlab_test_gold'], 'phyllotaxis_8_monicasgold')
% ReProjectAndSave(phyllotaxis_8_monicasgold, 'Images_n_spokes.mat', 'phyllotaxis_8_monicasgold',nc,k2im,csm,coil_rss)
% clear phyllotaxis_8_monicasgold
% load([datafolder 'matlab_test_gold'], 'phyllotaxis_8_monicasgold_init1')
% ReProjectAndSave(phyllotaxis_8_monicasgold_init1, 'Images_n_spokes.mat', 'phyllotaxis_8_monicasgold_init1',nc,k2im,csm,coil_rss)
% clear phyllotaxis_8_monicasgold_init1
% load([datafolder 'matlab_test_gold'], 'phyllotaxis_8_monicasgold_init1_yx')
% ReProjectAndSave(phyllotaxis_8_monicasgold_init1_yx, 'Images_n_spokes.mat', 'phyllotaxis_8_monicasgold_init1_yx',nc,k2im,csm,coil_rss)
% clear phyllotaxis_8_monicasgold_init1_yx
% load([datafolder 'matlab_test_gold'], 'phyllotaxis_8_monicasgold_init1_yx_otheraltern')
% ReProjectAndSave(phyllotaxis_8_monicasgold_init1_yx_otheraltern, 'Images_n_spokes.mat', 'phyllotaxis_8_monicasgold_init1_yx_otheraltern',nc,k2im,csm,coil_rss)
% clear phyllotaxis_8_monicasgold_init1_yx_otheraltern
load([datafolder 'matlab_test_gold2'], 'phyllotaxis_8_gold24')
ReProjectAndSave(phyllotaxis_8_gold24, 'Images_n_spokes.mat', 'phyllotaxis_8_gold24',nc,k2im,csm,coil_rss)
clear phyllotaxis_8_gold24
load([datafolder 'matlab_test_gold2'], 'phyllotaxis_8_goldnot24')
ReProjectAndSave(phyllotaxis_8_goldnot24, 'Images_n_spokes.mat', 'phyllotaxis_8_goldnot24',nc,k2im,csm,coil_rss)
clear phyllotaxis_8_goldnot24
% ReProjectAndSave(phyllotaxis_8_notfib_phillipsreadouts, 'Images_n_spokes.mat', 'phyllotaxis_8_notfib_phillipsreadouts',nc,k2im,csm,coil_rss)
% clear phyllotaxis_8_notfib_phillipsreadouts
% ReProjectAndSave(phyllotaxis_10, 'Images_n_spokes.mat', 'phyllotaxis_10_cardiacphase',nc,k2im,csm,coil_rss)
% clear phyllotaxis_10
% ReProjectAndSave(phyllotaxis_12, 'Images_n_spokes.mat', 'phyllotaxis_12_cardiacphase',nc,k2im,csm,coil_rss)
% clear phyllotaxis_12
% ReProjectAndSave(phyllotaxis_14, 'Images_n_spokes.mat', 'phyllotaxis_14_cardiacphase',nc,k2im,csm,coil_rss)
% clear phyllotaxis_14
% ReProjectAndSave(phyllotaxis_16, 'Images_n_spokes.mat', 'phyllotaxis_16',nc,k2im,csm,coil_rss)
% clear phyllotaxis_16
% ReProjectAndSave(stack_of_stars, 'Images_trajectories.mat', 'stack_of_stars',nc,k2im,csm,coil_rss)
% clear stack_of_stars
% 
% ReProjectAndSave(multigolden, 'Images_trajectories.mat', 'multigolden',nc,k2im,csm,coil_rss)
% clear multigolden

% %% recompute dcf
% [nx, ntviews, nz, nco] = size(phyllotaxis);
% phyllotaxis4dcf = reshape(phyllotaxis,[nx*ntviews*nz 3]);
% phyllotaxis4dcf = permute(phyllotaxis4dcf,[3 2 1]);
% verbose = 1; osf = 1.5; effMtx = 100; numIter = 10;
% mypre = sdc3_MAT(double(phyllotaxis4dcf),numIter,effMtx,verbose,osf);
% osf = 2.1; numIter = 30;
% myDCF = sdc3_MAT(double(phyllotaxis4dcf),numIter,effMtx,verbose,osf,mypre);

% %% Project back into k-space with inverse nufft, with own k distribution
% siz = [180 180 180];
% for coil = 1:nc
%     myim2k(:,coil) = nufft3_type2(double(phyllotaxis), double(k2im(:,:,:,coil)),-1,precision);
%     mytest_im(:,:,:,coil) = nufft3_type1(double(phyllotaxis), double(myim2k(:,coil).*myDCF), siz, +1,precision);
% end
% my_test_images = sum(mytest_im .* conj(csm),4)./coil_rss; imagine(my_test_images)

% %% recompute dcf
% [nx, ntviews, nz, nco] = size(stack_of_stars);
% stack_of_stars4dcf = reshape(stack_of_stars,[nx*ntviews*nz 3]);
% stack_of_stars4dcf = permute(stack_of_stars4dcf,[3 2 1]);
% verbose = 1; osf = 1.5; effMtx = 100; numIter = 10;
% mysospre = sdc3_MAT(double(stack_of_stars4dcf),numIter,effMtx,verbose,osf);
% osf = 2.1; numIter = 30;
% mysosDCF = sdc3_MAT(double(stack_of_stars4dcf),numIter,effMtx,verbose,osf,mysospre);

% %% Project back into k-space with inverse nufft, with own k distribution
% for coil = 1:nc
%     mysosim2k(:,coil) = nufft3_type2(double(stack_of_stars), double(k2im(:,:,:,coil)),-1,precision);
%     test_sosim(:,:,:,coil) = nufft3_type1(double(stack_of_stars), double(mysosim2k(:,coil).*mysosDCF), siz, +1,precision);
% end
% my_test_sosimages = sum(test_sosim .* conj(csm),4)./coil_rss; imagine(my_test_sosimages)
% % 
% % %% %% recompute dcf
% % [nx, ntviews, nco] = size(multigolden);
% % mg4dcf = reshape(multigolden,[nx*ntviews 3]);
% % mg4dcf = permute(mg4dcf,[3 2 1]);
% % verbose = 1; osf = 1.5; effMtx = 100; numIter = 10;
% % mymgpre = sdc3_MAT(double(mg4dcf),numIter,effMtx,verbose,osf);
% % osf = 2.1; numIter = 30;
% % mymgDCF = sdc3_MAT(double(mg4dcf),numIter,effMtx,verbose,osf,mymgpre);
% % 
% % %% Project back into k-space with inverse nufft, with own k distribution
% % for coil = 1:nc
% %     mymgim2k(:,coil) = nufft3_type2(double(multigolden), double(k2im(:,:,:,coil)),-1,precision);
% %     test_mgim(:,:,:,coil) = nufft3_type1(double(multigolden), double(mymgim2k(:,coil).*mymgDCF), siz, +1,precision);
% % end
% % my_test_mgimages = sum(test_mgim .* conj(csm),4)./coil_rss; imagine(my_test_mgimages)

