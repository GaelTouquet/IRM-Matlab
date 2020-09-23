% addpath('D:\Matlab_functions')
% addpath(('D:\Projects_D\SequenceDevelopement\CODES_2018\'))
addpath(genpath('.\CODES\'))
% addpath(('D:\Projects_D\SequenceDevelopement\CODES_2018\Code_bioeng022\'))
% addpath(('D:\Projects_D\SequenceDevelopement\CODES_2018\com_computation\'))
% addpath(genpath('D:\Projects_D\SequenceDevelopement\CODES_2018\MRecon-3.0.553\'))
% addpath(genpath('D:\Projects_D\SequenceDevelopement\CODES_2018\MReconKCL\'))
% addpath(genpath('D:\Projects_D\SequenceDevelopement\CODES_2018\imagine\'))
% addpath(genpath('D:\Projects_D\SequenceDevelopement\CODES_2018\gpuNUFFT-master'))

%% define data
datafolder = 'D:\Work\Post-doc MRI\DATA\';

filename = 'ph_11072019_2035184_5_2_wip_phyllo_classicV4';
load([datafolder filename '_coils']);
csm = csm/max(csm(:));% normalize the coil map intensities
coil_rss = (sum(csm.*conj(csm),4));
load([datafolder filename '_kdata']); load([datafolder filename '_channel_images_regrid']); 
load([datafolder filename '_k']); load([datafolder filename '_dcf'])
kdata_phyllo = kdata; k_phyllo = k; dcf_phyllo = dcf; chims_phyllo = squeeze(channel_images);
[nx,ntviews,nz,nc,npc] = size(kdata_phyllo);


k = reshape(k_phyllo,[nx*ntviews*nz 3]);
dcf = reshape(dcf_phyllo,[nx*ntviews*nz 1]);
precision = 1E-2;
siz = [180 180 180];

%% compute the phase of the central point of k-space
Display_phase(kdata_phyllo,0, false);

%% try the error model for each coil
[kdata_corr_global, coeffs3D_global] = Phase_correction(kdata_phyllo, false, true);
Display_phase(kdata_corr_global,0,false);

%% try the error model for each spoke group of each coil
[kdata_corr_spokes, coeffs3D_spokes] = Phase_correction(kdata_phyllo, true, true);
Display_phase(kdata_corr_spokes,0,false);

%% Image without correction

rawdata = reshape(kdata_phyllo,[nx*ntviews*nz nc]);
for coil = 1:nc
    k2im_uncorr(:,:,:,coil) = nufft3_type1(double(k), double(rawdata(:,coil).*dcf), siz, +1,precision);
end
uncorr_images = sum(k2im_uncorr .* conj(csm),4)./coil_rss;

%% Image with global correction

rawdata = reshape(kdata_corr_global,[nx*ntviews*nz nc]);
for coil = 1:nc
    k2im_global(:,:,:,coil) = nufft3_type1(double(k), double(rawdata(:,coil).*dcf), siz, +1,precision);
end
corr_global_images = sum(k2im_global .* conj(csm),4)./coil_rss;

%% Image with spoke-based correction

rawdata = reshape(kdata_corr_spokes,[nx*ntviews*nz nc]);
for coil = 1:nc
    k2im_spokes(:,:,:,coil) = nufft3_type1(double(k), double(rawdata(:,coil).*dcf), siz, +1,precision);
end
corr_spokes_images = sum(k2im_spokes .* conj(csm),4)./coil_rss;

%% Try Giese
figure;
Y = reshape(reshape(abs(central_amp(:,:,:,1)),[nx nz*ntviews 1]),[nx*nz*ntviews 1]);
x1 = reshape(repmat(permute(rem(reshape(phi_phyllo(:,:),[ntviews*ninterleaf 1]),pi*2),[3 1 2]),[272 1 1]),[nx*nz*ntviews 1]);
x2 = reshape(repmat(permute(rem(reshape(theta_phyllo(:,:),[ninterleaf*ntviews 1]),pi*2),[3 1 2]),[272 1 1]),[nx*nz*ntviews 1]);
scatter3(x1,x2,Y,0.5)
hold on;
% figure;
Y = reshape(reshape(abs(central_amp(136,:,:,1)),[nz*ntviews 1]),[nz*ntviews 1]);
x1 = reshape(repmat(permute(rem(reshape(phi_phyllo(:,:),[ninterleaf*ntviews 1]),pi*2),[3 1 2]),[1 1 1]),[nz*ntviews 1]);
x2 = reshape(repmat(permute(rem(reshape(theta_phyllo(:,:),[ninterleaf*ntviews 1]),pi*2),[3 1 2]),[1 1 1]),[nz*ntviews 1]);
scatter3(x1,x2,Y,1,'r')
hold on;
Y = reshape(reshape(abs(central_amp(135,:,:,1)),[nz*ntviews 1]),[nz*ntviews 1]);
x1 = reshape(repmat(permute(rem(reshape(phi_phyllo(:,:),[ninterleaf*ntviews 1]),pi*2),[3 1 2]),[1 1 1]),[nz*ntviews 1]);
x2 = reshape(repmat(permute(rem(reshape(theta_phyllo(:,:),[ninterleaf*ntviews 1]),pi*2),[3 1 2]),[1 1 1]),[nz*ntviews 1]);
scatter3(x1,x2,Y,1,'g')
hold off;
%% Backup


% corrfact = permute(corrfact,[3 1 2]);
% corrfact = repmat(corrfact,[272 1 1]);
% Scorr(:,ileaf,:) = kdata_phyllo(:,ileaf,:,coil).*exp(-i*corrfact);

% corrected image

% %% try the error model for each interleaf
% for ileaf = 1:8
%     x1 = phi_phyllo(ileaf,:); x2 = theta_phyllo(ileaf,:); x = [x1(:),x2(:)];
%     Y = squeeze(phase_phyllo(:,ileaf,:,1)); Y = Y(:);
%     
%     errmodel3D = @(coeffs,x) coeffs(1).*cos(x(:,1)).*sin(x(:,2)) + ...
%     coeffs(2).*sin(x(:,1)).*sin(x(:,2)) + coeffs(3).*cos(x(:,2)) + coeffs(4);
% 
%     opts = statset('nlinfit');
%     opts.RobustWgtFun = 'andrews';
%     opts.TolFun = 10^-12;
%     opts.MaxFunEvals = 10^5;
%     [coeffs3D(:,ileaf),R,J,CovB,MSE] = nlinfit(x,Y(:),errmodel3D,[0,0,0,0],opts);
%     
%     Yfit = errmodel3D(coeffs3D(:,ileaf),x);
% 
%     figure;scatter3(x(:,1),x(:,2),Y);
%     hold on
%     scatter3(x(:,1),x(:,2),Yfit,'*k');
%       
%     clear x Y Yfit
% end
% 
% for ileaf = 1:8
%     corrfact = coeffs3D(1,ileaf)*kpos_phyllo(ileaf,:,1) + coeffs3D(2,ileaf)*kpos_phyllo(ileaf,:,2) +coeffs3D(3,ileaf)*kpos_phyllo(ileaf,:,3);
%     corrfact = permute(corrfact,[3 1 2]);
%     corrfact = repmat(corrfact,[272 1 1]);
%     Scorr(:,ileaf,:) = kdata_phyllo(:,ileaf,:,coil).*exp(-i*corrfact);
%     clear corrfact
% end
% % 
% k = reshape(k_phyllo,[nx*ntviews*nz 3]);
% dcf = reshape(dcf_phyllo,[nx*ntviews*nz 1]);
% rawdata = kdata_phyllo; rawdata = reshape(kdata_phyllo,[nx*ntviews*nz nc]);
% precision = 1E-2;
% siz = [180 180 180];
% for coil = 1:nc
%     k2im(:,:,:,coil) = nufft3_type1(double(k), double(rawdata(:,coil).*dcf), siz, +1,precision);
% end
% images = sum(k2im .* conj(csm),4)./coil_rss;











% G = gpuNUFFT(double(k)',[],2,7,10,[180 180 180],[],true); 
% imags = G'*(rawdata.*repmat(dcf,[1 nc])); imags_phyllo = sqrt(sum(imags.*conj(imags),4)); imagine(imags_phyllo)
% 
% imags = G'*(rawdata(:,:,1).*dcf); 
% imagsScorr = G'*(Scorr(:).*dcf); 
% imagine(cat(2,imags,imagsScorr))




% %% sort the trajectory based on the smallest euclidian distance
% kpos_koosh = squeeze(k_koosh(1,:,:,:)/(.5)); kpos_koosh = reshape(kpos_koosh,[size(kpos_koosh,1)*size(kpos_koosh,2) 3]);
% kpos_phyllo = squeeze(k_phyllo(1,:,:,:)/(.5)); kpos_phyllo = reshape(kpos_phyllo,[size(kpos_phyllo,1)*size(kpos_phyllo,2) 3]);
% 
% kpos = kpos_koosh; 
% 
% A = kpos;
% 
% D = pdist(A);
% Z = squareform(D);  %// Get distance matrix
% 
% N = size(A,1);      %// Store the size of the input array for later usage
% Z(1:N+1:end) = Inf; %// Set diagonals as Infinites as we intend to find
%                     %// minimum along each row
% 
% %// Starting point and initialize an array to store the indices according
% %// to the sorted requirements set in the question
% idx = 1;
% out_idx = zeros(N,1);
% out_idx(1) = idx;
% 
% %// Perform an iterative search to look for nearest one starting from point-1
% for k = 2:N
%     start_ind = idx;
%     [~,idx] = min(Z(start_ind,:));
%     Z(:,start_ind) = Inf;
%     out_idx(k) = idx;
% end
% 
% %// Now that you have the list of indices based on the next closest one, 
% %// sort the input array based on those indices and have the desired output 
% out = A(out_idx,:);
%     
% coil=1;
% phase_phyllo_c = squeeze(phase_phyllo(:,:,:,coil));
% figure;scatter3(phi_phyllo(:),theta_phyllo(:),squeeze(phase_phyllo_c(:)))
% phase_koosh_c = squeeze(phase_koosh(:,:,:,coil)); phase_koosh_c = phase_koosh_c(:);
% phase_koosh_c = phase_koosh_c(out_idx);
% 
% x = 1:length(phase_koosh_c);
% figure;scatter(x(:),phase_koosh_c(:));hold on
% p = polyfit(x(:),phase_koosh_c(:),7);
% corrfact = polyval(p,x(:));
% plot(x(:),corrfact,'--')
% 
% 
% % [nx,ntviews,nz,nc] = size(kdata_koosh);
% % k = reshape(k_koosh,[nx*ntviews*nz 3]);
% % dcf = reshape(dcf_koosh,[nx*ntviews*nz 1]);
% % rawdata = kdata_koosh(:,:,:,coil,1); rawdata = reshape(rawdata,[nx*ntviews*nz 1]);
% % G = gpuNUFFT(double(k)',[],2,7,10,[220 220 220],[],true); 
% % imags = G'*(rawdata.*dcf); imagine(imags)
% % 
% % Scorr = squeeze(kdata_koosh(:,:,:,coil,1)); Scorr = reshape(Scorr,[nx ntviews*nz]).*exp(-i*repmat(corrfact',[nx 1]));
% % imagsS = G'*(Scorr(:).*dcf); imagine(imagsS)
% 
