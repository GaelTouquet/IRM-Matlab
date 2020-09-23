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
%% try with real data
% Fs = 272;
% T = 1/Fs;
% L = 272;

spoke1 = kdata(:,7,2483,1);
spoke2 = kdata(:,6,2448,1);

f1 = fft(abs(spoke1));
f2 = fft(flip(abs(spoke2)));

% f = Fs*(0:(L/2))/L;

g = f1.*conj(f2);
af1 = abs(fft(spoke1));
[M,I] = max(af1);

upval = M;
downval = M;
imin = I;
imax = I;

while downval>M*0.1 && upval>M*0.1
    imin = imin-1;
    imax = imax+1;
    downval = af1(imin);
    upval = af1(imax);
end

phases = angle(g(imin:imax));

dir = phases\reshape((imin:imax),size(phases));

linshift = dir * 0.5 / (272 * 2 * pi);

% gp2 = abs(g/L);
% gp1 = gp2(1:L/2+1);
% gp1(2:end-1) = 2*gp1(2:end-1);
% plot(f,gp1)
% 
% fp2 = abs(f1/L);
% plot(fp2)
% 
% cros = fft(g);
% 
% plot(abs(cros))


%% 
Fs = 100;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 10000;             % Length of signal
t = (0:L-1)*T;        % Time vector


bla = cos(2*pi*t);
blu = cos(2*pi*4.9*t);

fbla = fft(bla);
fblu = fft(blu);

blap2 = abs(fbla/L);
blap1 = blap2(1:L/2+1);
blap1(2:end-1) = 2*blap1(2:end-1);
blup2 = abs(fblu/L);
blup1 = blup2(1:L/2+1);
blup1(2:end-1) = 2*blup1(2:end-1);

f = Fs*(0:(L/2))/L;

g = fbla.*conj(fblu);
gp2 = abs(g/L);
gp1 = gp2(1:L/2+1);
gp1(2:end-1) = 2*gp1(2:end-1);
plot(f,gp1)
% 
% glu = fbla.*conj(fblu);
% cross = fft(abs(glu));

