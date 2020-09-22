function [ kdata_corr ] = gradient_delay_compensation( kdata )
%GRADIENT_DELAY_COMPENSATION Corrects radial data for gradient delays
%   from https://cds.ismrm.org/protected/11MProceedings/files/2816.pdf

% [nreadouts,nspokes,ninterleaves,ncoils,nphasecontrast] = size(kdata_phyllo);

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
% 
% S1 = abs(Spoke1);
% S2 = abs(flip(Spoke2));
% F1 = fft(S1);
% F2 = fft(S2);
% g = F1.*conj(F2);
% phase = angle(g);
% k = diff(phase);
% 
% 
% for i_coil = 1:n_coils:
%     spoke1 = 
% end

end

