function [ shift ] = compute_gradient_delay_shift( spoke1,spoke2 )
%COMPUTE_GRADIENT_DELAY_SHIFT Summary of this function goes here
%   Detailed explanation goes here

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
    if imin>1
        imin = imin-1;
    end
    if imax<size(af1,1)
        imax = imax+1;
    end
    downval = af1(imin);
    upval = af1(imax);
end

phases = angle(g(imin:imax));

dir = phases\reshape((imin:imax),size(phases));

shift = dir * 0.5 / (272 * 2 * pi);
end

