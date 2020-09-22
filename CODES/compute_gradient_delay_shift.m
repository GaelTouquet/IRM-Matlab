function [ shift ] = compute_gradient_delay_shift( spoke1,spoke2 )
%COMPUTE_GRADIENT_DELAY_SHIFT Summary of this function goes here
%   Detailed explanation goes here
S1 = spoke1;
S2 = flip(spoke2);
F1 = fft(S1);
F2 = fft(S2);
g = F1.*conj(F2);
phase = angle(g);
test1 = fft(g);
test2 = ifft(g);
k = diff(phase);

end

