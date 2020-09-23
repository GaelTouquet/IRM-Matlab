function [ shift ] = Phase_correction_HU( Spoke1, Spoke2 )
%PHASE_CORRECTION_HU Summary of this function goes here
%   Detailed explanation goes here
S1 = abs(Spoke1);
S2 = abs(flip(Spoke2));
F1 = fft(S1);
F2 = fft(S2);
g = F1.*conj(F2);
phase = angle(g);
k = diff(phase);

end

