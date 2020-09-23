function [ kdata_corr ] = gradient_delay_compensation_HU( kdata )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[nreadouts,nspokes,ninterleaves,ncoils] = size(kdata);

kdata_corr = zeros(size(kdata));

for icoil = 1:ncoils
   for iinterleaf = 1:ninterleaves
       for ispokes = 1:nspokes
           spoke1 = kdata(:,ispokes,iinterleaf,icoil);
           if ispokes == nspokes
               spoke2 = kdata(:,ispokes-1,iinterleaf,icoil);
           else
               spoke2 = kdata(:,ispokes+1,iinterleaf,icoil);
           end
           shift = compute_gradient_delay_shift( spoke1,spoke2 );
           kdata_corr(:,ispokes,iinterleaf,icoil) = fft(ifft(spoke1).*exp(-1*1i*shift));
       end
   end
end

end

