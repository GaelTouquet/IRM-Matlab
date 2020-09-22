function [ kdata_corrected ] = Phase_correction_linear( kdata_phyllo )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nreadout = size(kdata_phyllo,1);
nspokes = size(kdata_phyllo,2);
ninterleaf = size(kdata_phyllo,3);
ncoil = size(kdata_phyllo,4)

for coil = 1:ncoil
   for interleaf = 1:ninterleaf
       for spoke = 1:nspokes
           %derive the linear factor
           linear_factor = sum(kdata_phyllo()
       end
   end
end

end

