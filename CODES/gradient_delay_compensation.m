function [ kdata_corr ] = gradient_delay_compensation( kdata )
%GRADIENT_DELAY_COMPENSATION Corrects radial data for gradient delays
%   from https://cds.ismrm.org/protected/11MProceedings/files/2816.pdf

[nreadouts,nspokes,ninterleaves,ncoils] = size(kdata);

for i_coil = 1:ncoils

    spoke1 = kdata(:,7,2483,i_coil);
    spoke2 = kdata(:,6,2448,i_coil);

    xshift = compute_gradient_delay(spoke1,spoke2);


    spoke1 = kdata(:,7,2375,i_coil);
    spoke2 = kdata(:,6,2573,i_coil);

    yshift = compute_gradient_delay(spoke1,spoke2);


    spoke1 = kdata(:,1,2508,i_coil);
    spoke2 = kdata(:,2,158,i_coil);

    zshift = compute_gradient_delay(spoke1,spoke2);

end

end

