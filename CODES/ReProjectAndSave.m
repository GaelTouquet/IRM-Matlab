function ReProjectAndSave(kspace_points, the_path, the_name, nc, k2im,csm,coil_rss)
    %ReProjectAndSave - Re-projects an image onto the k-space points, then nufft back into an image, finally writes the image
    %
    % Syntax: ReProjectAndSave(kspace_points)
        %% recompute dcf
        number_of_elements = numel(kspace_points);
        kspace_points = reshape(kspace_points,[number_of_elements/3 3]);
        kspace_points4dcf = permute(kspace_points,[3 2 1]);
        verbose = 1; osf = 1.5; effMtx = 100; numIter = 10;
        predcf = sdc3_MAT(double(kspace_points4dcf),numIter,effMtx,verbose,osf);
        osf = 2.1; numIter = 30;
        dcf = sdc3_MAT(double(kspace_points4dcf),numIter,effMtx,verbose,osf,predcf);
        precision = 1E-2;
        siz = [180 180 180];
        %% Project back into k-space with inverse nufft, with own k distribution
        for coil = 1:nc
            reprojected_data(:,coil) = nufft3_type2(double(kspace_points), double(k2im(:,:,:,coil)),-1,precision);
            re_images(:,:,:,coil) = nufft3_type1(double(kspace_points), double(reprojected_data(:,coil).*dcf), siz, +1,precision);
        end
        re_decoiled_images = sum(re_images .* conj(csm),4)./coil_rss;
        S.(the_name) = re_decoiled_images;
        save(the_path,'-struct','S','-append')
end

