function [ kdata_corrected, coeffs3D ] = Phase_correction( kdata_phyllo, spokes, display )
%PHASE_CORRECTION Correction from central phase in phyllotaxis distributions
%   3D generalisation from MOUSSAVI 2014
phase_phyllo = angle(kdata_phyllo(size(kdata_phyllo,1)/2,:,:,:));

[phi_phyllo, theta_phyllo, kpos_phyllo] = compute_phyllo_angles(kdata_phyllo,0);

nreadout = size(kdata_phyllo,1);
nspokes = size(kdata_phyllo,2);
ninterleaf = size(kdata_phyllo,3);
ncoil = size(kdata_phyllo,4);
% TODO if elements of phase array both close to -pi et pi (phase plane goes
% through modulus) then add 2pi to negative elements. Then when correcting,
% make sure that if correction is lower than -pi, add 2pi to correction

errmodel3D = @(coeffs,x) coeffs(1).*sin(x(:,1)).*sin(x(:,2)) + ...
    coeffs(2).*cos(x(:,1)).*sin(x(:,2)) + ...
    coeffs(3).*cos(x(:,2)) + coeffs(4);
opts = statset('nlinfit');
opts.RobustWgtFun = 'andrews';
opts.TolFun = 10^-12;
opts.MaxFunEvals = 10^5;


if spokes == 0
    x1 = rem(reshape(phi_phyllo,[nspokes*ninterleaf 1]),pi*2);
    x2 = reshape(theta_phyllo,[nspokes*ninterleaf 1]);
    x = [x1(:),x2(:)];
    for coil = 1:ncoil
        Y(:,coil) = reshape(squeeze(phase_phyllo(:,:,:,coil)),[nspokes*ninterleaf 1]);

        % check wether the phase plane goes through the +pi plane by checking
        % if there are points that are above 3 and also points below -3
        close_to_top = false;
        close_to_bottom = false;
        for i = 1:nspokes*ninterleaf
            if Y(i,coil) < -3
                close_to_bottom = true;
            end
        end
        for i = 1:nspokes*ninterleaf
            if Y(i,coil) > 3
                close_to_top = true;
            end
        end
        if close_to_top && close_to_bottom
            % put everthing that is negative above zero
            for i = 1:nspokes*ninterleaf
                if Y(i,coil) < 0
                    Y(i,coil) = Y(i,coil) + 2*pi;
                end
            end
        end


        [coeffs3D(:,coil),R,J,CovB,MSE] = nlinfit(x,Y(:,coil),errmodel3D,[0,0,0,0],opts);
        corrfact = coeffs3D(1,coil)*kpos_phyllo(:,:,1) + ...
            coeffs3D(2,coil)*kpos_phyllo(:,:,2) + ...
            coeffs3D(3,coil)*kpos_phyllo(:,:,3) + ...
            coeffs3D(4,coil);
        corrfact = permute(corrfact,[3 1 2]);
        corrfact = repmat(corrfact,[272 1 1]);
        kdata_corrected(:,:,:,coil) = kdata_phyllo(:,:,:,coil).*exp(-1*1i*corrfact);
    end
else
    for coil = 1:ncoil
        for spoke = 1:nspokes
            x1 = rem(reshape(phi_phyllo(spoke,:),[ninterleaf 1]),pi*2);
            x2 = reshape(theta_phyllo(spoke,:),[ninterleaf 1]);
            x = [x1(:),x2(:)];
            Y(:,coil) = reshape(phase_phyllo(:,spoke,:,coil),[ninterleaf 1]);

            % check wether the phase plane goes through the +pi plane by checking
            % if there are points that are above 3 and also points below -3
            close_to_top = false;
            close_to_bottom = false;
            for i = 1:ninterleaf
                if Y(i,coil) < -3
                    close_to_bottom = true;
                end
            end
            for i = 1:ninterleaf
                if Y(i,coil) > 3
                    close_to_top = true;
                end
            end
            if close_to_top && close_to_bottom
                % put everthing that is negative above zero
                for i = 1:ninterleaf
                    if Y(i,coil) < 0
                        Y(i,coil) = Y(i,coil) + 2*pi;
                    end
                end
            end


            [coeffs3D(:,coil,spoke),R,J,CovB,MSE] = nlinfit(x,Y(:,coil),errmodel3D,[0,0,0,0],opts);
            corrfact = coeffs3D(1,coil,spoke)*kpos_phyllo(spoke,:,1) + ...
                coeffs3D(2,coil,spoke)*kpos_phyllo(spoke,:,2) + ...
                coeffs3D(3,coil,spoke)*kpos_phyllo(spoke,:,3) + ...
                coeffs3D(4,coil,spoke); 
            corrfact = permute(corrfact,[3 1 2]);
            corrfact = repmat(corrfact,[272 1 1]);
            kdata_corrected(:,spoke,:,coil) = kdata_phyllo(:,spoke,:,coil).*exp(-1*1i*corrfact);
        end
end

if display
    Display_phase( kdata_phyllo, 0, coeffs3D )
end
    

end

