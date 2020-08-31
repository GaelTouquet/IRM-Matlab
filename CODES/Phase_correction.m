function [ kdata_corrected ] = Phase_correction( kdata_phyllo )
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
x1 = rem(reshape(phi_phyllo,[nspokes*ninterleaf 1]),pi*2);
x2 = reshape(theta_phyllo,[nspokes*ninterleaf 1]);
x = [x1(:),x2(:)];

for coil = 1:ncoil
    Y(:,coil) = reshape(squeeze(phase_phyllo(:,:,:,coil)),[nspokes*ninterleaf 1]);
    [coeffs3D(:,coil),R,J,CovB,MSE] = nlinfit(x,Y(:,coil),errmodel3D,[0,0,0,0],opts);
    Yfit(:,coil) = errmodel3D(coeffs3D(:,coil),x);
    corrfact = coeffs3D(1,coil)*kpos_phyllo(:,:,1) + ...
        coeffs3D(2,coil)*kpos_phyllo(:,:,2) + ...
        coeffs3D(3,coil)*kpos_phyllo(:,:,3) + ...
        coeffs3D(4,coil);
    corrfact = permute(corrfact,[3 1 2]);
    corrfact = repmat(corrfact,[272 1 1]);
    kdata_corrected(:,:,:,coil) = kdata_phyllo(:,:,:,coil).*exp(-i*corrfact);
end

end

