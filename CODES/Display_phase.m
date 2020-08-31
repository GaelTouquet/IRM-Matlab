function Display_phase( kdata_phyllo, coil )
%DISPLAY Displays the central phase of a phyllotaxis distribution.
phase_phyllo = angle(kdata_phyllo(size(kdata_phyllo,1)/2,:,:,:));

[phi_phyllo, theta_phyllo, kpos_phyllo] = compute_phyllo_angles(kdata_phyllo,0);

nreadout = size(kdata_phyllo,1);
nspokes = size(kdata_phyllo,2);
ninterleaf = size(kdata_phyllo,3);
ncoil = size(kdata_phyllo,4);

figure;
x1 = rem(reshape(phi_phyllo,[nspokes*ninterleaf,1]),pi*2);
x2 = reshape(theta_phyllo,[nspokes*ninterleaf,1]);
x = [x1(:),x2(:)];

if coil == 0
    for icoil = 1:ncoil
        Y(:,icoil) = reshape(squeeze(phase_phyllo(:,:,:,icoil)),[nspokes*ninterleaf 1]);
        subplot(3,9,icoil); scatter3(x(:,1),x(:,2),Y(:,icoil),1,Y(:,icoil));
        xlabel('phi (rad)');ylabel('theta (rad)');zlabel('phase (rad)'); title(['Phyllo - Coil #' num2str(coil)]); % caxis([-3.1416/4 3.1415/4]);
    end
else
    Y(:,coil) = reshape(squeeze(phase_phyllo(:,:,:,coil)),[nspokes*ninterleaf 1]);
    scatter3(x(:,1),x(:,2),Y(:,coil),1,Y(:,coil));
    xlabel('phi (rad)');ylabel('theta (rad)');zlabel('phase (rad)'); title(['Phyllo - Coil #' num2str(coil)]); % caxis([-3.1416/4 3.1415/4]);
end

