%%%%
function Da = abs_phase_correction(Dp, reconpars,Phase0)

Dc = ktoi(Dp,1); 
Da = zeros(size(Dc));

for nc = 1:reconpars.nc,
    aux = squeeze(Dc(:,:,nc));
   
    for j = 1:reconpars.projections,
        S = aux(:,j);
        absPhase = spatial_absolute_phase(S,1,2*reconpars.nr-1);

         Da(:,j,nc) = complex_mul_absolute_phase(S,2*reconpars.nr, -absPhase + Phase0(:,nc));
    end
end
Da = itok(Da,1);
%%%%%