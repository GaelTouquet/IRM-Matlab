function DataCorrected  = absolute_phase_correction(Data, dim)

% warning('Changed this file to make 2D work but have not tested 3D since change (UTE paper used previous version)')

dims = size(Data);
if numel(dims) < 4
    dims(4) = 1;
end
% 
% if dims(4) > 1 && ~(matlabpool('size'))
%     try
%         matlabpool
%     end
% end

if dim ==2
    reconpars.nz = size(Data,3);
elseif dim ==3
    % Data = reshape(Data, [dims(1), dims(2)*dims(3), dims(4)]); 
    Data = reshape(Data, [dims(1), dims(2)*dims(3), 1, dims(4)]); 
    reconpars.nz = 1;
else
    error('dim must be 2 or 3')
end

reconpars.nr = size(Data,1)/2;
reconpars.projections = size(Data,2);
% reconpars.nc = size(Data,3);
reconpars.nc = size(Data,4);

DataCorrected = zeros(size(Data ));

reconpars.nz;

for k = 1:reconpars.nz,
    
    if dims(4) > 1;
%         [absPhase] = unwrap(abs_phase_compute(squeeze(Data(:,:,:,k)),reconpars));
        [absPhase] = unwrap(abs_phase_compute(squeeze(Data(:,:,k,:)),reconpars));
        Phase0 = mean(squeeze(absPhase),1);
    else
        Phase0 = 0;
    end
%     DataCorrected(:,:,k,:) = abs_phase_correction(squeeze(Data(:,:,:,k)),reconpars,Phase0);
    DataCorrected(:,:,k,:) = abs_phase_correction(squeeze(Data(:,:,k,:)),reconpars,Phase0);
end



DataCorrected = reshape(DataCorrected,dims);



