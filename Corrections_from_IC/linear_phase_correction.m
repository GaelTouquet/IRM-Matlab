function Dp = linear_phase_correction(Data, dim, golden,koosh)

if nargin < 4
    koosh = false;
end

if nargin < 3
    golden = 0;
end

dims = size(Data);
if numel(dims) <4
    dims(4) = 1;
end

if koosh
    Dc = ktoi(Data,1);
    D_aux = zeros(size(Data));
    for nc = 1:size(Data,4)
        aux = Dc(:,:,:,nc);
    for int_no = 1:size(Data,3)/2
        for spoke_no = 1:size(Data,2);
            
            S = squeeze(aux(:,spoke_no, 2*(int_no-1)+1));
            Sr = squeeze(aux(:,spoke_no, 2*int_no));
            
            linPhase = spatial_linear_phase(S,1,size(Data,1)-1);
            linPhase_alt = spatial_linear_phase(Sr,1,size(Data,1)-1);
            linPhase_avr = (linPhase + linPhase_alt)/2;
            phaseOff = -size(Data,1)/2 * linPhase_avr;
            aux(:,spoke_no, 2*(int_no-1)+1) = complex_mul_linear_phase(S,size(Data,1),linPhase_avr,phaseOff);
            aux(:,spoke_no, 2*int_no) = complex_mul_linear_phase(Sr,size(Data,1),linPhase_avr,phaseOff);
            
        end
    end
    D_aux(:,:,:,nc) = aux;
    end
    
else % 2D
if golden
    if golden == 1
        angles = (0:dims(2)-1)'*(pi/180)*(180*0.618034);
    elseif golden == 2
        angles = (0:dims(2)-1)'*(pi/180)*(180*0.1312674636);
    end
    
    %     pairs = zeros(floor(dims(2)/2), 2);
    %     jj = 1;
    %
    %     for ii = 1:floor(dims(2)/2)
    %
    %         while isnan(angles(jj))
    %             jj = jj+1;
    %         end
    %
    %
    %             pairs(ii,1) = jj;
    %             cur_angle = angles(jj);
    %
    %             angles(jj) = NaN;
    %             angle_diff = mod(angles - cur_angle + pi,2*pi);
    %             [~, I] = min(angle_diff);
    %             angles(I) = NaN;
    %             pairs(ii,2) = I;
    %
    %             jj = jj+1;
    %
    %     end
    %
    %     close all
    %     figure
    %             cm = lines;
    %
    %     angles = (0:dims(2)-1)'*(pi/180)*(180*0.618034);
    %     angle_diffs = abs(mod(angles(pairs(:,1)) - angles(pairs(:,2)) + pi,2*pi)-2*pi)*180/pi
    %
    %         for ii = 1:floor(dims(2)/2)
    % hold off
    %             plot([0 1] * cos(angles(pairs(ii,1))), [0 1] * sin(angles(pairs(ii,1))),'Color',cm(ii,:));
    %            hold on
    %             plot([0 1] * cos(angles(pairs(ii,2))), [0 1] * sin(angles(pairs(ii,2))),'Color',cm(ii,:));
    %     xlim([-1 1])
    %     ylim([-1 1])
    %             waitforbuttonpress
    %         end
    
    
    [angles1, angles2] = ndgrid(angles);
    angle_diff = angles1 - angles2 +pi;
    angle_diff = mod(angle_diff,2*pi);
    %    angle_diff = angle_diff + triu(inf(size(angle_diff)));
    [~, I] = min(angle_diff);
    
end
% 
% if dims(4) > 1 && ~(matlabpool('size'))
%     try
%         matlabpool
%     end
% end


if dim == 2
    reconpars.nz = dims(3);
elseif dim ==3
    % Data = reshape(Data, [dims(1), dims(2)*dims(3), dims(4)]);
    Data = reshape(Data, [dims(1), dims(2)*dims(3), 1, dims(4)]);
    reconpars.nz = 1;
else
    error('dim must be 2 or 3')
end


reconpars.nc = dims(4);
reconpars.projections = dims(2);
reconpars.nr = dims(1)/2;

% Dc = ktoi(squeeze(Data),1);
Dc = ktoi(Data,1);
D_aux = zeros(size(Data));

for nz = 1:reconpars.nz;
    
    for nc = 1:reconpars.nc,
        
        aux = Dc(:,:,nz,nc);
        
        if ~golden
            for j = 1:reconpars.projections/2,
                
                S = squeeze(aux(:,2*(j-1)+1));
                Sr = squeeze(aux(:,2*j));
                
                linPhase = spatial_linear_phase(S,1,2*reconpars.nr-1);
                linPhase_alt = spatial_linear_phase(Sr,1,2*reconpars.nr-1);
                linPhase_avr = (linPhase + linPhase_alt)/2;
                phaseOff = -(reconpars.nr) * linPhase_avr;
                aux(:,2*(j-1)+1) = complex_mul_linear_phase(S,2*reconpars.nr,linPhase_avr,phaseOff);
                aux(:,2*j) = complex_mul_linear_phase(Sr,2*reconpars.nr,linPhase_avr,phaseOff);
                
            end
        else % golden
            iscorrected = false(reconpars.projections);
            linPhase_list = zeros(1,reconpars.projections);
            for j = 1:reconpars.projections,
                
                if ~ iscorrected(j)
                    S = squeeze(aux(:,j));
                    Sr = squeeze(aux(:,I(j)));
                    
                    linPhase = spatial_linear_phase(S,1,2*reconpars.nr-1);
                    linPhase_alt = spatial_linear_phase(Sr,1,2*reconpars.nr-1);
                    linPhase_avr = (linPhase + linPhase_alt)/2;
                    linPhase_list(j) = linPhase_avr; % for debug only
                    phaseOff = -(reconpars.nr) * linPhase_avr;
                    aux(:,j) = complex_mul_linear_phase(S,2*reconpars.nr,linPhase_avr,phaseOff);
                    aux(:,I(j)) = complex_mul_linear_phase(Sr,2*reconpars.nr,linPhase_avr,phaseOff);
                    
                    iscorrected(j) = true;
                    iscorrected(I(j)) = true;
                end
            end
            
            
        end
        
        D_aux(:,:,nz,nc) = aux;
        
    end
end

end

Dp = itok(D_aux,1);


