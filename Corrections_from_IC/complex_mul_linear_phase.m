%%%%%%%%%%%
function Scorr = complex_mul_linear_phase(S,nx,linPhase,phaseOff );
% phase_corrected = complex_mul_linear_phase(S,nx,linPhase,phaseOff );

%Sreal    = real(S);
%Simag    = imag(S);
%for x=0:nx-1,
%phase    = x*linPhase + phaseOff;
%  sinus    = sin(phase);
%  cosinus  = cos(phase);
%  Scorr_real(x+1) = Sreal(x+1)*cosinus - Simag(x+1)*sinus;
%  Scorr_imag(x+1) = Sreal(x+1)*sinus   + Simag(x+1)*cosinus;
%end
%Scorr = complex(Scorr_real,Scorr_imag);

x = 0:nx-1;
phase = x*linPhase + phaseOff;
Scorr = S(:) .* exp(1i*phase(:));

