function linPhase = spatial_linear_phase( S, nFirst, nLast );
% linPhase = spatial_linear_phase( S, nFirst, nLast )
%
% uses S = sum_{k} S_{k} * conj(S_{k+1}) to estimate the linear phase
%
% from: spatialpc_ts.c

e = sum(S(nFirst:nLast-1).* conj(S(nFirst+1:nLast)));
mag = abs(e);

%ei = 0;
%er = 0;
%Sreal = real(S);
%Simag = imag(S);
%for l=nFirst:nLast,
%  er = er + Sreal(l)*Sreal(l+1) + Simag(l)*Simag(l+1);
%  ei = ei + Simag(l)*Sreal(l+1) - Sreal(l)*Simag(l+1);
%end
%mag = abs(complex(er,ei));

if( mag ~= 0 )
  %ei = ei/mag;
  %er = er/mag;
  %linPhase = atan2(ei, er);
  e = e/mag;
  linPhase = atan2(imag(e), real(e));
else
  linPhase = 0;
end

