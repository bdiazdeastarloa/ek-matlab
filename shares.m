function pival = shares(p,w)

global g t theta b d n

aux1 = g*(t.*(w.^(-theta*b)).*(p.^(1-b)));
aux2 = repmat(aux1,1,n)';
pinv = p.^(-1);
pinvvec = repmat(pinv,1,n);

% Transpose of the trade shares matrix
pival = (d.*aux2.*pinvvec)';  
              
