function [fval,fjac] = prices(p,w)

global b theta t d n g 


aux1 = g*t.*(w.^(-theta*b));
aux2 = repmat(aux1,1,n)';
G = d.*aux2;
fval = inv(G)*p -p.^(1-b);
fjac = inv(G) - (1-b)*diag(p.^(-b));