function [c,ceq]=constr_counter(eqp)

global a b theta t d n L Yo g

p = eqp(1:n,1);
w = eqp(n+1:2*n,1);

%Prices
aux1 = g*t.*w.^(-theta*b);
aux2 = repmat(aux1,1,n)';
G = d.*aux2;
q = p - G*(p.^(1-b)); 


% Shares
pvec = repmat(p,1,n)';
aux3 = aux1.*(p.^(1-b));
aux4 = repmat(aux3,1,n);
Pi = d'.*aux4.*(pvec.^(-1));                      

% Wages adjustment
% Now using non-manufacturing income Yo
H = eye(n)-(1-b+a*b)*Pi;
h = w.*L - a*b*(inv(H))*Pi*Yo;

f = [q;h];

e = 10e-12;
ceq = []; 
%ceq = f;
c = [f-e*ones(2*n,1);-e*ones(2*n,1)-f];
%c = [];
