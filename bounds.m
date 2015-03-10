%% Calculate price bounds implied by the model

function pbounds = bounds(w)

global b theta t d n g 

temp1 = g*t.*(w.^(-theta*b));
temp2 = repmat(temp1,1,n)';
G = d.*temp2;
Gmax = max(sum(G,2));
pmax = repmat(Gmax^(1/b),n,1);
pmin = (diag(G)).^(1/b);

pbounds = [pmax,pmin];
