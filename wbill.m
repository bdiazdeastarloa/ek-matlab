function ylval = wbill(Pi)

global b n a Y

H = eye(n)-(1-b)*Pi;

ylval = a*b*(inv(H))*Pi*Y;
