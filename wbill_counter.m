function ylval = wbill_counter(Pi)

global b n a Yo

H = eye(n)-(1-b+a*b)*Pi;
ylval = a*b*(inv(H))*Pi*Yo;