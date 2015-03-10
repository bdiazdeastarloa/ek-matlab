%% Get data from EK files
clc
clear all
load DATA_EK

% First index is importing country (column 1)
% Second is exporting country (column 2)

% 19 countries

n = 19;

% Data from 1990:

n1 = 6860;
nn = 7220;

%% Trade shares from data on exports and domestic sales

xnimat = datafe(:,25:43);
xrate = datafe(:,2);                % Exchange rates
xratemat = xrate*ones(1,n);
xnimat = 1000000*xnimat;
xnimat = xnimat./xratemat;
xn = sum(xnimat,2);                 % Total spending on manufactures
xnmat = xn*ones(1,n);
Pi_in = xnimat./xnmat;

% Note: on the diagonal of Pi we see domestic sales as a fraction of total
% spending, Xnn/Xn

xnn = diag(xnimat);                   % Domestic spending on manufactures
exps = sum(xnimat)' - xnn;            % Exports 
imps = sum(xnimat,2) - xnn;           % Imports

%% Geographic barriers matrix  

D = datafe(:,6:24);
D = exp(D);                 % Data is ln(dni^(-theta))

%% Labor and wages

l = datafe(:,3);
l = 2080*l;                 % Labor in hours
aw = datafe(:,4);           % Absolute wage
aw = aw/2080;               % Adjust for hours

% Adjust labor and wages for skills

% Human capital
trade2 = trade2(n1:nn,:);   % Data for 1990
hk = trade2(1:n,6);

% Adjust labor and wages for human capital
l = l.*(exp(0.06*hk));
aw = aw.*(exp(-0.06*hk));

%% Income and shares parameters

y = 1000000*(datafe(:,1)./xrate);   % GDP
yl = aw.*l;                         % Manufacturing labor income
yo = y - yl;                        % Nonmanufacturing income

% Check value of beta
betacheck = yl./(xnn+exps);

% Check alpha
myalpha = (yl + (imps - exps))./y;
weight = y/sum(y);
myalpha = weight'*myalpha;

%% Competitiveness

S = datafe(:,5);                     % Source-country competitiveness
S = S - S(n,1);                      % Relative to US


%% Clear unused data

clear datafe trade1 trade2 trade3 pppdat1 getdist

save DATA