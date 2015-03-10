%%
clc
diary('ek_counter.out');
%% ECON 597: Empirical Methods
% 
% Empirical Methods Project - Spring 2011
% Revisiting computational issues of Eaton and Kortum (2002)
% 
% By Bernardo Diaz de Astarloa

% Counterfactuals analysis:

% - Increase in US level of technology.
% - Decrease in trade costs.


%% Arrange data

clear all
run getdata

%% Parameters 

clear all
clc

% Global variables are only parameters which do not change
% Now use non-manufacturing income Yo
global a b theta t d n L Yo g 

load DATA
load baseline

n = 19;                              % Number of countries
a = myalpha;                
b  = 0.21221;
theta = 8.28;                        % Shape of Pareto dist
relaw = aw./aw(n,1);                 % Relative wage to US
L = l;                               % Labor
Yo = y - basew(:,1).*l;              % Non-man income using baseline wages
elast = 0.1;                         % Elasticity of wage adjustment

% Gamma constant (as in the original code)
g = (b^(-b))*((1-b)^(-(1-b)));                    
g = g^(theta);                     

%% Counterfactuals

disp('Counterfactual analysis')

% Base and counterfactual technology
lambdat = 1;                         % Baseline
%lambdat = 1.2;                       % Counterfactual

t = exp( b*(S+theta*(log(relaw))) );  % Technology: see Table VI
t(n) = lambdat*t(n);                  % Change US level

disp('Technology: lambdat=1 is baseline')
disp(lambdat)

% Base and counterfactual barriers
%lambdag = 1;                          % Baseline
lambdag = 0.69;                       % Counterfactual

D = D.^(-1/theta)-ones(n,n);         % dni - 1
D = lambdag*D;                        
d = (D+ones(n,n)).^(-theta);         % Geographic barriers dni^(-theta)
           
disp('Geographic barriers: lambdag=1 is baseline')
disp(lambdag)

%%  Original EK (2002) approach:
%   Solve for the equilibrium using an iterative procedure with a
%   Newton-Raphson intermediate step.
disp('Newton iterative procedure starts here')

tic
w = aw;                     % Set initial wages at data
mytolw = 1e-10;              
myiter = 200;
tolw = 1;
i = 1;

while (tolw>mytolw && i<myiter)
    
    i
    
    % Set bounds (implied by the model) for prices
    pbounds = bounds(w);                   
    p0 = 0.5*(pbounds(:,1)+pbounds(:,2));  % Half way between min and max

%    % Alternative initial price
%    p0 = 2*ones(n,1);
 
    % Given wages, solve for prices using Newton's method
    optset('newton','tol',1e-12);
    optset('newton','showiters',1);
    optset('newton','maxit',200);    
    
    [pval,fval] = newton(@(p) prices(p,w),p0);

    p = pval;                       % Prices
    
    Pi = shares(p,w);               % Get (transpose of) trade shares
    
    yl = wbill_counter(Pi);                 % Get wage bill
    
    exl = (yl./w - l)./l;           % Implied excess labor demand
   
    % Update wages using excess demand for labor
    wnew = w.*(ones(n,1) + exl*elast);      
    w = wnew;
    tolw = max(abs(exl)); 
    i = i+1;
    
    [p,w,exl]
    
end

disp('--- Exit flag ---')
if i == myiter;
    disp('Maximum number of iterations reached');
else
    disp('Loop converged');
end

disp('----------------------')
disp('Time consumed: ')
toc
disp('----------------------')


eqw = w;                                % Equilibrium wages
eqrw = eqw./eqw(n,1);                   % Relative to the US

pbounds = bounds(eqw);
p0 = 0.5*(pbounds(:,1)+pbounds(:,2));
p = newton(@(p) prices(p,eqw),p0);

eqp = p.^(-1/theta);                    % Manufacturing prices
eqrp = eqp./eqp(n,1);                   % Relative to the US

realw = eqw./(eqp.^a);                  % Real wages

tradeshares = shares(p,eqw)';           % Predicted trade shares

neww = [eqw,realw];
newp = [eqp,eqrp];

% Welfare
newwelfare = (yl+Yo).*(eqp.^(-a));

% Change in welfare
gain = 100*(newwelfare - welfare)./welfare;
rgain = 100*gain/gain(n,1);
lgain = 100*log(newwelfare./welfare);

disp('Gain as a % of US gain')
disp(rgain)

disp('% log difference in welfare')
disp(lgain)


%%  KNITRO as alternative algorithm:
%   Use KNITRO package to solve as a constrained optimization problem.
%   Constant objective function and constraints are just the equations
%   we want to solve. 

%   Constraints are stated in constr_counter.m

disp('KNITRO starts here')

tic
w0 = aw;                                       % Set initial wages at data
wmax = inf;
wmin = 0.0000000001*ones(n,1);

% Set bounds (implied by the model) for prices
pbounds = bounds(w0);                   

pmax = pbounds(:,1);
%pmax = inf;
%pmin = pbounds(:,2);
pminalt = 0.00000001*ones(n,1);
pmin = pminalt;
%pmin = -inf;                   

%p0 = 0.5*(pbounds(:,1)+pbounds(:,2));  % Half way between max and min
p0 = 0.5*(pbounds(:,1)+pminalt);
%p0 = 10*ones(n,1);

% Now we solve for prices and wages jointly, so we need an initial vector
eqp0 = [p0;w0];
eqpmax = [pmax;wmax];
eqpmin = [pmin;wmin];

% Fmincon options
% options=optimset('Algorithm','interior-point','MaxIter',2000,'MaxFunEval'...
%      ,10000,'Tolcon',1E-15,'Display','iter', 'TolX',1E-15);

% User supplied KNITRO options in 'knitropts.opt'
options=optimset('Tolcon',10E-15, 'TolX',10E-15);
  [pval,fval,exitflag,output,lambda] =...
     ktrlink(@obj,eqp0,[],[],[],[],eqpmin,eqpmax,...
     @(eqp) constr_counter(eqp),options,'knitropts.opt');

disp('----------------------')
disp('Time consumed: ')
toc
disp('----------------------')

eqw = pval(n+1:2*n,1);                  % Equilibrium wages
eqrw = eqw./eqw(n,1);                   % Relative to the US

p = pval(1:n,1);                        % Equilibrium prices^(-theta)

eqp = p.^(-1/theta);                    % Manufacturing prices
eqrp = eqp./eqp(n,1);                   % Relative to the US

realw = eqw./(p.^a);                    % Real wages

tradeshares = shares(p,eqw)';           % Predicted trade shares

nitroneww = [eqw,realw,eqrw];
nitronewp = [eqp,eqrp];

yl = wbill_counter(tradeshares');

% Welfare
nitronewwelfare = (yl+Yo).*(eqp.^(-a));

% Change in welfare
gain = 100*(newwelfare - welfare)./welfare;
rgain = 100*gain/gain(n,1);
lgain = 100*log(newwelfare./welfare);

disp('Gain as a % of US gain')
disp(rgain)

disp('% log difference in welfare')
disp(lgain)







diary('off');
