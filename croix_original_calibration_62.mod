%FAIR WAGE MODEL 

var
c    	    % consumption  
invest      % capital investment
e           % worker effort
q           % ratio mu/lamb or real value of capital
k    	    % capital stock
rk          % rental price capital
R           % nominal gross quarterly interest rate
lamb  	    % Lagrange multiplier household budget constraint
u    	    % capacity utilization
s           % price dispersion
y    	    % output/aggdemand
x1          % term for recursive expression for price 
x2          % term for recursive expression for price
pi   	    % gross quarterly inflation
w   	    % real wage
wtil        % real optimal wage
z    	    % real marginal cost
ptil        % relative optimal price
N	        % aggregate employment rate
f1          % term for recursive expression for optimal wage
f2          % term for recursive expression for optimal wage
g           % government spending
yd          % deviation of y from new steady state               
pi_y        % annual inflation (%)
pi_q        % quarterly inflation (%)  
                        
; 

varexo
pi_tar		% gross quarterly inflation target
R_tar		% nominal interest rate target
phi0            % effort function parameter
phi1            % effort function parameter
;

%	PARAMETERS

parameters 
beta		    % discount rate
h       	    % habit parameter
sigmac          % Consumption CRRA parameter
phi2            % effort function parameter
phi3            % effort function parameter
psi             % effort function parameter
theta           % final firms substitution parameter
alpha           % intermediate production function cobb douglas
Phi             % Fixed Cost of production
xip             % calvo probability price
xiw             % calvo probability wage
gammap          % price index
gammaw          % wage index
tau 		    % capital depreciation rate 
delta           % investment adj cost function
gam1		    % cap utilization function
gam2		    % cap utilization function
rho             % monetary policy parameter
g_ss            %government spending steady state
;

%parameter calibration

% 1) PREFERENCES
h = 0.445;              % degree of habit persistence
beta = 1.03^(-1/4);	    % subjective discount factor
sigmac= 1.755;          %CRRA parameter

%EFFORT
psi = 0.358;                 %substitution between different types of effort incentives
phi2 = 0.004;              %effect of employment on effort                
phi3 = 0.785;              %effect of wages on effort

% 2) ELASTICITIES OF SUBSTITUTION
theta = 5/6;            % price-elasticity of demand for a differianted good

% 3) TECHNOLOGY
tau = 0.025;		% depreciation rate
alpha = 0.76;		% share of labor (it works for 0,76 +-0,02) 
delta = 5.59;       % investment adj cost function from De la Croix
gam1 = 0.032417;	% cap utilization function = 1/beta - (1-tau)
gam2 = 0.000324;	% cap utilization function
Phi = 0.25;  % firm fixed cost 

% 3) CALVO PARAMETERS
xip = 0.892;	 % on prices
xiw = 0.822;      % on wages
gammap = 0.9;  % price index to past inflation ----> Change here for different indexation degree
gammaw = 0.9; % wage index to past inflation ----> Change here for different indexation degree

% 5) MONETARY POLICY
rho = 1.5;		% monetary policy parameter. 
% rho = 3;      % change here for the rho=3 case

% STEADY STATE
g_ss = 0.18;

%-------------------------------------------MODEL--------------------------------------------------------------%

model;

%capital accumulation Eq (6)
k = (1-tau)*k(-1) + invest - invest*(delta/2)*((invest/invest(-1))-1)^2;

%----HOUSEHOLD FIRST ORDER CONDITIONS for consumption, effort and bond----%

%household euler Eq (8) 
lamb= (c-h*c(-1))^(-sigmac)-h*beta*(c(+1)-h*c)^(-sigmac);

% household effort Eq (9) (now that we consider calvo in wages, effort function may be different, considering that xiw families receive w^* and 1-xiw receive w)
e = phi1*((w^psi)*(1-phi3)-phi2*((1/(1-N))^psi)-(phi0-phi2-phi3))/psi;

%household bond demand Eq (10)
lamb = beta*R*lamb(+1)/(pi(+1));

%----HOUSEHOLD FIRST ORDER CONDITIONS for capital, investment and utilization rate----%

%household capital Eq (11)
q*lamb=beta*lamb(+1)*((1-tau)*q(+1) + rk(+1)*u(+1) - gam1*(u(+1)-1) - (gam2/2)*((u(+1)-1)^2));

%household investment Eq (12)
lamb = lamb*q*(1 - (delta/2)*(((invest/invest(-1))-1)^2) - delta*(invest/invest(-1))*((invest/invest(-1))-1)) + beta*lamb(+1)*q(+1)*delta*((invest(+1)/invest)^2)*((invest(+1)/invest)-1);

%household utilization rate Eq (13)
rk=gam1+gam2*(u-1);


%----FIRMS DECISIONS----%

%INTERMEDIATE RETAIL FIRMS-----------------------------------------------------------------------------------
%Optimal price setting 
x1=z*y*(ptil^(1/(theta-1)-1))+beta*xip*(lamb(+1)/lamb)*((ptil/ptil(+1))^(1/(theta-1)-1))*(((pi^gammap)*(pi_tar^(1-gammap))/pi(+1))^(1/(theta-1)))*x1(+1);
x2=y*(ptil^(1/(theta-1)))+beta*xip*(lamb(+1)/lamb)*((ptil/ptil(+1))^(1/(theta-1)))*(((pi^gammap)*(pi_tar^(1-gammap))/pi(+1))^(theta/(theta-1)))*x2(+1);
theta*x2=x1;

%aggregate price law of motion 
1=xip*(pi^(theta/(1-theta)))*((pi(-1)^(gammap))*(pi_tar^(1-gammap)))^(theta/(theta-1))+(1-xip)*(ptil)^(theta/(theta-1));

%price dispersion 
s=(1-xip)*(ptil)^(1/(theta-1))+xip*(pi/(pi(-1)^(gammap)*pi_tar^(1-gammap)))^(1/(1-theta))*s(-1);

%INTERMEDIATE PRODUCERS-----------------------------------------------------------------------------------

%fraction of 1 hour demand
w/e=z*alpha*(((e*N)^(alpha-1))*((u*k(-1))^(1-alpha)));

%capital stock demand (problem)
rk=z*(1-alpha)*((e*N)^(alpha)*((u*k(-1))^(-alpha)));

%optimal wage setting
f1=z*((s*y+Phi)/e)*(wtil^psi)+beta*xiw*((wtil/wtil(+1))^psi)*(lamb(+1)/lamb)*(((pi^gammaw)*(pi_tar^(1-gammaw)))^psi)*f1(+1);
f2=wtil*N + beta*xiw*(wtil/wtil(+1))*(lamb(+1)/lamb)*(((pi^gammaw)*(pi_tar^(1-gammaw)))/pi(+1))*f2(+1);
alpha*phi1*f1=f2;
w=xiw*w(-1)*((pi(-1)^gammaw)*(pi_tar^(1-gammaw))/pi)+(1-xiw)*wtil;

%----AGGREGATE ECONOMY----%

%aggregate output
s*y = ((e*N)^(alpha))*((u*k(-1))^(1-alpha))-Phi;

%aggregate demand 
y=c+invest+g+(gam1*(u-1)+(gam2/2)*(u-1)^2)*k(-1);

%government 
g=g_ss;

%----MONETARY POLICY AND EXOGENOUS VARIABLES----%
%MONETARY POLICY Eq (35)
(R/R_tar)=(pi/pi_tar)^(rho);

%variables generated for the simulations
pi_y = 100*((pi^4)-1);
pi_q = 100*(pi-1);
yd = 100*(y-1.43763)/1.43763;     % steady state output is 1.43763

end;

%-------------------------------------------END OF MODEL--------------------------------------------------------------%


%------------------------------------STEADY STATE VALUES (initial)----------------------------------------------------%

initval;
%STEADY STATE FOR EXOGENOUS VARIABLES 
%pi_tar = 1.04^0.25;		% old inflation target=4%. Change here for the cold turkey exercise
pi_tar = 1.06^0.25;		% inflation target=6%. Change here for the cold turkey exercise
%pi_tar = 1.08^0.25;		    % inflation target=8%. Change here for the cold turkey exercise
R_tar = pi_tar/beta;	    % nominal interest rate target    
phi0= 0.619285104970736;   %effort function scale parameter for N=0.95 and e=1 (comes from finding w_ss and using this in consumer optimal effort) -----> when inflation changes, dont forget to update this
phi1= 0.936460522539359;   %effort function scale parameter for N=0.95 and e=1 (comes from finding w_ss and using this in firm optimal effort)-----> when inflation changes, dont forget to update this

% Nominal gross interest rate---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
R = R_tar;    

% EXACT STEADY STATE VALUES
% gross inflation---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pi = pi_tar;  

% capacity utilization----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
u = 1;              

%real value of capital----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
q=1;

%capital rental price------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rk=gam1;

%optimal relative price from aggregate price movement----------------------------------------------------------------------------------------------------------------------------------------------------------------------
ptil=1;

%real marginal cost from theta*x2=x1---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
z=theta;   

%price dispersion applying above conditions--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
s=1;

%wage in steady state from alpha*phi1*f1=f2--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
w = ((alpha*phi1*z*(1-beta*xiw)/(1-beta*xiw*(pi_tar^psi)))*((((z*(1-alpha))/rk))^((1-alpha)/alpha)))^(1/(1-psi));

%optimal wage and real wage in steady state---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
wtil= w;

%effort from firm CPO-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
e = (w/(alpha*z))*(((z*(1-alpha))/rk)^((alpha-1)/alpha));

%N in steady state from effort equation------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
N = 1 - ((1/phi2)*((1-phi3)*(w^psi) - psi*e/phi1 - phi0+phi2+phi3))^(-1/psi);

%capital in steady state found from the firm CPO---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
k=e*N*((z*(1-alpha)/rk)^(1/alpha));    

%aggregate output found by using k and N in the aggregate production function----------------------------------------------------------------------------------------------------------------------------------------------
y= (1/s)*(((e*N)^alpha)*((u*k)^(1-alpha))-Phi); 

%recursive expression price------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
x1=z*y/(1-beta*xip);          

%recursive expression price------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
x2 = y/(1-beta*xip); 

%recursive expression wage-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
f1=z*N*(((z*(1-alpha))/rk)^((1-alpha)/alpha))*(w^psi)/(1-beta*xiw*(pi_tar^psi));

%recursive expression wage-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
f2=w*N/(1-beta*xiw);

%investment---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
invest=tau*k; 

%government found as resid of y-c-i=g from aggregate demand-----------------------------------------------------------------------------------------------------------------------------------------------------------------
g=g_ss;  

%consumption in steady state------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
c = y - invest - g;

%lagrange multiplier in steady steate-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
lamb=(c*(1-h))^(-sigmac)*(1-h*beta);

%simulation variables
pi_y = 100*((pi^4)-1);
pi_q = 100*(pi-1);
end;

resid;
steady;
check; 


%------------------------------------STEADY STATE VALUES (endval)----------------------------------------------------%
endval;
% STEADY STATE FOR EXOGENOUS VARIABLES 
pi_tar = 1.02^0.25;		    % new inflation target=2%.
R_tar = pi_tar/beta;	    % nominal interest rate target    
phi0= 0.625186643687715;   %effort function scale parameter for N=0.95 and e=1 (comes from finding w_ss and using this in consumer optimal effort) -----> when inflation changes, dont forget to update this
phi1= 0.951143624313913;   %effort function scale parameter for N=0.95 and e=1 (comes from finding w_ss and using this in firm optimal effort)-----> when inflation changes, dont forget to update this

% Nominal gross interest rate---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
R = R_tar;    

% EXACT STEADY STATE VALUES
% gross inflation---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pi = pi_tar;  

% capacity utilization----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
u = 1;              

%real value of capital----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
q=1;

%capital rental price------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rk=gam1;

%optimal relative price from aggregate price movement----------------------------------------------------------------------------------------------------------------------------------------------------------------------
ptil=1;

%real marginal cost from theta*x2=x1---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
z=theta;   

%price dispersion applying above conditions--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
s=1;

%wage in steady state from alpha*phi1*f1=f2--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
w = ((alpha*phi1*z*(1-beta*xiw)/(1-beta*xiw*(pi_tar^psi)))*((((z*(1-alpha))/rk))^((1-alpha)/alpha)))^(1/(1-psi));

%optimal wage and real wage in steady state---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
wtil= w;

%effort from firm CPO-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
e = (w/(alpha*z))*(((z*(1-alpha))/rk)^((alpha-1)/alpha));

%N in steady state from effort equation------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
N = 1 - ((1/phi2)*((1-phi3)*(w^psi) - psi*e/phi1 - phi0+phi2+phi3))^(-1/psi);

%capital in steady state found from the firm CPO---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
k=e*N*((z*(1-alpha)/rk)^(1/alpha));    

%aggregate output found by using k and N in the aggregate production function----------------------------------------------------------------------------------------------------------------------------------------------
y= (1/s)*(((e*N)^alpha)*((u*k)^(1-alpha))-Phi); 

%recursive expression price------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
x1=z*y/(1-beta*xip);          

%recursive expression price------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
x2 = y/(1-beta*xip); 

%recursive expression wage-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
f1=z*N*(((z*(1-alpha))/rk)^((1-alpha)/alpha))*(w^psi)/(1-beta*xiw*(pi_tar^psi));

%recursive expression wage-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
f2=w*N/(1-beta*xiw);

%investment---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
invest=tau*k; 

%government found as resid of y-c-i=g from aggregate demand-----------------------------------------------------------------------------------------------------------------------------------------------------------------
g=g_ss;  

%consumption in steady state------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
c = y - invest - g;

%lagrange multiplier in steady steate-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
lamb=(c*(1-h))^(-sigmac)*(1-h*beta);

%simulation variables
pi_y = 100*((pi^4)-1);
pi_q = 100*(pi-1); 

end;

resid;
steady;

check; 


%------------------------------------PERFECT FORESIGHT SIMULATION-------------------------------------------------%
perfect_foresight_setup(periods=100);
perfect_foresight_solver;
