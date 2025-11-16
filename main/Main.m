% B-spline based growth transition matrix for capelin in Barents sea
% ** Estimated run time of this code is less than 1 minutes **

clear; clc;
rng('default');

%% Initialization
l_min = 5;
l_max = 21;

year_start = 1988;                     
year_finish = 2017;

years = year_start:year_finish-1;

year = year_start;                    % Starting point of initial data
degree = 2;                           % Degree of the B-spline polynomials
n_knote = 31;                         % Number of B-spline basis functions
Lf = 14;                              % Fecundity length
dt = 1;                               % Discritization of the time
t = 1:dt:(year_finish-year_start)+1;    % Time period of simulation 
step = 0.5;                           % Discritization of the length

l_mid_min = l_min+step/2;
l_mid_max = l_max-step/2;

x = l_min:step:l_max;
x_mid = l_mid_min:step:l_mid_max;
nx = length(x_mid);
nt = length(t);
x_knotes = linspace(min(x),max(x),n_knote);

global model
model.nknote = n_knote;
model.knot = knots(x_knotes,n_knote,degree);
model.degree = degree;
model.x_mid = x_mid;
model.x_edges= x;
model.Lf = Lf;
model.nx = nx;
model.Lf = Lf;
model.year_start = year_start;

N3=zeros(nx,nt);
M2=zeros(3,nt-1);
P2=zeros(2,nt-1);
GT2=zeros(32,32,nt-1); 

Fval = zeros(1,nt-1);
Flag = zeros(1,nt-1);

%X0 = zeros(2,nt-1);
load('initial_point.mat')

I1 = Bspline_int1(x_mid,step);
I2 = Bspline_int2(x_mid,step);
model.I1 = I1;
model.I2 = I2;

%% Parameters
%fa = exp(-0.2049);
%fb = 3.4871;
%f = fecundity(x_mid,fa,fb);    % fecundity
%v = valnerability(x_mid);      % valnerability
%matu = maturity(x_mid);        % maturity

%% Dynamics

% Data
[a1, a2, a3, a4] = Initial_Data(year);
N3(:,1) = a3; 

for j = 2:nt
    [a1, a2, a3, a4] = Initial_Data(year+j-2);

    % Fishing from Jan to April of year j
    [A1,A2,A3,A4,C1,C2,C3,C4] = Initial_Data_new(year+j-1);
       
    % Stock assessment in the 1 of October of year j:  
    % Dynamic of the age class 3

    [G, p1,p2,m1,m2,m3,fval,exitflag] = GTMfun(a2,A3,C2);

    Survived = a2.*(1-maturity(p1,p2));
    Mortality = mortality(m1,m2,m3);
    N3(:,j) = G*(dt*exp(-12*Mortality).*Survived-exp(-6*Mortality).*C2);

    M2(:,j-1) = [m1,m2,m3];
    P2(:,j-1) = [p1,p2];
    GT2(:,:,j-1) = G;

    Fval(j-1) = fval;
    Flag(j-1) = exitflag;
       
end

%  Write_result(B,C,VB);                    % Documenting the results
%save('initial_point.mat', 'X0');   % to write initial parameter valus
%% Figures

%for i=2:nt
   %testGPT(year,year+i,x_mid,N3)           % Comparing empirical...
%end                                               % data with simulations
idx = k_mean_GTM2(GT2,years);
[N_l, N_h, A3] = monte_carlo_compare(idx, GT2);
%[N_l,N_h,A3] = compare(idx, GT2);
%[N_l,N_h,A3] = simulation(idx, GT2);
%MortalityPlot(years,M2)
%SumMortalityPlot(years,M2)
%MovingAverageMortalityPlot1(years, M2, 3)
%MaturityPlot(years,P2)
