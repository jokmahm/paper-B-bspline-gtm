% In this program I test new growth transition matrix

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

Fval = zeros(1,nt-1);
Flag = zeros(1,nt-1);

%X0 = zeros(2,nt-1);
load('initial_point.mat')

I1 = Bspline_int1(x_mid,step);
I2 = Bspline_int2(x_mid,step);
model.I1 = I1;
model.I2 = I2;

% Data
[a1, a2, a3, a4] = Initial_Data(year);
N3(:,1) = a3; 

% ... (everything above unchanged) ...

GT2_bs  = zeros(nx,nx,nt-1);   % B-spline
GT2_nA  = zeros(nx,nx,nt-1);   % Normal-A (Chen)

Nhat_bs = zeros(nx,nt); Nhat_bs(:,1)=N3(:,1);
Nhat_nA = zeros(nx,nt); Nhat_nA(:,1)=N3(:,1);

Err_bs  = zeros(1,nt-1);
Err_nA  = zeros(1,nt-1);

M2_bs = zeros(3,nt-1);  P2_bs = zeros(2,nt-1);
M2_nA = zeros(3,nt-1);  P2_nA = zeros(2,nt-1);

% --- add collectors next to the others ---
GT2_gam = zeros(nx,nx,nt-1);
Nhat_gam = zeros(nx,nt); Nhat_gam(:,1) = N3(:,1);
Err_gam  = zeros(1,nt-1);
M2_gam   = zeros(3,nt-1);  P2_gam = zeros(2,nt-1);

% --- add alongside other collectors ---
GT2_lgn = zeros(nx,nx,nt-1);
Nhat_lgn = zeros(nx,nt); Nhat_lgn(:,1) = N3(:,1);
Err_lgn  = zeros(1,nt-1);
M2_lgn   = zeros(3,nt-1);  P2_lgn = zeros(2,nt-1);

% --- for runtime
t_logn  = zeros(1,nt-1);
t_norm  = zeros(1,nt-1);
t_bs    = zeros(1,nt-1);
t_gamma = zeros(1,nt-1);

for j = 2:nt
    [a1, a2, a3, a4] = Initial_Data(year+j-2);
    [A1,A2,A3,A4,C1,C2,C3,C4] = Initial_Data_new(year+j-1);

    % ---------- Log-Normal (new) ----------
    tic;
    [G_ln, p1_l,p2_l,m1_l,m2_l,m3_l, Linf_l,K_l,CV_l, fval_l,ef_l] = ...
        GTMfun_Lognormal(a2, A3, C2, model.x_mid(:), model.x_edges(:));
    t_logn(j-1) = toc;

    Survived_l = a2.*(1-maturity(p1_l,p2_l));
    Mort_l     = mortality(m1_l,m2_l,m3_l);
    Nhat_lgn(:,j)= G_ln*(dt*exp(-12*Mort_l).*Survived_l - exp(-6*Mort_l).*C2);
    GT2_lgn(:,:,j-1) = G_ln;
    M2_lgn(:,j-1) = [m1_l,m2_l,m3_l];  P2_lgn(:,j-1) = [p1_l,p2_l];
    Err_lgn(j-1)  = norm(Nhat_lgn(:,j) - A3);
    

    % ---------- Normal (Chen-style) ----------
    tic;
    [G_nA, p1_a,p2_a,m1_a,m2_a,m3_a, Linf_a,K_a,sL_a,sK_a,rho_a, fval_a,ef_a] = ...
        GTMfun_Normal(a2, A3, C2, model.x_mid(:), model.x_edges(:));
    t_norm(j-1) = toc;

    Survived_a = a2.*(1-maturity(p1_a,p2_a));
    Mort_a     = mortality(m1_a,m2_a,m3_a);
    Nhat_nA(:,j)= G_nA*(dt*exp(-12*Mort_a).*Survived_a - exp(-6*Mort_a).*C2);
    GT2_nA(:,:,j-1) = G_nA;
    M2_nA(:,j-1) = [m1_a,m2_a,m3_a];  P2_nA(:,j-1) = [p1_a,p2_a];
    Err_nA(j-1)  = norm(Nhat_nA(:,j) - A3);

    % ---------- B-spline (original) ----------
    tic;
    [G_bs, p1_bs,p2_bs,m1_bs,m2_bs,m3_bs,fval_bs,ef_bs] = GTMfun(a2, A3, C2);
    t_bs(j-1) = toc;

    Survived_bs = a2.*(1-maturity(p1_bs,p2_bs));
    Mort_bs     = mortality(m1_bs,m2_bs,m3_bs);
    Nhat_bs(:,j)= G_bs*(dt*exp(-12*Mort_bs).*Survived_bs - exp(-6*Mort_bs).*C2);
    GT2_bs(:,:,j-1) = G_bs;
    M2_bs(:,j-1) = [m1_bs,m2_bs,m3_bs];  P2_bs(:,j-1) = [p1_bs,p2_bs];
    Err_bs(j-1)  = norm(Nhat_bs(:,j) - A3);

    % ---------- Gamma (traditional) ----------
    tic;
    [G_g, p1_g,p2_g,m1_g,m2_g,m3_g, Linf_g,K_g,CV_g, fval_g,ef_g] = ...
        GTMfun_Gamma(a2, A3, C2, model.x_mid(:), model.x_edges(:));
    t_gamma(j-1) = toc;

    Survived_g = a2.*(1-maturity(p1_g,p2_g));
    Mort_g     = mortality(m1_g,m2_g,m3_g);
    Nhat_gam(:,j)= G_g*(dt*exp(-12*Mort_g).*Survived_g - exp(-6*Mort_g).*C2);
    GT2_gam(:,:,j-1) = G_g;
    M2_gam(:,j-1) = [m1_g,m2_g,m3_g];  P2_gam(:,j-1) = [p1_g,p2_g];
    Err_gam(j-1)  = norm(Nhat_gam(:,j) - A3);

end

% --------- summary ----------
fprintf('Median ||N - G*Nprev||:  B-spline=%.4g  |  Normal(Chen)=%.4g  |  Gamma=%.4g  |  LogNorm=%.4g\n', ...
    median(Err_bs), median(Err_nA), median(Err_gam), median(Err_lgn));
fprintf('Mean   ||N - G*Nprev||:  B-spline=%.4g  |  Normal(Chen)=%.4g  |  Gamma=%.4g  |  LogNorm=%.4g\n', ...
    mean(Err_bs), mean(Err_nA), mean(Err_gam), mean(Err_lgn));

% Build an N_obs matrix aligned with your predictions; here we stack A3 used each j
N_obs = zeros(nx, nt); N_obs(:,1)=N3(:,1);   % or leave first col zero/NaN if no obs at j=1
for j=2:nt
    [~,~,A3,~] = Initial_Data(year+j-2); % or cache A3 in the loop
    N_obs(:,j) = A3;
end

% plot (bars = observed; solid=B-spline; dashed=Normal; dotted=Gamma)
% compare_and_plot_all(N_obs, Nhat_bs, Nhat_nA, Nhat_gam, Nhat_lgn);

% Build per-year LaTeX table
%tex = make_error_table_tex(year_start, year_finish, Err_bs, Err_nA, Err_gam, Err_lgn, ...
   % 'filename','gtm_error_table.tex', ...
   % 'decimals',4, ...
   % 'caption','Per-year projection error (||N_y - G N_{y-1}||_2) for each GTM', ...
   % 'label','tab:gtm-errors', ...
   % 'booktabs', true);

% (Optional) peek at the first few lines in MATLAB
%disp(extractBetween(string(tex), 1, min(strlength(string(tex)), 400)));

fprintf("\nAverage run time (seconds) per estimation:\n");
fprintf("   B-spline:    %.4f sec\n", mean(t_bs));
fprintf("   Normal:      %.4f sec\n", mean(t_norm));
fprintf("   Gamma:       %.4f sec\n", mean(t_gamma));
fprintf("   Log-normal:  %.4f sec\n", mean(t_logn));

compare_and_plot_grid(N_obs, Nhat_bs, Nhat_nA, Nhat_gam, Nhat_lgn)
plot_error_timeseries(Err_bs, Err_nA, Err_gam, Err_lgn)
