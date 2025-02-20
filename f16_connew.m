%% FUNCTION 16  (with different initial points)- with exact derivatives and finite differences

sparse=true;

F = @(x) F79(x);  % Defining F16 as function handle
JF_gen = @(x,exact,fin_dif2,h) JF79(x,exact,fin_dif2,h); % Defining JF16 as function handle ;
HF_gen= @(x,exact,fin_dif2,h) HF79New(x,sparse,exact,fin_dif2,h,JF); % Defining HF16 as function handle (sparse version)
load forcing_terms.mat % possible terms for adaptive tolerance 
format long 
%% n=10^3 (1e3)

rng(345989);

n=1e3
kmax=1.5e3; % maximum number of iterations of Newton method
tolgrad=5e-7; % tolerance on gradient norm

cg_maxit=50; % maximum number of iterations of coniugate gradient method (for the linear system)
z0=zeros(n,1); % initial point of coniugate gradient method (for the linear system)

% Backtracking parameters
c1=1e-4;
rho=0.50;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

x0 = ones(n, 1); % initial point
N=10; % number of initial points to be generated

% Initial points:
Mat_points=repmat(x0,1,N+1); 
rand_mat=2*rand(n, N)-1;
Mat_points(:,2:end)=Mat_points(:,2:end) + rand_mat; % matrix with columns=initial points


JF=@(x) JF_gen(x,true,false,1e-4);
HF=@(x) HF79New(x,true,false,false,1e-4,JF);
[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1_df1,flag1, converged1, violations1] = truncated_newton(Mat_points(:,8), F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
disp(f1)