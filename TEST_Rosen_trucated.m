%% ROSENBROCK FUNCTION

f_Rosen = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2 ;
gradf_Rosen= @(x) [400*x(1)^3+(2-400*x(2))-2;
                   200*(x(2)-x(1)^2)];
Hessf_Rosen=@(x) [1200*x(1)^2-400*x(2)+2, -400*x(1) ;
                    -400*x(1), 200];
x0_a=[1.2;1.2];
x0_b=[-1.2;1];

load forcing_terms.mat

%% TEST of Trucated Newton without preconditioning with exact derivatives with con x0=[1.2;1.2]

kmax=500;
tolgrad=1e-8;
cg_maxit=50;

z0=[0;0];
c1=1e-4;
rho=0.9;
btmax=200; % compatible with rho

% con x0=[1.2;1.2]
[x1, f1, gradf1_norm, k1, xseq1, btseq1,cgiterseq1,flag1,v1] = truncated_newton(x0_a, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
% anche con rho=0.9 e btmax=200 termina alla terza iterazione per Armijo (prima avevo messo rho piu piccolo)
[x2, f2, gradf2_norm, k2, xseq2, btseq2,cgiterseq2,flag2,v2] = truncated_newton(x0_a, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
% ottengo gli stessi risultati variando termine per tolleranza adattiva

%% TEST of Trucated Newton without preconditioning with exact derivatives with con x0=[-1.2;1]
kmax=500;
tolgrad=1e-8;
cg_maxit=50;

z0=[0;0];
c1=1e-4;
rho=0.8;
btmax=130;

[x3, f3, gradf3_norm, k3, xseq3, btseq3,cgiterseq3,flag3,v3] = truncated_newton(x0_b, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);

[x4, f4, gradf4_norm, k4, xseq4, btseq4,cgiterseq4,flag4,v4] = truncated_newton(x0_b, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);

% con c1=1e-4
% in entrambe con kmax=100 fa tutte le iterazioni e cg in 2 passi circa, mentre backtracking si assesta su 43/ 44 iter con rho=0.9
% con kmax=500 : ancora kmax
% con rho=0.5 o 0.7 bloccato per Armijo al secondo passo entrambe
% con rho= 0.8 arriva a kmax=500. MA si bt=20/21---> migliore finora

%con c1=0.1: Armijo dopo 7 iterazioni NON buon risultato

flag3
flag4
%% sezione tuning on x0=[-1.2;1]
kmax=500;
tolgrad=1e-8;
cg_maxit=50;

z0=[0;0];
c1=1e-3;
rho=0.8;
btmax=130;
[x3, f3, gradf3_norm, k3, xseq3, btseq3,cgiterseq3,flag3,v3] = truncated_newton(x0_b, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
x3
f3
gradf3_norm
k3
last_bt=btseq3(k3-10:k3)
last_cg=cgiterseq3(k3-10:k3)
flag3
v3
%migliore combo kmax=500

plot_iterative_optimization_results(f_Rosen,xseq3, btseq3);

