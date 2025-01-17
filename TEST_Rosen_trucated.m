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

kmax=1500;
tolgrad=1e-5;
cg_maxit=50;

z0=[0;0];
c1=1e-4;
rho=0.5;
btmax=50; % compatible with rho


% con x0=[1.2;1.2]
[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,c1,flag1,v1] = truncated_newton(x0_a, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag1 %max num it 
f1
gradf_norm1
last_bt=btseq1(k1-10:k1) % bt=49 mi sto spostando di poco
last_cg=cgiterseq1(k1-10:k1)% violato sempre curvature condition (tranne in 1), quasi sempre steepest descent max(cgiterseq1)=1
%conv_ord1(k1-10:k1) 
%violations1

%%

[x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,c2,flag2,v2] = truncated_newton(x0_a, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
flag2
f2
gradf_norm2
last_bt=btseq2(end-10:end) 
last_cg=cgiterseq2(end-10:end)
%conv_ord2(end-10:end) 
%violations2
%%%%%%%%%%%%%%%% ottengo gli stessi risultati variando termine per tolleranza adattiva

%% TEST of Trucated Newton without preconditioning with exact derivatives with con x0=[-1.2;1]
kmax=500;
tolgrad=1e-5;
cg_maxit=50;

z0=[0;0];
c1=1e-4;
rho=0.5;
btmax=50;
% rho=0.8;
% btmax=130;

[x3, f3, gradf_norm3, k3, xseq3, btseq3,cgiterseq3,c3,flag3,v3] = truncated_newton(x0_b, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag3 %armijo non soddsfatto a k3=1
f3
gradf_norm3
last_bt=btseq3(end-10:end) 
last_cg=cgiterseq3(end-10:end)
%conv_ord3(end-10:end) 
%violations3

%[x4, f4, gradf4_norm, k4, xseq4, btseq4,cgiterseq4,c4,flag4,v4] = truncated_newton(x0_b, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);

% con c1=1e-4
% in entrambe con kmax=100 fa tutte le iterazioni e cg in 2 passi circa, mentre backtracking si assesta su 43/ 44 iter con rho=0.9
% con kmax=500 : ancora kmax
% con rho=0.5 o 0.7 bloccato per Armijo al secondo passo entrambe
% con rho= 0.8 arriva a kmax=500. MA si bt=20/21---> migliore finora

%con c1=0.1: Armijo dopo 7 iterazioni NON buon risultato

flag3
flag4
%% sezione tuning on x0=[-1.2;1]
kmax=1500;
rho=0.85;
btmax=220;
[x3, f3, gradf_norm3, k3, xseq3, btseq3,cgiterseq3,c3,flag3,v3] = truncated_newton(x0_b, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag3 % CONVERGE in 868 step
f3
gradf_norm3
last_bt=btseq3(end-10:end) % alla fine passi ragionevli bt=28
last_cg=cgiterseq3(end-10:end) % sempre tra 0 e 2
%c3(end-10:end) %% oscilla!!
%v3 % 1 SOLA VIOLAZIONE

%plot_iterative_optimization_results(f_Rosen,xseq3, btseq3);

[x4, f4, gradf_norm4, k4, xseq4, btseq4,cgiterseq4,c4,flag4,v4] = truncated_newton(x0_b, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
flag4 % CONVERGE in 1018 step
f4
gradf_norm4
last_bt=btseq4(end-10:end) % alla fine passi ragionevli bt=28
last_cg=cgiterseq4(end-10:end) % sempre tra 0 e 2
%c4(end-10:end) %% oscilla!!
%v4 % 2 SOLA VIOLAZIONE

% vedi se diminuire rho tra 0.5 e 0.85