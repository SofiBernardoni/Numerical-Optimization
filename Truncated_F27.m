%% FUNCTION 27 - PENALTY FUNCTION 1

rng(345989);

F = @(x) F27(x);  % Definizione della funzione F27 come handle
JF = @(x) JF27(x); % Definizione della funzione JF27 come handle
HF= @(x) HF27(x); % Definizione della funzione HF27 come handle

load forcing_terms.mat % termini per tolleranza adattiva

%% n=10^3 (1e3)

rng(345989);

n=1e3; 
x0=(1:n)';  % Punto iniziale

kmax=500;
tolgrad=1e-8;
cg_maxit=50;

z0=zeros(n,1);
c1=1e-4;
rho=0.5;
btmax=47; % compatible with rho (with alpha0=1 you get min_step 1e-14)

[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1, violations1] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag1 % MAX IT 
f1 
gradf_norm1 %%%%%%%%%%%%%%%% FA SCHIFO 1E-4
%k1
last_bt=btseq1(k1-10:k1) %%% VIENE BT=46
last_cg=cgiterseq1(k1-10:k1) %%% VIENE CG=1 (sempre curvatura negativa violata!!)
%conv_ord1(k1-10:k1)
%violations1

%% CAMBIO KMAX PER VEDERE DA QUANTO STA LI----> è lì da iter 200

kmax=200;

[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1, violations1] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag1 % MAX IT 
f1 
gradf_norm1 %%%%%%%%%%%%%%%% FA SCHIFO 1E-4
%k1
last_bt=btseq1(k1-10:k1) %%% VIENE BT=46
last_cg=cgiterseq1(k1-10:k1) 