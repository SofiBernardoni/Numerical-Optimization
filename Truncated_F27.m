%% FUNCTION 27 - PENALTY FUNCTION 1

rng(345989);

F = @(x) F27(x);  % Definizione della funzione F27 come handle
JF = @(x) JF27(x); % Definizione della funzione JF27 come handle
HF= @(x) HF27(x,true); % Definizione della funzione HF27 come handle  % check if sparsity is ok

load forcing_terms.mat % termini per tolleranza adattiva

%% n=10^3 (1e3)

rng(345989);

n=1e3; 
x0=(1:n)';  % Punto iniziale

kmax=1e3;
tolgrad=1e-5;
cg_maxit=50;

z0=zeros(n,1);
c1=1e-4;
rho=0.5;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1, violations1] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag1 %  CONVERGE IN 39 ITERAZIONI
f1 
gradf_norm1 
last_bt=btseq1(k1-10:k1) 
last_cg=cgiterseq1(k1-10:k1)
%conv_ord1(k1-10:k1)
%violations1 % non ce ne sono

%% cambio tolgrad

tolgrad=5e-7;

[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1, violations1] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag1 %  CONVERGE IN 41 ITERAZIONI
f1
gradf_norm1
last_bt=btseq1(k1-10:k1)  
last_cg=cgiterseq1(k1-10:k1) 
%conv_ord1(k1-10:k1)

%%  n=10^4 (1e4)

rng(345989);

n=1e4; 
x0=(1:n)';  % Punto iniziale

kmax=1e3;
tolgrad=1e-5;
cg_maxit=50;

z0=zeros(n,1);
c1=1e-4;
rho=0.5;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

[x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2,flag2, violations2] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag2 %  CONVERGE IN 46 ITERAZIONI
f2
gradf_norm2 
last_bt=btseq2(k2-10:k2) 
last_cg=cgiterseq2(k2-10:k2)
%conv_ord2(k2-10:k2)
%violations2  % non ce ne sono

%% cambio tolgrad

tolgrad=5e-7;
[x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2,flag2, violations2] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag2 %  CONVERGE IN 47 ITERAZIONI
f2
gradf_norm2 
last_bt=btseq2(k2-10:k2) 
last_cg=cgiterseq2(k2-10:k2)
%conv_ord2(k2-10:k2)
%violations2 


%%  n=10^5 (1e5)
%v ATTENZIONE: ANCORA PROBLEMA SPARSITà

rng(345989);

n=1e5; 
x0=(1:n)';  % Punto iniziale

kmax=1e3;
tolgrad=1e-5;
cg_maxit=50;

z0=zeros(n,1);
c1=1e-4;
rho=0.5;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

[x3, f3, gradf_norm3, k3, xseq3, btseq3,cgiterseq3,conv_ord3,flag3, violations3] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag3 
f3
gradf_norm3
last_bt=btseq3(end-10:end) 
last_cg=cgiterseq3(end-10:end)
%conv_ord3(end-10:end)
%violations3  
