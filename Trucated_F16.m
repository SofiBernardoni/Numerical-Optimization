%% FUNCTION 16 - BANDED TRIGONOMETRIC PROBLEM

rng(345989);

F = @(x) F16(x);  % Definizione della funzione F16 come handle
JF = @(x) JF16(x); % Definizione della funzione JF16 come handle
HF= @(x) HF16(x); % Definizione della funzione HF16 come handle

load forcing_terms.mat % termini per tolleranza adattiva

%% n=10^3 (1e3)

rng(345989);

n=1e3; 
x0 = ones(n, 1);  % Punto iniziale

kmax=1e3;
tolgrad=1e-8;
cg_maxit=50;

z0=zeros(n,1);
c1=1e-4;
rho=0.5;
btmax=47; % compatible with rho (with alpha0=1 you get min_step 1e-14)

[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1, violations1] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag1 % MAX IT (stessi risultati con kmax=500 o kmax=1e3)
f1
gradf_norm1
%k1
last_bt=btseq1(k1-10:k1)
last_cg=cgiterseq1(k1-10:k1)
%conv_ord1(k1-10:k1)
%violations1

% btseq=45 per tutti gli ultimi (come nel caso con btmax=45). Anche alzando btmax=47 (ancora compatibile)
% ctseq=25 gli ultimi

%fterms_quad --> armijo not satisfied

%% cambio rho e btmax

rho=0.3;
btmax=28; 
[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1, violations1] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag1 %MAX IT
f1
gradf_norm1

last_bt=btseq1(k1-10:k1)  %26 per tutti gli ultimi
last_cg=cgiterseq1(k1-10:k1) %sempre 25

%% aumento c1=1e-1

% ottengo ancora max it con bt=45 e cg=25

kmax=1e3;
tolgrad=1e-8;
cg_maxit=50;

z0=zeros(n,1);
c1=0.1; %%% parametro variato %%
rho=0.5;
btmax=47; % compatible with rho (with alpha0=1 you get min_step 1e-14)

[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1, violations1] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag1 % MAX IT (stessi risultati con kmax=500 o kmax=1e3)
f1
gradf_norm1
%k1
last_bt=btseq1(k1-10:k1)
last_cg=cgiterseq1(k1-10:k1)

%% dimuinuisco c1=1e-10

%stesso risultato

kmax=1e3;
tolgrad=1e-8;
cg_maxit=50;

z0=zeros(n,1);
c1=1e-10; %%% parametro variato %%
rho=0.5;
btmax=47; % compatible with rho (with alpha0=1 you get min_step 1e-14)

[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1, violations1] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag1 % MAX IT (stessi risultati con kmax=500 o kmax=1e3)
f1
gradf_norm1
%k1
last_bt=btseq1(k1-10:k1)
last_cg=cgiterseq1(k1-10:k1)





%% n=10^4 (1e4)

rng(345989);

n=1e4; 
x0 = ones(n, 1);  % Punto iniziale

kmax=1e3;
tolgrad=1e-8;
cg_maxit=50;

z0=zeros(n,1);
c1=1e-4;
rho=0.5;
btmax=47; % compatible with rho (with alpha0=1 you get min_step 1e-14)

[x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2,flag2, violations2] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag2 % armijo not satis k2=29
f2
gradf_norm2
%k2
last_bt=btseq2(k2-10:k2) % non si fissa
last_cg=cgiterseq2(k2-10:k2) % alla fine 18
%conv_ord1(k1-10:k1)
%violations1
%% change rho=0.7 btmax=95

% come prima armijo not satis k2=54 
c1=1e-4;
rho=0.7;
btmax=95; % compatible with rho (with alpha0=1 you get min_step 1e-14)

[x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2,flag2, violations2] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag2 % armijo not satis k2=54 
f2
gradf_norm2
%k2
last_bt=btseq2(k2-10:k2) % non si fissa
last_cg=cgiterseq2(k2-10:k2) % alla fine 18
%conv_ord1(k1-10:k1)

%% change rho=0.85 btmax=200

% come prima armijo not satis k2=39. gli altri bt=187--> si muove di 1e-14
% (si muove lentissimo)

c1=1e-4;
rho=0.85;
btmax=200; % compatible with rho (with alpha0=1 you get min_step 1e-14)

[x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2,flag2, violations2] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag2 % armijo not satis k2=54 
f2
gradf_norm2
%k2
last_bt=btseq2(k2-10:k2) % non si fissa
last_cg=cgiterseq2(k2-10:k2) % alla fine 18
%conv_ord1(k1-10:k1)
