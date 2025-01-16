%% FUNCTION 79 - PENALTY FUNCTION 1

rng(345989);

F = @(x) F79(x);  % Definizione della funzione F79 come handle
JF = @(x) JF79(x); % Definizione della funzione JF79 come handle
HF= @(x) HF79(x); % Definizione della funzione HF79 come handle

load forcing_terms.mat % termini per tolleranza adattiva

%% n=10^3 (1e3)

rng(345989);

n=1e3;
x0=-1*ones(n,1);

kmax=500;
tolgrad=1e-8;
cg_maxit=50;

z0=zeros(n,1);
c1=1e-4;
rho=0.90;
btmax=320; % compatible with rho (with alpha0=1 you get min_step 1e-14)

[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1, violations1] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag1 
f1 
gradf_norm1 
%k1
last_bt=btseq1(k1-10:k1) 
last_cg=cgiterseq1(k1-10:k1)
%conv_ord1(k1-10:k1)
%violations1

%NOTA: con parametri standard armijo non soddisfatta
%rho=0.90;
%btmax=320;
% --> k=8 armijo lo blocca
