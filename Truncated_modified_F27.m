%% FUNCTION 27 - PENALTY FUNCTION 1 (with specific truncated newton function useful for big dimensions)

% NON STA FUNZIONANDO

rng(345989);

F = @(x) F27(x);  % Definizione della funzione F27 come handle
JF = @(x) JF27(x,true,0); % Definizione della funzione JF27 come handle (derivata esatta)

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

[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1,converged1, violations1] = truncated_newton_27(x0, F, JF, true, 0, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag1 %  CONVERGE IN 39 ITERAZIONI
f1 
gradf_norm1 
last_bt=btseq1(k1-10:k1) 
last_cg=cgiterseq1(k1-10:k1)
%conv_ord1(k1-10:k1)
%violations1 % non ce ne sono

%% cambio tolgrad

tolgrad=5e-7;

[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1,converged1, violations1] = truncated_newton_27(x0, F, JF, true, 0, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag1 %  CONVERGE IN 41 ITERAZIONI
f1
gradf_norm1
last_bt=btseq1(k1-10:k1)  
last_cg=cgiterseq1(k1-10:k1) 
%conv_ord1(k1-10:k1)
