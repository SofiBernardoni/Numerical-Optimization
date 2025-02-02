%% FUNCTION 79 - PENALTY FUNCTION 1

rng(345989);

F = @(x) F79(x);  % Definizione della funzione F79 come handle
JF = @(x) JF79(x,true,false,0); % Definizione della funzione JF79 come handle (derivata esatta)
HF= @(x) HF79(x,true,true,false,0); % Definizione della funzione HF79 come handle (derivata esatta)  % check if sparsity is ok

load forcing_terms.mat % termini per tolleranza adattiva

%% n=10^3 (1e3)

rng(345989);

n=1e3;
x0=-1*ones(n,1);

kmax=1e3;
tolgrad=1e-5;
cg_maxit=50;

z0=zeros(n,1);
c1=1e-4;
rho=0.50;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1,converged1, violations1] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag1 %CONVERGE IN 928 STEP
f1 
gradf_norm1 
last_bt=btseq1(end-10:end) %alla fine sempre passo unitarrio. IN REALTA LO FA SEMPRE (lo vedo con max(btseq1)=0)
last_cg=cgiterseq1(end-10:end) % fa sempre steepest descent max(cgiterseq1)=0. vedi violazione da ogni passo
conv_ord1(end-10:end) % ORDINE DI CONVERGENZA 1(ci oscilla intorno alla fine)
%violations1 % VIOLAZIONI AD OGNI PASSO (steepest)

%% abbasso tolgrad=5e-7 e kmax=1500

kmax=1500;
tolgrad=5e-7;

[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1, converged1, violations1] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag1 %CONVERGE IN 1095 STEP
f1 
gradf_norm1 
%k1
last_bt=btseq1(k1-10:k1) %max(btseq1)=0 --> sempre step=1
last_cg=cgiterseq1(k1-10:k1) %max(cgiterseq1)=0 --> sempre steepest descent (violations1=k2)
conv_ord1(k1-10:k1) % ORDINE DI CONVERGENZA 1
%violations1


%% n=10^4 (1e4)

rng(345989);

n=1e4;
x0=-1*ones(n,1);

kmax=1500;
tolgrad=1e-3;
cg_maxit=50;

z0=zeros(n,1);
c1=1e-4;
rho=0.5;
btmax=50; %  compatible with rho (with alpha0=1 you get min_step 8.8e-16)

[x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2,flag2, converged2, violations2] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag2 %%% ARMIJO not satis k2=38
f2
gradf_norm2 
last_bt=btseq2(k2-10:k2) 
last_cg=cgiterseq2(k2-10:k2)
%conv_ord2(k2-10:k2)
%violations2


%% rho=0.9, btmax=340

% stessa situa (anche per scelte intermedie)
rho=0.90;
btmax=340; % min_step=2,76e-16

[x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2,flag2, converged2, violations2] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag2 %%% ARMIJO not satis k2=325 con grad_norm=19
f2
gradf_norm2 
last_bt=btseq2(end-10:end) % sifanno passi piccolissimi bt=326 1.2e-15
last_cg=cgiterseq2(end-10:end)
% No violazioni

%% CAMBIO I PARAM di ARmijo --> diminuisco c=1e-8
c1=1e-8;
rho=0.5;
btmax=50; % min_step= 8,8e-16

[x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2,flag2,converged2, violations2] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag2 %%% ARMIJO not satis k2=38
f2
gradf_norm2 
last_bt=btseq2(end-10: end) % sempre max 
last_cg=cgiterseq2(end-10:end) %% è 2 inizialmente



%%
c1=1e-2;
rho=0.5;
btmax=50;

[x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2,flag2,converged2, violations2] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag2 %%% ARMIJO not satis k2=38
f2
gradf_norm2 
last_bt=btseq2(k2-10:k2) % sempre max 
last_cg=cgiterseq2(k2-10:k2) %% è 2 inizialmente

%% n=10^5 (1e5)

rng(345989);

n=1e5;
x0=-1*ones(n,1);

kmax=1500;
tolgrad=1e-3;
cg_maxit=50;

z0=zeros(n,1);
c1=1e-4;
rho=0.5;
btmax=50; %  compatible with rho (with alpha0=1 you get min_step 8.8e-16)

[x3, f3, gradf_norm3, k3, xseq3, btseq3,cgiterseq3,conv_ord3,flag3, converged3, violations3] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag3 % armijo non soddisfatto con param standard. (all'iteraz 1043)
f3
gradf_norm3 
last_bt=btseq3(end-10:end) % si arriva al max dei bt
last_cg=cgiterseq3(end-10:end)
%conv_ord3(end-10:end)
%violations3

%NB: anche con c1=1e-8

%% cambio param backtracking rho=0.8
rho=0.8;
btmax=158; %  compatible with rho (with alpha0=1 you get min_step 4.878e-16)

[x3, f3, gradf_norm3, k3, xseq3, btseq3,cgiterseq3,conv_ord3,flag3, converged3, violations3] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag3 % armijo non soddisfatto (all'iteraz 1348). gradiente ancora alto (57). Nello stesso punto di prima
f3
gradf_norm3 
last_bt=btseq3(end-10:end) % si arriva a 154. SI MUOVE DI POCHISSIMO 1e-15 (in 500 passi sarebbe movimento complessivo di 1e-13)
last_cg=cgiterseq3(end-10:end)
%conv_ord3(end-10:end)
%violations3

%% cambio param backtracking c1=1e-2

c1=1e-2;
rho=0.5;
btmax=50; %  compatible with rho (with alpha0=1 you get min_step 8.8e-16)

[x3, f3, gradf_norm3, k3, xseq3, btseq3,cgiterseq3,conv_ord3,flag3,converged3, violations3] = truncated_newton(x0, F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
flag3 % armijo non soddisfatto con param standard. (all'iteraz 1043)
f3
gradf_norm3 
last_bt=btseq3(end-10:end) % si arriva al max dei bt
last_cg=cgiterseq3(end-10:end)
%conv_ord3(end-10:end)
%violations3
