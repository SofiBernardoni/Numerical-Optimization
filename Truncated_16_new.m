%% FUNCTION 16 - BANDED TRIGONOMETRIC PROBLEM (with different initial points)

% Utilizza un unico file alla fine

rng(345989);

F = @(x) F16(x);  % Definizione della funzione F16 come handle
JF = @(x) JF16(x,true,0); % Definizione della funzione JF16 come handle (derivata esatta)
HF= @(x) HF16(x,true,true,0); % Definizione della funzione HF16 come handle (derivata esatta) % Sparse version

load forcing_terms.mat % termini per tolleranza adattiva

%% n=10^3 (1e3)

rng(345989);

n=1e3; 

kmax=1.5e3;
tolgrad=5e-7; 
cg_maxit=50;

z0=zeros(n,1);
c1=1e-4;
rho=0.5;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

x0 = ones(n, 1);  % Initial point
N=10; % number of initial points to be generated

% Initial points:
Mat_points=repmat(x0,1,N+1); 
rand_mat=2*rand(n, N)-1;
Mat_points(:,2:end)=Mat_points(:,2:end) + rand_mat; % matrix with columns=initial points

vec_times1=zeros(1,N+1); % vector with execution times
vec_val1=zeros(1,N+1); %vector with minimal values found
vec_grad1=zeros(1,N+1); %vector with final gradient
vec_iter1=zeros(1,N+1); %vector with number of iterations 

% INSERIRE ORDINE CONVERGENZA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
vec_converged1=zeros(1,N+1); %vector of booleans (true if it has converged) 

vec_violations1=zeros(1,N+1); % per vedere quante violazioni in gradiente coniugato

for j =1:N+1
    
    tic;
    [x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1, converged1, violations1] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
    vec_times1=toc;

    disp(['Tentativo n. ',num2str(j),': ',flag1]) % introdurre conteggio fallimenti/successi
    vec_converged1(j)=converged1;
    last_bt=btseq1(end-10:end) ; % salvare??
    last_cg=cgiterseq1(end-10:end) ; % salvare??
    %conv_ord1(end-10:end) %aggiustare
    vec_val1(j)=f1;
    vec_grad1(j)=gradf_norm1;
    vec_iter1(j)=k1;
    vec_violations1(j)=violations1;
end

% Forse ha senso plottare poi solo i risultati delle convergenze
% per confrontare i metodi sulle varie dimensioni e enlle varianti ha senso
% usare per esempio la media e le statistiche sui vari successi ottenuti (tipo la media delle iterazioni e del tempo)

% INSERIRE TABELLA
% INSERIRE GRAFICI

%% n=10^4 (1e4)

rng(345989);

n=1e4; 

kmax=1.5e3;
tolgrad=5e-7; 
cg_maxit=50;

z0=zeros(n,1);
c1=1e-4;
rho=0.5;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

x0 = ones(n, 1);  % Initial point
N=10; % number of initial points to be generated

% Initial points:
Mat_points=repmat(x0,1,N+1); 
rand_mat=2*rand(n, N)-1;
Mat_points(:,2:end)=Mat_points(:,2:end) + rand_mat; % matrix with columns=initial points

vec_times2=zeros(1,N+1); % vector with execution times
vec_val2=zeros(1,N+1); %vector with minimal values found
vec_grad2=zeros(1,N+1); %vector with final gradient
vec_iter2=zeros(1,N+1); %vector with number of iterations 

% INSERIRE ORDINE CONVERGENZA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
vec_converged2=zeros(1,N+1); %vector of booleans (true if it has converged) 

vec_violations2=zeros(1,N+1); % per vedere quante violazioni in gradiente coniugato

for j =1:N+1
    
    tic;
    [x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2,flag2, converged2, violations2] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
    vec_times2=toc;

    disp(['Tentativo n. ',num2str(j),': ',flag2]) % introdurre conteggio fallimenti/successi
    vec_converged2(j)=converged2;
    last_bt=btseq2(end-10:end) ; % salvare??
    last_cg=cgiterseq2(end-10:end) ; % salvare??
    %conv_ord2(end-10:end) %aggiustare
    vec_val2(j)=f2;
    vec_grad2(j)=gradf_norm2;
    vec_iter2(j)=k2;
    vec_violations2(j)=violations2;
end

% Forse ha senso plottare poi solo i risultati delle convergenze
% per confrontare i metodi sulle varie dimensioni e enlle varianti ha senso
% usare per esempio la media e le statistiche sui vari successi ottenuti (tipo la media delle iterazioni e del tempo)

% INSERIRE TABELLA
% INSERIRE GRAFICI

