%% FUNCTION 79  (with different initial points)- with exact derivatives and finite differences- QUADRATIC TERM OF CONVERGENCE

sparse=true;

F = @(x) F79(x);  % Defining F79 as function handle
JF_gen = @(x,exact,fin_dif2,h) JF79(x,exact,fin_dif2,h); % Defining JF79 as function handle 
HF_gen= @(x,exact,fin_dif2,h) HF79(x,sparse,exact,fin_dif2,h); % Defining HF79 as function handle (sparse version)

load forcing_terms.mat % possible terms for adaptive tolerance 

%% n=10^3 (1e3)

rng(345989);

n=1e3; 

kmax=1.5e3; % maximum number of iterations of Newton method
tolgrad=5e-7; % tolerance on gradient norm

cg_maxit=50; % maximum number of iterations of coniugate gradient method (for the linear system)
z0=zeros(n,1); % initial point of coniugate gradient method (for the linear system)

% Backtracking parameters
c1=1e-4;
rho=0.50;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

x0=-1*ones(n,1);  % initial point
N=10; % number of initial points to be generated

% Initial points:
Mat_points=repmat(x0,1,N+1); 
rand_mat=2*rand(n, N)-1;
Mat_points(:,2:end)=Mat_points(:,2:end) + rand_mat; % matrix with columns=initial points

% Structure for EXACT derivatives
vec_times1_ex=zeros(1,N+1); % vector with execution times 
vec_val1_ex=zeros(1,N+1); %vector with minimal values found
vec_grad1_ex=zeros(1,N+1); %vector with final gradient
vec_iter1_ex=zeros(1,N+1); %vector with number of iterations 
% INSERIRE ORDINE CONVERGENZA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
vec_converged1_ex=zeros(1,N+1); % vector of booleans (true if it has converged) 
vec_violations1_ex=zeros(1,N+1); % vector with number of violations of curvature condition in Newton method

JF_ex = @(x) JF_gen(x,true,false,0);
HF_ex = @(x) HF_gen(x,true,false,0);

% Structure for derivatives approximated with FINITE DIFFERENCES (classical version) 
mat_times1_fd1=zeros(6,N+1); % matrix with execution times 
mat_val1_fd1=zeros(6,N+1); %matrix with minimal values found
mat_grad1_fd1=zeros(6,N+1); %matrix with final gradient
mat_iter1_fd1=zeros(6,N+1); %matrix with number of iterations 
% INSERIRE ORDINE CONVERGENZA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mat_converged1_fd1=zeros(6,N+1); % matrix of booleans (true if it has converged) 
mat_violations1_fd1=zeros(6,N+1); % matrix with number of violations of curvature condition in Newton method

JF_fd1 = @(x,h) JF_gen(x,false,false,h);
HF_fd1 = @(x,h) HF_gen(x,false,false,h);

% Structure for derivatives approximated with FINITE DIFFERENCES (version with h=h*abs(x_j) as increment) 
mat_times1_fd2=zeros(6,N+1); % matrix with execution times 
mat_val1_fd2=zeros(6,N+1); %matrix with minimal values found
mat_grad1_fd2=zeros(6,N+1); %matrix with final gradient
mat_iter1_fd2=zeros(6,N+1); %matrix with number of iterations 
% INSERIRE ORDINE CONVERGENZA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mat_converged1_fd2=zeros(6,N+1); % matrix of booleans (true if it has converged) 
mat_violations1_fd2=zeros(6,N+1); % matrix with number of violations of curvature condition in Newton method

JF_fd2 = @(x,h) JF_gen(x,false,true,h);
HF_fd2 = @(x,h) HF_gen(x,false,true,h);

for j =1:N+1
    disp(['Condizione iniziale n. ',num2str(j)])

    % EXACT DERIVATIVES
    tic;
    [x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1, converged1, violations1] = truncated_newton(Mat_points(:,j), F, JF_ex, HF_ex, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
    vec_times1_ex(j)=toc;

    disp(['Exact derivatives: ',flag1]) 
    vec_converged1_ex(j)=converged1;
    %conv_ord1(end-10:end) %aggiustare
    vec_val1_ex(j)=f1;
    vec_grad1_ex(j)=gradf_norm1;
    vec_iter1_ex(j)=k1;
    vec_violations1_ex(j)=violations1;
    
    for i=2:2:12
    h=10^(-i);
    
    % FINITE DIFFERENCES 1
    JF=@(x)JF_fd1(x,h);
    HF=@(x)HF_fd1(x,h);
    tic;
    [x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1, converged1, violations1] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
    mat_times1_fd1(i/2,j)=toc;

    disp(['Finite differences (classical version) with h=1e-',num2str(i),' : ',flag1]) 
    mat_converged1_fd1(i/2,j)=converged1;
    %conv_ord1(end-10:end) %aggiustare
    mat_val1_fd1(i/2,j)=f1;
    mat_grad1_fd1(i/2,j)=gradf_norm1;
    mat_iter1_fd1(i/2,j)=k1;
    mat_violations1_fd1(i/2,j)=violations1;


    % FINITE DIFFERENCES 2
    JF=@(x) JF_fd2(x,h);
    HF=@(x) HF_fd2(x,h);
    tic;
    [x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1, converged1, violations1] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
    mat_times1_fd2(i/2,j)=toc;

    disp(['Finite differences (new version) with h=1e-',num2str(i),' : ',flag1]) 
    mat_converged1_fd2(i/2,j)=converged1;
    %conv_ord1(end-10:end) %aggiustare
    mat_val1_fd2(i/2,j)=f1;
    mat_grad1_fd2(i/2,j)=gradf_norm1;
    mat_iter1_fd2(i/2,j)=k1;
    mat_violations1_fd2(i/2,j)=violations1;

    end
end

% Forse ha senso plottare poi solo i risultati delle convergenze
% per confrontare i metodi sulle varie dimensioni e enlle varianti ha senso
% usare per esempio la media e le statistiche sui vari successi ottenuti (tipo la media delle iterazioni e del tempo)

% INSERIRE TABELLA
% INSERIRE GRAFICI

%% n=10^4 (1e4)

rng(345989);

n=1e4; 

kmax=1.5e3; % maximum number of iterations of Newton method
tolgrad=5e-7; % tolerance on gradient norm

cg_maxit=100; % maximum number of iterations of coniugate gradient method (for the linear system)
z0=zeros(n,1); % initial point of coniugate gradient method (for the linear system)

% Backtracking parameters
c1=1e-4;
rho=0.50;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

x0=-1*ones(n,1);  % initial point
N=10; % number of initial points to be generated

% Initial points:
Mat_points=repmat(x0,1,N+1); 
rand_mat=2*rand(n, N)-1;
Mat_points(:,2:end)=Mat_points(:,2:end) + rand_mat; % matrix with columns=initial points

% Structure for EXACT derivatives
vec_times2_ex=zeros(1,N+1); % vector with execution times 
vec_val2_ex=zeros(1,N+1); %vector with minimal values found
vec_grad2_ex=zeros(1,N+1); %vector with final gradient
vec_iter2_ex=zeros(1,N+1); %vector with number of iterations 
% INSERIRE ORDINE CONVERGENZA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
vec_converged2_ex=zeros(1,N+1); % vector of booleans (true if it has converged) 
vec_violations2_ex=zeros(1,N+1); % vector with number of violations of curvature condition in Newton method

JF_ex = @(x) JF_gen(x,true,false,0);
HF_ex = @(x) HF_gen(x,true,false,0);

% Structure for derivatives approximated with FINITE DIFFERENCES (classical version) 
mat_times2_fd1=zeros(6,N+1); % matrix with execution times 
mat_val2_fd1=zeros(6,N+1); %matrix with minimal values found
mat_grad2_fd1=zeros(6,N+1); %matrix with final gradient
mat_iter2_fd1=zeros(6,N+1); %matrix with number of iterations 
% INSERIRE ORDINE CONVERGENZA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mat_converged2_fd1=zeros(6,N+1); % matrix of booleans (true if it has converged) 
mat_violations2_fd1=zeros(6,N+1); % matrix with number of violations of curvature condition in Newton method

JF_fd1 = @(x,h) JF_gen(x,false,false,h);
HF_fd1 = @(x,h) HF_gen(x,false,false,h);

% Structure for derivatives approximated with FINITE DIFFERENCES (version with h=h*abs(x_j) as increment) 
mat_times2_fd2=zeros(6,N+1); % matrix with execution times 
mat_val2_fd2=zeros(6,N+1); %matrix with minimal values found
mat_grad2_fd2=zeros(6,N+1); %matrix with final gradient
mat_iter2_fd2=zeros(6,N+1); %matrix with number of iterations 
% INSERIRE ORDINE CONVERGENZA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mat_converged2_fd2=zeros(6,N+1); % matrix of booleans (true if it has converged) 
mat_violations2_fd2=zeros(6,N+1); % matrix with number of violations of curvature condition in Newton method

JF_fd2 = @(x,h) JF_gen(x,false,true,h);
HF_fd2 = @(x,h) HF_gen(x,false,true,h);

for j =1:N+1
    disp(['Condizione iniziale n. ',num2str(j)])

    % EXACT DERIVATIVES
    tic;
    [x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2,flag2, converged2, violations2] = truncated_newton(Mat_points(:,j), F, JF_ex, HF_ex, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
    vec_times2_ex(j)=toc;

    disp(['Exact derivatives: ',flag2]) 
    vec_converged2_ex(j)=converged2;
    %conv_ord2(end-10:end) %aggiustare
    vec_val2_ex(j)=f2;
    vec_grad2_ex(j)=gradf_norm2;
    vec_iter2_ex(j)=k2;
    vec_violations2_ex(j)=violations2;
    
    for i=2:2:12
    h=10^(-i);
    
    % FINITE DIFFERENCES 1
    JF=@(x)JF_fd1(x,h);
    HF=@(x)HF_fd1(x,h);
    tic;
    [x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2,flag2, converged2, violations2] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
    mat_times2_fd1(i/2,j)=toc;

    disp(['Finite differences (classical version) with h=1e-',num2str(i),' : ',flag2]) 
    mat_converged2_fd1(i/2,j)=converged2;
    %conv_ord2(end-10:end) %aggiustare
    mat_val2_fd1(i/2,j)=f2;
    mat_grad2_fd1(i/2,j)=gradf_norm2;
    mat_iter2_fd1(i/2,j)=k2;
    mat_violations2_fd1(i/2,j)=violations2;


    % FINITE DIFFERENCES 2
    JF=@(x) JF_fd2(x,h);
    HF=@(x) HF_fd2(x,h);
    tic;
    [x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2,flag2, converged2, violations2] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
    vec_times2_fd2(i/2,j)=toc;

    disp(['Finite differences (new version) with h=1e-',num2str(i),' : ',flag2]) 
    vec_converged2_fd2(i/2,j)=converged2;
    %conv_ord2(end-10:end) %aggiustare
    vec_val2_fd2(i/2,j)=f2;
    vec_grad2_fd2(i/2,j)=gradf_norm2;
    vec_iter2_fd2(i/2,j)=k2;
    vec_violations2_fd2(i/2,j)=violations2;

    end
end

% Forse ha senso plottare poi solo i risultati delle convergenze
% per confrontare i metodi sulle varie dimensioni e enlle varianti ha senso
% usare per esempio la media e le statistiche sui vari successi ottenuti (tipo la media delle iterazioni e del tempo)

% INSERIRE TABELLA
% INSERIRE GRAFICI


%% n=10^5 (1e5)

rng(345989);

n=1e5; 

kmax=1.5e3; % maximum number of iterations of Newton method
tolgrad=5e-7; % tolerance on gradient norm

cg_maxit=100; % maximum number of iterations of coniugate gradient method (for the linear system)
z0=zeros(n,1); % initial point of coniugate gradient method (for the linear system)

% Backtracking parameters
c1=1e-4;
rho=0.50;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

x0=-1*ones(n,1);  % initial point
N=10; % number of initial points to be generated

% Initial points:
Mat_points=repmat(x0,1,N+1); 
rand_mat=2*rand(n, N)-1;
Mat_points(:,2:end)=Mat_points(:,2:end) + rand_mat; % matrix with columns=initial points

% Structure for EXACT derivatives
vec_times3_ex=zeros(1,N+1); % vector with execution times 
vec_val3_ex=zeros(1,N+1); %vector with minimal values found
vec_grad3_ex=zeros(1,N+1); %vector with final gradient
vec_iter3_ex=zeros(1,N+1); %vector with number of iterations 
% INSERIRE ORDINE CONVERGENZA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
vec_converged3_ex=zeros(1,N+1); % vector of booleans (true if it has converged) 
vec_violations3_ex=zeros(1,N+1); % vector with number of violations of curvature condition in Newton method

JF_ex = @(x) JF_gen(x,true,false,0);
HF_ex = @(x) HF_gen(x,true,false,0);

% Structure for derivatives approximated with FINITE DIFFERENCES (classical version) 
mat_times3_fd1=zeros(6,N+1); % matrix with execution times 
mat_val3_fd1=zeros(6,N+1); %matrix with minimal values found
mat_grad3_fd1=zeros(6,N+1); %matrix with final gradient
mat_iter3_fd1=zeros(6,N+1); %matrix with number of iterations 
% INSERIRE ORDINE CONVERGENZA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mat_converged3_fd1=zeros(6,N+1); % matrix of booleans (true if it has converged) 
mat_violations3_fd1=zeros(6,N+1); % matrix with number of violations of curvature condition in Newton method

JF_fd1 = @(x,h) JF_gen(x,false,false,h);
HF_fd1 = @(x,h) HF_gen(x,false,false,h);

% Structure for derivatives approximated with FINITE DIFFERENCES (version with h=h*abs(x_j) as increment) 
mat_times3_fd2=zeros(6,N+1); % matrix with execution times 
mat_val3_fd2=zeros(6,N+1); %matrix with minimal values found
mat_grad3_fd2=zeros(6,N+1); %matrix with final gradient
mat_iter3_fd2=zeros(6,N+1); %matrix with number of iterations 
% INSERIRE ORDINE CONVERGENZA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mat_converged3_fd2=zeros(6,N+1); % matrix of booleans (true if it has converged) 
mat_violations3_fd2=zeros(6,N+1); % matrix with number of violations of curvature condition in Newton method

JF_fd2 = @(x,h) JF_gen(x,false,true,h);
HF_fd2 = @(x,h) HF_gen(x,false,true,h);

for j =1:N+1
    disp(['Condizione iniziale n. ',num2str(j)])

    % EXACT DERIVATIVES
    tic;
    [x3, f3, gradf_norm3, k3, xseq3, btseq3,cgiterseq3,conv_ord3,flag3, converged3, violations3] = truncated_newton(Mat_points(:,j), F, JF_ex, HF_ex, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
    vec_times3_ex(j)=toc;

    disp(['Exact derivatives: ',flag3]) 
    vec_converged3_ex(j)=converged3;
    %conv_ord3(end-10:end) %aggiustare
    vec_val3_ex(j)=f3;
    vec_grad3_ex(j)=gradf_norm3;
    vec_iter3_ex(j)=k3;
    vec_violations3_ex(j)=violations3;
    
    for i=2:2:12
    h=10^(-i);
    
    % FINITE DIFFERENCES 1
    JF=@(x)JF_fd1(x,h);
    HF=@(x)HF_fd1(x,h);
    tic;
    [x3, f3, gradf_norm3, k3, xseq3, btseq3,cgiterseq3,conv_ord3,flag3, converged3, violations3] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
    mat_times2_fd1(i/2,j)=toc;

    disp(['Finite differences (classical version) with h=1e-',num2str(i),' : ',flag3]) 
    mat_converged3_fd1(i/2,j)=converged3;
    %conv_ord3(end-10:end) %aggiustare
    mat_val3_fd1(i/2,j)=f3;
    mat_grad3_fd1(i/2,j)=gradf_norm3;
    mat_iter3_fd1(i/2,j)=k3;
    mat_violations3_fd1(i/2,j)=violations3;


    % FINITE DIFFERENCES 2
    JF=@(x) JF_fd2(x,h);
    HF=@(x) HF_fd2(x,h);
    tic;
    [x3, f3, gradf_norm3, k3, xseq3, btseq3,cgiterseq3,conv_ord3,flag3, converged3, violations3] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
    vec_times3_fd2(i,j)=toc;

    disp(['Finite differences (new version) with h=1e-',num2str(i),' : ',flag3]) 
    vec_converged3_fd2(i,j)=converged3;
    %conv_ord3(end-10:end) %aggiustare
    vec_val3_fd2(i,j)=f3;
    vec_grad3_fd2(i,j)=gradf_norm3;
    vec_iter3_fd2(i,j)=k3;
    vec_violations3_fd2(i,j)=violations3;

    end
end

% Forse ha senso plottare poi solo i risultati delle convergenze
% per confrontare i metodi sulle varie dimensioni e enlle varianti ha senso
% usare per esempio la media e le statistiche sui vari successi ottenuti (tipo la media delle iterazioni e del tempo)

% INSERIRE TABELLA
% INSERIRE GRAFICI