addpath('\Users\sofia\Documents\Numerical-Optimization\F79_truncated');
addpath('\Users\sofia\Documents\Numerical-Optimization\F27_truncated');
addpath('\Users\sofia\Documents\Numerical-Optimization\F16_truncated');
addpath('\Users\sofia\Documents\Numerical-Optimization\Nelder');
addpath('\Users\sofia\Documents\Numerical-Optimization\Rosenbrock');


%% FUNCTION 16  (with different initial points)- with exact derivatives and finite differences

sparse=true;

F = @(x) F16(x);  % Defining F16 as function handle
JF_gen = @(x,exact,fin_dif2,h) JF16(x,exact,fin_dif2,h); % Defining JF16 as function handle 
HF_gen= @(x,exact,fin_dif2,h) HF16(x,sparse,exact,fin_dif2,h); % Defining HF16 as function handle (sparse version)

load forcing_terms.mat % possible terms for adaptive tolerance 


%% n=10^5 (1e5)

rng(345989);

n=1e5; 

kmax=1.5e3; % maximum number of iterations of Newton method
tolgrad=5e-4; % tolerance on gradient norm

cg_maxit=100; % maximum number of iterations of coniugate gradient method (for the linear system)
z0=zeros(n,1); % initial point of coniugate gradient method (for the linear system)

% Backtracking parameters
c1=1e-4;
rho=0.50;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

x0 = ones(n, 1);  % initial point
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
mat_conv3_ex=zeros(15:N+1);
vec_converged3_ex=zeros(1,N+1); % vector of booleans (true if it has converged) 
vec_violations3_ex=zeros(1,N+1); % vector with number of violations of curvature condition in Newton method

JF_ex = @(x) JF_gen(x,true,false,0);
HF_ex = @(x) HF_gen(x,true,false,0);

% Structure for derivatives approximated with FINITE DIFFERENCES (classical version) 
mat_times3_fd1=zeros(6,N+1); % matrix with execution times 
mat_val3_fd1=zeros(6,N+1); %matrix with minimal values found
mat_grad3_fd1=zeros(6,N+1); %matrix with final gradient
mat_iter3_fd1=zeros(6,N+1); %matrix with number of iterations 
mat_conv3_fd1=cell(6,N+1);
mat_converged3_fd1=zeros(6,N+1); % matrix of booleans (true if it has converged) 
mat_violations3_fd1=zeros(6,N+1); % matrix with number of violations of curvature condition in Newton method

JF_fd1 = @(x,h) JF_gen(x,false,false,h);
HF_fd1 = @(x,h) HF_gen(x,false,false,h);

% Structure for derivatives approximated with FINITE DIFFERENCES (version with h=h*abs(x_j) as increment) 
mat_times3_fd2=zeros(6,N+1); % matrix with execution times 
mat_val3_fd2=zeros(6,N+1); %matrix with minimal values found
mat_grad3_fd2=zeros(6,N+1); %matrix with final gradient
mat_iter3_fd2=zeros(6,N+1); %matrix with number of iterations 
mat_conv3_fd2=cell(6,N+1);
mat_converged3_fd2=zeros(6,N+1); % matrix of booleans (true if it has converged) 
mat_violations3_fd2=zeros(6,N+1); % matrix with number of violations of curvature condition in Newton method

JF_fd2 = @(x,h) JF_gen(x,false,true,h);
HF_fd2 = @(x,h) HF_gen(x,false,true,h);

for j =1:N+1
    disp(['Condizione iniziale n. ',num2str(j)])

    % EXACT DERIVATIVES
    tic;
    [x3, f3, gradf_norm3, k3, xseq3, btseq3,cgiterseq3,conv_ord3_ex,flag3, converged3, violations3] = truncated_newton(Mat_points(:,j), F, JF_ex, HF_ex, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
    vec_times3_ex(j)=toc;

    disp(['Exact derivatives: ',flag3]) 
    vec_converged3_ex(j)=converged3;
    %conv_ord3(end-10:end) %aggiustare
    vec_val3_ex(j)=f3;
    vec_grad3_ex(j)=gradf_norm3;
    vec_iter3_ex(j)=k3;
    vec_violations3_ex(j)=violations3;
    
    last_vals = conv_ord3_ex(max(end-14,1):end);
    mat_conv3_ex(:, j) = last_vals;
    
    for i=2:2:12
    h=10^(-i);
    
    % FINITE DIFFERENCES 1
    JF=@(x)JF_fd1(x,h);
    HF=@(x)HF_fd1(x,h);
    tic;
    [x3, f3, gradf_norm3, k3, xseq3, btseq3,cgiterseq3,conv_ord3_df1,flag3, converged3, violations3] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
    mat_times2_fd1(i/2,j)=toc;

    disp(['Finite differences (classical version) with h=1e-',num2str(i),' : ',flag3]) 
    mat_converged3_fd1(i/2,j)=converged3;
    %conv_ord3(end-10:end) %aggiustare
    mat_val3_fd1(i/2,j)=f3;
    mat_grad3_fd1(i/2,j)=gradf_norm3;
    mat_iter3_fd1(i/2,j)=k3;
    mat_violations3_fd1(i/2,j)=violations3;
    last_vals = conv_ord3_df1(max(end-14,1):end);
    mat_conv3_fd1(i/2, j) = {last_vals};


    % FINITE DIFFERENCES 2
    JF=@(x) JF_fd2(x,h);
    HF=@(x) HF_fd2(x,h);
    tic;
    [x3, f3, gradf_norm3, k3, xseq3, btseq3,cgiterseq3,conv_ord3_df2,flag3, converged3, violations3] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
    mat_times3_fd2(i/2,j)=toc;

    disp(['Finite differences (new version) with h=1e-',num2str(i),' : ',flag3]) 
    mat_converged3_fd2(i/2,j)=converged3;
    %conv_ord3(end-10:end) %aggiustare
    mat_val3_fd2(i/2,j)=f3;
    mat_grad3_fd2(i/2,j)=gradf_norm3;
    mat_iter3_fd2(i/2,j)=k3;
    mat_violations3_fd2(i/2,j)=violations3;
    last_vals = conv_ord3_df2(max(end-14,1):end);
    mat_conv3_fd2(i/2, j) = {last_vals};


    end
end

% Forse ha senso plottare poi solo i risultati delle convergenze
% per confrontare i metodi sulle varie dimensioni e enlle varianti ha senso
% usare per esempio la media e le statistiche sui vari successi ottenuti (tipo la media delle iterazioni e del tempo)

% INSERIRE TABELLA
% INSERIRE GRAFICI
%%
num_initial_points = N + 1;

% Crea una figura
figure;
hold on;


% Plot per ciascuna condizione iniziale
for j = 1:num_initial_points
    % Estrai l'ordine di convergenza per la j-esima condizione iniziale
    conv_ord_ex = mat_conv3_ex(:,j); % Derivate esatte

    plot(1:15,conv_ord_ex, 'Color', 'b', 'LineWidth', 1.5);
    hold on;
    for i =1:6 
        conv_ord_fd1 = mat_conv3_fd1{i, j}; % Differenze finite classiche

       
        conv_ord_fd2 = mat_conv3_fd2{i, j}; % Differenze finite adattative
        plot(1:15,conv_ord_fd1, '-', 'Color', 'r', 'LineWidth', 1.5);

        hold on;
        plot(1:15,conv_ord_fd2, '-o', 'Color', 'g', 'LineWidth', 1.5);
        hold on;
    end

   
end

% Aggiungi titolo e legenda
title('Ordine di Convergenza per Tutte le Condizioni Iniziali');
xlabel('Iterazione');
ylabel('Ordine di Convergenza');
legend({'Exact Derivatives', 'dif fin_1', 'dif fin_2'}, 'Location', 'Best');
grid on;
hold off;


%% Calcolo delle medie considerando solo le esecuzioni convergenti

% Metodo Exact (derivate esatte) - media unica
vec_times_ex_clean = vec_times3_ex; % copia dei tempi
vec_times_ex_clean(vec_converged3_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_t3 = mean(vec_times_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_times_fd1_clean = mat_times3_fd1;
mat_times_fd1_clean(mat_converged3_fd1 == 0) = NaN;
avg_fd1 = mean(mat_times_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)
%mat_times3_fd2(2:2:end,: )

% Metodo FD2 (Finite Differences nuove)
mat_times_fd2_clean = mat_times3_fd2;
mat_times_fd2_clean(mat_converged3_fd2 == 0) = NaN;
avg_fd2 = mean(mat_times_fd2_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Creazione delle etichette per i valori di h
h_exponents = [2, 4, 6, 8, 10, 12];  % Solo valori di h (senza h=0)
h_labels = arrayfun(@(e) sprintf('h=1e-%d', e), h_exponents, 'UniformOutput', false);

% Preparazione dei dati per la tabella
% FD1 e FD2 hanno le medie per ogni h, mentre Exact è ripetuto in tutte le colonne
%exact_vals = [avg_exact, avg_exact]; % Esatto in tutte le colonne
fd1_vals = avg_fd1';  % Trasposto per allineare con le colonne
fd2_vals = avg_fd2';  % Trasposto per allineare con le colonne

% Costruzione della tabella
rowNames = {'FD1', 'FD2'};
columnNames = [ h_labels,'Exact']; % Prima colonna "Exact", poi gli h
data = [ fd1_vals, avg_exact_t3; fd2_vals, avg_exact_t3;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T7 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation times table (only for successful runs): F16, n=10^5, superlinear');
disp(T7);

% (Opzionale) Salva la tabella in un file CSV
%writetable(T, 'Time_tabel_f16_3_suplin.csv', 'WriteRowNames', true);


%% Calcolo delle medie considerando solo le esecuzioni convergenti

% Metodo Exact (derivate esatte) - media unica
vec_times_ex_clean = vec_iter3_ex; % copia dei tempi
vec_times_ex_clean(vec_converged3_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_i3 = mean(vec_times_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_times_fd1_clean = mat_iter3_fd1;
mat_times_fd1_clean(mat_converged3_fd1 == 0) = NaN;
avg_fd1 = mean(mat_times_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_times_fd2_clean = mat_iter3_fd2;
mat_times_fd2_clean(mat_converged3_fd2 == 0) = NaN;
avg_fd2 = mean(mat_times_fd2_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Creazione delle etichette per i valori di h
h_exponents = [2, 4, 6, 8, 10, 12];  % Solo valori di h (senza h=0)
h_labels = arrayfun(@(e) sprintf('h=1e-%d', e), h_exponents, 'UniformOutput', false);

% Preparazione dei dati per la tabella
% FD1 e FD2 hanno le medie per ogni h, mentre Exact è ripetuto in tutte le colonne
%exact_vals = [avg_exact, avg_exact]; % Esatto in tutte le colonne
fd1_vals = avg_fd1';  % Trasposto per allineare con le colonne
fd2_vals = avg_fd2';  % Trasposto per allineare con le colonne

% Costruzione della tabella
rowNames = {'FD1', 'FD2'};
columnNames = [ h_labels,'Exact']; % Prima colonna "Exact", poi gli h
data = [ fd1_vals, avg_exact_i3; fd2_vals, avg_exact_i3;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T8 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation iteration table (only for successful runs): F16, n=10^5, superlinear');
disp(T8);

% (Opzionale) Salva la tabella in un file CSV
%writetable(T, 'Iteration_tabel_f16_3_suplin.csv', 'WriteRowNames', true);

%% Calcolo delle medie considerando solo le esecuzioni convergenti

% Metodo Exact (derivate esatte) - media unica
vec_times_ex_clean = vec_val3_ex; % copia dei tempi
vec_times_ex_clean(vec_converged3_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_f3 = mean(vec_times_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_times_fd1_clean = mat_val3_fd1;
mat_times_fd1_clean(mat_converged3_fd1 == 0) = NaN;
avg_fd1 = mean(mat_times_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_times_fd2_clean = mat_val3_fd2;
mat_times_fd2_clean(mat_converged3_fd2 == 0) = NaN;
avg_fd2 = mean(mat_times_fd2_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Creazione delle etichette per i valori di h
h_exponents = [2, 4, 6, 8, 10, 12];  % Solo valori di h (senza h=0)
h_labels = arrayfun(@(e) sprintf('h=1e-%d', e), h_exponents, 'UniformOutput', false);

% Preparazione dei dati per la tabella
% FD1 e FD2 hanno le medie per ogni h, mentre Exact è ripetuto in tutte le colonne
%exact_vals = [avg_exact, avg_exact]; % Esatto in tutte le colonne
fd1_vals = avg_fd1';  % Trasposto per allineare con le colonne
fd2_vals = avg_fd2';  % Trasposto per allineare con le colonne

% Costruzione della tabella
rowNames = {'FD1', 'FD2'};
columnNames = [ h_labels,'Exact']; % Prima colonna "Exact", poi gli h
data = [ fd1_vals, avg_exact_f3; fd2_vals, avg_exact_f3;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T9 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation fmin value table (only for successful runs): F16, n=10^5, superlinear');
disp(T9);

% (Opzionale) Salva la tabella in un file CSV
%writetable(T, 'Fminvalue_tabel_f16_3_suplin.csv', 'WriteRowNames', true);

writetable(T7, 'results_f16_suplin.xlsx', 'Sheet', 'time5','WriteRowNames', true);
writetable(T8, 'results_f16_suplin.xlsx', 'Sheet', 'niter_5','WriteRowNames', true);
writetable(T9, 'results_f16_suplin.xlsx', 'Sheet', 'f_val_5','WriteRowNames', true);


% %%
% % Creazione della tabella
% data = [avg_exact_t1, avg_exact_t2, avg_exact_t3;
%         avg_exact_i1, avg_exact_i2, avg_exact_i3;
%         avg_exact_f1, avg_exact_f2, avg_exact_f3];
% 
% % Definizione delle intestazioni
% rowNames = {'Avergae Time', 'Avergae Iter', 'Average fval'};
% columnNames = {'n=10^3', 'n=10^4', 'n=10^5'};
% 
% 
% 
% % Creazione tabella MATLAB
% T_compare = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);
% disp(T_compare)
% 
% % Salvataggio su Excel
% writetable(T_compare, 'results_f16_suplin.xlsx', 'Sheet', 'ExactComparison', 'WriteRowNames', true);
