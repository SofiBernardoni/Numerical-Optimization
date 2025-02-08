%% FUNCTION 16  (with different initial points)- with exact derivatives and finite differences

sparse=true;

F = @(x) F16(x);  % Defining F16 as function handle
JF_gen = @(x,exact,fin_dif2,h) JF16(x,exact,fin_dif2,h); % Defining JF16 as function handle 
HF_gen= @(x,exact,fin_dif2,h) HF16(x,sparse,exact,fin_dif2,h); % Defining HF16 as function handle (sparse version)

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

x0 = ones(n, 1); % initial point
N=10; % number of initial points to be generated

% Initial poin7ts:
Mat_points=repmat(x0,1,N+1); 
rand_mat=2*rand(n, N)-1;
Mat_points(:,2:end)=Mat_points(:,2:end) + rand_mat; % matrix with columns=initial points

% Structure for EXACT derivatives
vec_times1_ex=zeros(1,N+1); % vector with execution times 
vec_val1_ex=zeros(1,N+1); %vector with minimal values found
vec_grad1_ex=zeros(1,N+1); %vector with final gradient
vec_iter1_ex=zeros(1,N+1); %vector with number of iterations 
vec_cg_iter1_ex=zeros(1,N+1); %vector with mean number of inner iterations
vec_bt1_ex=zeros(1,N+1); %vector with mean number of backtracking iterations
mat_conv1_ex=zeros(12, N+1);
vec_converged1_ex=zeros(1,N+1); % vector of booleans (true if it has converged) 
vec_violations1_ex=zeros(1,N+1); % vector with number of violations of curvature condition in Newton method

JF_ex = @(x) JF_gen(x,true,false,0);
HF_ex = @(x) HF_gen(x,true,false,0);

% Structure for derivatives approximated with FINITE DIFFERENCES (classical version) 
mat_times1_fd1=zeros(6,N+1); % matrix with execution times 
mat_val1_fd1=zeros(6,N+1); %matrix with minimal values found
mat_grad1_fd1=zeros(6,N+1); %matrix with final gradient
mat_iter1_fd1=zeros(6,N+1); %matrix with number of iterations 
mat_cg_iter1_fd1=zeros(6,N+1); %matrix with mean number of inner iterations
mat_bt1_fd1=zeros(6,N+1); %matrix with mean number of backtracking iterations
mat_conv1_fd1=cell(6, N+1);
mat_converged1_fd1=zeros(6,N+1); % matrix of booleans (true if it has converged) 
mat_violations1_fd1=zeros(6,N+1); % matrix with number of violations of curvature condition in Newton method

JF_fd1 = @(x,h) JF_gen(x,false,false,h);
HF_fd1 = @(x,h) HF_gen(x,false,false,h);

% Structure for derivatives approximated with FINITE DIFFERENCES (version with h=h*abs(x_j) as increment) 
mat_times1_fd2=zeros(6,N+1); % matrix with execution times 
mat_val1_fd2=zeros(6,N+1); %matrix with minimal values found
mat_grad1_fd2=zeros(6,N+1); %matrix with final gradient
mat_iter1_fd2=zeros(6,N+1); %matrix with number of iterations 
mat_cg_iter1_fd2=zeros(6,N+1); %matrix with mean number of inner iterations
mat_bt1_fd2=zeros(6,N+1); %matrix with mean number of backtracking iterations
mat_conv1_fd2=cell(6,N+1);
mat_converged1_fd2=zeros(6,N+1); % matrix of booleans (true if it has converged) 
mat_violations1_fd2=zeros(6,N+1); % matrix with number of violations of curvature condition in Newton method

JF_fd2 = @(x,h) JF_gen(x,false,true,h);
HF_fd2 = @(x,h) HF_gen(x,false,true,h);

for j =1:N+1
    disp(['Condizione iniziale n. ',num2str(j)])

    % EXACT DERIVATIVES
    tic;
    [x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1_ex,flag1, converged1, violations1] = truncated_newton(Mat_points(:,j), F, JF_ex, HF_ex, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
    vec_times1_ex(j)=toc;

    disp(['Exact derivatives: ',flag1]) 
    vec_converged1_ex(j)=converged1;
    vec_val1_ex(j)=f1;
    vec_grad1_ex(j)=gradf_norm1;
    vec_iter1_ex(j)=k1;
    vec_cg_iter1_ex(j)=sum(cgiterseq1)/k1; 
    vec_bt1_ex(j)=sum(btseq1)/k1; 
    vec_violations1_ex(j)=violations1;
    last_vals = conv_ord1_ex(max(end-11,1):end);
    mat_conv1_ex(:, j) = last_vals;

    
    for i=2:2:12
    h=10^(-i);
    
    % FINITE DIFFERENCES 1
    JF=@(x)JF_fd1(x,h);
    HF=@(x)HF_fd1(x,h);
    tic;
    [x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1_df1,flag1, converged1, violations1] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
    mat_times1_fd1(i/2,j)=toc;

    disp(['Finite differences (classical version) with h=1e-',num2str(i),' : ',flag1]) 
    mat_converged1_fd1(i/2,j)=converged1;
    mat_val1_fd1(i/2,j)=f1;
    mat_grad1_fd1(i/2,j)=gradf_norm1;
    mat_iter1_fd1(i/2,j)=k1;
    mat_cg_iter1_fd1(i/2,j)=sum(cgiterseq1)/k1; 
    mat_bt1_fd1(i/2,j)=sum(btseq1)/k1;
    mat_violations1_fd1(i/2,j)=violations1;
    last_vals = conv_ord1_df1(max(end-11,1):end);
    mat_conv1_fd1(i/2, j) = {last_vals};


    % FINITE DIFFERENCES 2
    JF=@(x) JF_fd2(x,h);
    HF=@(x) HF_fd2(x,h);
    tic;
    [x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1_df2,flag1, converged1, violations1] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
    mat_times1_fd2(i/2,j)=toc;

    disp(['Finite differences (new version) with h=1e-',num2str(i),' : ',flag1]) 
    mat_converged1_fd2(i/2,j)=converged1;
    mat_val1_fd2(i/2,j)=f1;
    mat_grad1_fd2(i/2,j)=gradf_norm1;
    mat_iter1_fd2(i/2,j)=k1;
    mat_cg_iter1_fd2(i/2,j)=sum(cgiterseq1)/k1; 
    mat_bt1_fd2(i/2,j)=sum(btseq1)/k1;
    mat_violations1_fd2(i/2,j)=violations1;
    last_vals = conv_ord1_df2(max(end-11,1):end);
    mat_conv1_fd2(i/2, j) = {last_vals};

    end
end



%% GRAFICO CON INNER ITERATION
num_initial_points = N + 1;

% Crea una figura
figure;
hold on;


% Plot per ciascuna condizione iniziale
for j = 1:num_initial_points
    % Estrai l'ordine di convergenza per la j-esima condizione iniziale
    conv_ord_ex = mat_conv1_ex(:,j); % Derivate esatte

    plot(1:12,conv_ord_ex, 'Color', 'b', 'LineWidth', 1.5);
    hold on;
    for i =1:6 
        conv_ord_fd1 = mat_conv1_fd1{i, j}; % Differenze finite classiche

       
        conv_ord_fd2 = mat_conv1_fd2{i, j}; % Differenze finite adattative
        plot(1:12,conv_ord_fd1, '-', 'Color', 'r', 'LineWidth', 1.5);

        hold on;
        plot(1:12,conv_ord_fd2, '-o', 'Color', 'g', 'LineWidth', 1.5);
        hold on;
    end

   
end

% Aggiungi titolo e legenda
title('F16, n=10^3, superlinear');
xlabel('Iterazione');
ylabel('Ordine di Convergenza');
legend({'Exact Derivatives', 'dif fin_1', 'dif fin_2'}, 'Location', 'Best');
grid on;
hold off;


%% Calcolo delle medie considerando solo le esecuzioni convergenti

% Metodo Exact (derivate esatte) - media unica
vec_times_ex_clean = vec_times1_ex; % copia dei tempi
vec_times_ex_clean(vec_converged1_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_t1 = mean(vec_times_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_times_fd1_clean = mat_times1_fd1;
mat_times_fd1_clean(mat_converged1_fd1 == 0) = NaN;
avg_fd1 = mean(mat_times_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_times_fd2_clean = mat_times1_fd2;
mat_times_fd2_clean(mat_converged1_fd2 == 0) = NaN;
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
data = [ fd1_vals, avg_exact_t1; fd2_vals, avg_exact_t1;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T1 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation times table (only for successful runs): F16, n=10^3, superlinear');
disp(T1);

% (Opzionale) Salva la tabella in un file CSV
%writetable(T, 'Time_tabel_f16_3_suplin.csv', 'WriteRowNames', true);


%% Calcolo delle medie considerando solo le esecuzioni convergenti

% Metodo Exact (derivate esatte) - media unica
vec_times_ex_clean = vec_iter1_ex; % copia dei tempi
vec_times_ex_clean(vec_converged1_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_i1 = mean(vec_times_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_times_fd1_clean = mat_iter1_fd1;
mat_times_fd1_clean(mat_converged1_fd1 == 0) = NaN;
avg_fd1 = mean(mat_times_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_times_fd2_clean = mat_iter1_fd2;
mat_times_fd2_clean(mat_converged1_fd2 == 0) = NaN;
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
data = [ fd1_vals, avg_exact_i1; fd2_vals, avg_exact_i1;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T2 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation iteration table (only for successful runs): F16, n=10^3, superlinear');
disp(T2);

% (Opzionale) Salva la tabella in un file CSV
%writetable(T, 'Iteration_tabel_f16_3_suplin.csv', 'WriteRowNames', true);





%% Calcolo delle medie considerando solo le esecuzioni convergenti

% Metodo Exact (derivate esatte) - media unica
vec_times_ex_clean = vec_val1_ex; % copia dei tempi
vec_times_ex_clean(vec_converged1_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_f1 = mean(vec_times_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_times_fd1_clean = mat_val1_fd1;
mat_times_fd1_clean(mat_converged1_fd1 == 0) = NaN;
avg_fd1 = mean(mat_times_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_times_fd2_clean = mat_val1_fd2;
mat_times_fd2_clean(mat_converged1_fd2 == 0) = NaN;
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
data = [ fd1_vals, avg_exact_f1; fd2_vals, avg_exact_f1;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T3 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation fmin value table (only for successful runs): F16, n=10^3, superlinear');
disp(T3);




%% VIOLATION

% Metodo Exact (derivate esatte) - media unica
vec_times_ex_clean = vec_violations1_ex; % copia dei tempi
vec_times_ex_clean(vec_converged1_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_v1 = mean(vec_times_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_times_fd1_clean = mat_violations1_fd1;
mat_times_fd1_clean(mat_converged1_fd1 == 0) = NaN;
avg_fd1 = mean(mat_times_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_times_fd2_clean = mat_violations1_fd2;
mat_times_fd2_clean(mat_converged1_fd2 == 0) = NaN;
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
data = [ fd1_vals, avg_exact_v1; fd2_vals, avg_exact_v1;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T10 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation violation  table (only for successful runs): F16, n=10^3, superlinear');
disp(T10);


%% BT-SEQ
% Metodo Exact (derivate esatte) - media unica
vec_bt_ex_clean = vec_bt1_ex; % copia dei tempi
vec_bt_ex_clean(vec_converged1_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_bt1 = mean(vec_bt_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_bt_fd1_clean = mat_bt1_fd1;
mat_bt_fd1_clean(mat_converged1_fd1 == 0) = NaN;
avg_fd1 = mean(mat_bt_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_bt_fd2_clean = mat_bt1_fd2;
mat_bt_fd2_clean(mat_converged1_fd2 == 0) = NaN;
avg_fd2 = mean(mat_bt_fd2_clean, 2, 'omitnan'); % media per ogni h (6x1)

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
data = [ fd1_vals, avg_exact_bt1; fd2_vals, avg_exact_bt1;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T11 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation bt iteration table (only for successful runs): F16, n=10^3, superlinear');
disp(T11);



%% CG-SEQ
% Metodo Exact (derivate esatte) - media unica
vec_bt_ex_clean = vec_cg_iter1_ex; % copia dei tempi
vec_bt_ex_clean(vec_converged1_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_cg1 = mean(vec_bt_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_bt_fd1_clean = mat_cg_iter1_fd1;
mat_bt_fd1_clean(mat_converged1_fd1 == 0) = NaN;
avg_fd1 = mean(mat_bt_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_bt_fd2_clean = mat_cg_iter1_fd2;
mat_bt_fd2_clean(mat_converged1_fd2 == 0) = NaN;
avg_fd2 = mean(mat_bt_fd2_clean, 2, 'omitnan'); % media per ogni h (6x1)

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
data = [ fd1_vals, avg_exact_cg1; fd2_vals, avg_exact_cg1;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T12 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation cg iteration table (only for successful runs): F16, n=10^3, superlinear');
disp(T12);



%% Calcolo quanti a convergenza

% Creazione delle etichette per i valori di h
h_exponents = [2, 4, 6, 8, 10, 12];  % Solo valori di h (senza h=0)
h_labels = arrayfun(@(e) sprintf('h=1e-%d', e), h_exponents, 'UniformOutput', false);

% Preparazione dei dati per la tabella
% FD1 e FD2 hanno le medie per ogni h, mentre Exact è ripetuto in tutte le colonne
%exact_vals = [avg_exact, avg_exact]; % Esatto in tutte le colonne
fd1_vals = sum(mat_converged1_fd1,2)';  % Trasposto per allineare con le colonne
fd2_vals = sum(mat_converged1_fd2,2)';  % Trasposto per allineare con le colonne

% Costruzione della tabella
rowNames = {'FD1', 'FD2'};
columnNames = [ h_labels,'Exact']; % Prima colonna "Exact", poi gli h
data = [ fd1_vals, sum(vec_converged1_ex); fd2_vals, sum(vec_converged1_ex);];

% Creiamo la tabella con i nomi delle colonne e delle righe
T13 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Number of converged : F16, n=10^3, superlinear');
disp(T13);

writetable(T1, 'results_f16_suplin.xlsx', 'Sheet', 'time_3','WriteRowNames', true);
writetable(T2, 'results_f16_suplin.xlsx', 'Sheet', 'niter_3','WriteRowNames', true);
writetable(T3, 'results_f16_suplin.xlsx', 'Sheet', 'f_val_3','WriteRowNames', true);
writetable(T10, 'results_f16_suplin.xlsx', 'Sheet', 'viol_3','WriteRowNames', true);
writetable(T11, 'results_f16_suplin.xlsx', 'Sheet', 'bt_3','WriteRowNames', true);
writetable(T12, 'results_f16_suplin.xlsx', 'Sheet', 'cg_3','WriteRowNames', true);
writetable(T13, 'results_f16_suplin.xlsx', 'Sheet', 'n_conv3','WriteRowNames', true);



%% n=10^4 (1e4)

rng(345989);

n=1e4; 

kmax=1.5e3; % maximum number of iterations of Newton method
tolgrad=1e-5; % tolerance on gradient norm %%%%%%%%%%%%%%%%%%%% decide if we want to keep the tolerance 5e-7

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
vec_times2_ex=zeros(1,N+1); % vector with execution times 
vec_val2_ex=zeros(1,N+1); %vector with minimal values found
vec_grad2_ex=zeros(1,N+1); %vector with final gradient
vec_iter2_ex=zeros(1,N+1); %vector with number of iterations 
vec_cg_iter2_ex=zeros(1,N+1); %vector with mean number of inner iterations
vec_bt2_ex=zeros(1,N+1); %vector with mean number of backtracking iterations
mat_conv2_ex=zeros(12,N+1);
vec_converged2_ex=zeros(1,N+1); % vector of booleans (true if it has converged) 
vec_violations2_ex=zeros(1,N+1); % vector with number of violations of curvature condition in Newton method

JF_ex = @(x) JF_gen(x,true,false,0);
HF_ex = @(x) HF_gen(x,true,false,0);

% Structure for derivatives approximated with FINITE DIFFERENCES (classical version) 
mat_times2_fd1=zeros(6,N+1); % matrix with execution times 
mat_val2_fd1=zeros(6,N+1); %matrix with minimal values found
mat_grad2_fd1=zeros(6,N+1); %matrix with final gradient
mat_iter2_fd1=zeros(6,N+1); %matrix with number of iterations 
mat_cg_iter2_fd1=zeros(6,N+1); %matrix with mean number of inner iterations
mat_bt2_fd1=zeros(6,N+1); %matrix with mean number of backtracking iterations
mat_conv2_fd1=cell(6,N+1);
mat_converged2_fd1=zeros(6,N+1); % matrix of booleans (true if it has converged) 
mat_violations2_fd1=zeros(6,N+1); % matrix with number of violations of curvature condition in Newton method

JF_fd1 = @(x,h) JF_gen(x,false,false,h);
HF_fd1 = @(x,h) HF_gen(x,false,false,h);

% Structure for derivatives approximated with FINITE DIFFERENCES (version with h=h*abs(x_j) as increment) 
mat_times2_fd2=zeros(6,N+1); % matrix with execution times 
mat_val2_fd2=zeros(6,N+1); %matrix with minimal values found
mat_grad2_fd2=zeros(6,N+1); %matrix with final gradient
mat_iter2_fd2=zeros(6,N+1); %matrix with number of iterations 
mat_cg_iter2_fd2=zeros(6,N+1); %matrix with mean number of inner iterations
mat_bt2_fd2=zeros(6,N+1); %matrix with mean number of backtracking iterations
mat_conv2_fd2=cell(6,N+1);
mat_converged2_fd2=zeros(6,N+1); % matrix of booleans (true if it has converged) 
mat_violations2_fd2=zeros(6,N+1); % matrix with number of violations of curvature condition in Newton method

JF_fd2 = @(x,h) JF_gen(x,false,true,h);
HF_fd2 = @(x,h) HF_gen(x,false,true,h);

for j =1:N+1
    disp(['Condizione iniziale n. ',num2str(j)])

    % EXACT DERIVATIVES
    tic;
    [x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2_ex,flag2, converged2, violations2] = truncated_newton(Mat_points(:,j), F, JF_ex, HF_ex, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
    vec_times2_ex(j)=toc;

    disp(['Exact derivatives: ',flag2]) 
    vec_converged2_ex(j)=converged2;
    vec_val2_ex(j)=f2;
    vec_grad2_ex(j)=gradf_norm2;
    vec_iter2_ex(j)=k2;
    vec_cg_iter2_ex(j)=sum(cgiterseq2)/k2; 
    vec_bt2_ex(j)=sum(btseq2)/k2; 
    vec_violations2_ex(j)=violations2;
    last_vals = conv_ord2_ex(max(end-11,1):end);
    mat_conv2_ex(:, j) = last_vals;
    
    for i=2:2:12
    h=10^(-i);
    
    % FINITE DIFFERENCES 1
    JF=@(x)JF_fd1(x,h);
    HF=@(x)HF_fd1(x,h);
    tic;
    [x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2_df1,flag2, converged2, violations2] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
    mat_times2_fd1(i/2,j)=toc;

    disp(['Finite differences (classical version) with h=1e-',num2str(i),' : ',flag2]) 
    mat_converged2_fd1(i/2,j)=converged2;
    mat_val2_fd1(i/2,j)=f2;
    mat_grad2_fd1(i/2,j)=gradf_norm2;
    mat_iter2_fd1(i/2,j)=k2;
    mat_cg_iter2_fd1(i/2,j)=sum(cgiterseq2)/k2; 
    mat_bt2_fd1(i/2,j)=sum(btseq2)/k2;
    mat_violations2_fd1(i/2,j)=violations2;
    last_vals = conv_ord2_df1(max(end-11,1):end);
    mat_conv2_fd1(i/2, j) = {last_vals};


    % FINITE DIFFERENCES 2
    JF=@(x) JF_fd2(x,h);
    HF=@(x) HF_fd2(x,h);
    tic;
    [x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2_df2,flag2, converged2, violations2] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
    mat_times2_fd2(i/2,j)=toc;

    disp(['Finite differences (new version) with h=1e-',num2str(i),' : ',flag2]) 
    mat_converged2_fd2(i/2,j)=converged2;
    mat_val2_fd2(i/2,j)=f2;
    mat_grad2_fd2(i/2,j)=gradf_norm2;
    mat_iter2_fd2(i/2,j)=k2;
    mat_cg_iter2_fd2(i/2,j)=sum(cgiterseq2)/k2; 
    mat_bt2_fd2(i/2,j)=sum(btseq2)/k2;
    mat_violations2_fd2(i/2,j)=violations2;
    last_vals = conv_ord2_df2(max(end-11,1):end);
    mat_conv2_fd2(i/2, j) = {last_vals};

    end
end

%%
num_initial_points = N + 1;

% Crea una figura
figure;
hold on;


% Plot per ciascuna condizione iniziale
for j = 1:num_initial_points
    % Estrai l'ordine di convergenza per la j-esima condizione iniziale
    conv_ord_ex = mat_conv2_ex(:,j); % Derivate esatte

    plot(1:12,conv_ord_ex, 'Color', 'b', 'LineWidth', 1.5);
    hold on;
    for i =1:6 
        conv_ord_fd1 = mat_conv2_fd1{i, j}; % Differenze finite classiche

       
        conv_ord_fd2 = mat_conv2_fd2{i, j}; % Differenze finite adattative
        plot(1:12,conv_ord_fd1, '-', 'Color', 'r', 'LineWidth', 1.5);

        hold on;
        plot(1:12,conv_ord_fd2, '-o', 'Color', 'g', 'LineWidth', 1.5);
        hold on;
    end

   
end

% Aggiungi titolo e legenda
title('F16, n=10^4, superlinear');
xlabel('Iterazione');
ylabel('Ordine di Convergenza');
legend({'Exact Derivatives', 'dif fin_1', 'dif fin_2'}, 'Location', 'Best');
grid on;
hold off;


%% Calcolo delle medie considerando solo le esecuzioni convergenti

% Metodo Exact (derivate esatte) - media unica
vec_times_ex_clean = vec_times2_ex; % copia dei tempi
vec_times_ex_clean(vec_converged2_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_t2 = mean(vec_times_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_times_fd1_clean = mat_times2_fd1;
mat_times_fd1_clean(mat_converged2_fd1 == 0) = NaN;
avg_fd1 = mean(mat_times_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_times_fd2_clean = mat_times2_fd2;
mat_times_fd2_clean(mat_converged2_fd2 == 0) = NaN;
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
data = [ fd1_vals, avg_exact_t2; fd2_vals, avg_exact_t2;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T4 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation times table (only for successful runs): F16, n=10^4, superlinear');
disp(T4);

% (Opzionale) Salva la tabella in un file CSV
%writetable(T, 'Time_tabel_f16_3_suplin.csv', 'WriteRowNames', true);


%% Calcolo delle medie considerando solo le esecuzioni convergenti

% Metodo Exact (derivate esatte) - media unica
vec_times_ex_clean = vec_iter2_ex; % copia dei tempi
vec_times_ex_clean(vec_converged2_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_i2 = mean(vec_times_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_times_fd1_clean = mat_iter2_fd1;
mat_times_fd1_clean(mat_converged2_fd1 == 0) = NaN;
avg_fd1 = mean(mat_times_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_times_fd2_clean = mat_iter2_fd2;
mat_times_fd2_clean(mat_converged2_fd2 == 0) = NaN;
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
data = [ fd1_vals, avg_exact_i2; fd2_vals, avg_exact_i2;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T5 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation iteration table (only for successful runs): F16, n=10^4, superlinear');
disp(T5);

% (Opzionale) Salva la tabella in un file CSV
%writetable(T, 'Iteration_tabel_f16_3_suplin.csv', 'WriteRowNames', true);

%% Calcolo delle medie considerando solo le esecuzioni convergenti

% Metodo Exact (derivate esatte) - media unica
vec_times_ex_clean = vec_val2_ex; % copia dei tempi
vec_times_ex_clean(vec_converged2_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_f2 = mean(vec_times_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_times_fd1_clean = mat_val2_fd1;
mat_times_fd1_clean(mat_converged2_fd1 == 0) = NaN;
avg_fd1 = mean(mat_times_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_times_fd2_clean = mat_val2_fd2;
mat_times_fd2_clean(mat_converged2_fd2 == 0) = NaN;
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
data = [ fd1_vals, avg_exact_f2; fd2_vals, avg_exact_f2;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T6 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation fmin value table (only for successful runs): F16, n=10^4, superlinear');
disp(T6);

% (Opzionale) Salva la tabella in un file CSV
%writetable(T, 'Fminvalue_tabel_f16_3_suplin.csv', 'WriteRowNames', true);

%% VIOLATION

% Metodo Exact (derivate esatte) - media unica
vec_times_ex_clean = vec_violations2_ex; % copia dei tempi
vec_times_ex_clean(vec_converged2_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_v2 = mean(vec_times_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_times_fd1_clean = mat_violations2_fd1;
mat_times_fd1_clean(mat_converged2_fd1 == 0) = NaN;
avg_fd1 = mean(mat_times_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_times_fd2_clean = mat_violations2_fd2;
mat_times_fd2_clean(mat_converged2_fd2 == 0) = NaN;
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
data = [ fd1_vals, avg_exact_v2; fd2_vals, avg_exact_v2;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T14 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation violation  table (only for successful runs): F16, n=10^4, suplin');
disp(T14);


%% BT-SEQ
% Metodo Exact (derivate esatte) - media unica
vec_bt_ex_clean = vec_bt2_ex; % copia dei tempi
vec_bt_ex_clean(vec_converged2_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_bt2 = mean(vec_bt_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_bt_fd1_clean = mat_bt2_fd1;
mat_bt_fd1_clean(mat_converged2_fd1 == 0) = NaN;
avg_fd1 = mean(mat_bt_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_bt_fd2_clean = mat_bt2_fd2;
mat_bt_fd2_clean(mat_converged2_fd2 == 0) = NaN;
avg_fd2 = mean(mat_bt_fd2_clean, 2, 'omitnan'); % media per ogni h (6x1)

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
data = [ fd1_vals, avg_exact_bt2; fd2_vals, avg_exact_bt2;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T15 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation bt iteration table (only for successful runs): F16, n=10^4, suplin');
disp(T15);



%% CG-SEQ
% Metodo Exact (derivate esatte) - media unica
vec_bt_ex_clean = vec_cg_iter2_ex; % copia dei tempi
vec_bt_ex_clean(vec_converged2_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_cg2 = mean(vec_bt_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_bt_fd1_clean = mat_cg_iter2_fd1;
mat_bt_fd1_clean(mat_converged2_fd1 == 0) = NaN;
avg_fd1 = mean(mat_bt_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_bt_fd2_clean = mat_cg_iter2_fd2;
mat_bt_fd2_clean(mat_converged2_fd2 == 0) = NaN;
avg_fd2 = mean(mat_bt_fd2_clean, 2, 'omitnan'); % media per ogni h (6x1)

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
data = [ fd1_vals, avg_exact_cg2; fd2_vals, avg_exact_cg2;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T16 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation cg iteration table (only for successful runs): F16, n=10^4, suplin');
disp(T16);



%% Calcolo quanti a convergenza


% Creazione delle etichette per i valori di h
h_exponents = [2, 4, 6, 8, 10, 12];  % Solo valori di h (senza h=0)
h_labels = arrayfun(@(e) sprintf('h=1e-%d', e), h_exponents, 'UniformOutput', false);

% Preparazione dei dati per la tabella
% FD1 e FD2 hanno le medie per ogni h, mentre Exact è ripetuto in tutte le colonne
%exact_vals = [avg_exact, avg_exact]; % Esatto in tutte le colonne
fd1_vals = sum(mat_converged2_fd1,2)';  % Trasposto per allineare con le colonne
fd2_vals = sum(mat_converged2_fd2,2)';  % Trasposto per allineare con le colonne

% Costruzione della tabella
rowNames = {'FD1', 'FD2'};
columnNames = [ h_labels,'Exact']; % Prima colonna "Exact", poi gli h
data = [ fd1_vals, sum(vec_converged2_ex); fd2_vals, sum(vec_converged2_ex);];

% Creiamo la tabella con i nomi delle colonne e delle righe
T17 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Number of converged : F16, n=10^4, suplin');
disp(T17);

writetable(T4, 'results_f16_suplin.xlsx', 'Sheet', 'time_4','WriteRowNames', true);
writetable(T5, 'results_f16_suplin.xlsx', 'Sheet', 'niter_4','WriteRowNames', true);
writetable(T6, 'results_f16_suplin.xlsx', 'Sheet', 'f_val_4','WriteRowNames', true);
writetable(T14, 'results_f16_suplin.xlsx', 'Sheet', 'viol_4','WriteRowNames', true);
writetable(T15, 'results_f16_suplin.xlsx', 'Sheet', 'bt_4','WriteRowNames', true);
writetable(T16, 'results_f16_suplin.xlsx', 'Sheet', 'cg_4','WriteRowNames', true);
writetable(T17, 'results_f16_suplin.xlsx', 'Sheet', 'n_conv4','WriteRowNames', true);


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
vec_cg_iter3_ex=zeros(1,N+1); %vector with mean number of inner iterations
vec_bt3_ex=zeros(1,N+1); %vector with mean number of backtracking iterations
mat_conv3_ex=zeros(12:N+1);
vec_converged3_ex=zeros(1,N+1); % vector of booleans (true if it has converged) 
vec_violations3_ex=zeros(1,N+1); % vector with number of violations of curvature condition in Newton method

JF_ex = @(x) JF_gen(x,true,false,0);
HF_ex = @(x) HF_gen(x,true,false,0);

% Structure for derivatives approximated with FINITE DIFFERENCES (classical version) 
mat_times3_fd1=zeros(6,N+1); % matrix with execution times 
mat_val3_fd1=zeros(6,N+1); %matrix with minimal values found
mat_grad3_fd1=zeros(6,N+1); %matrix with final gradient
mat_iter3_fd1=zeros(6,N+1); %matrix with number of iterations 
mat_cg_iter3_fd1=zeros(6,N+1); %matrix with mean number of inner iterations
mat_bt3_fd1=zeros(6,N+1); %matrix with mean number of backtracking iterations
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
mat_cg_iter3_fd2=zeros(6,N+1); %matrix with mean number of inner iterations
mat_bt3_fd2=zeros(6,N+1); %matrix with mean number of backtracking iterations
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
    vec_val3_ex(j)=f3;
    vec_grad3_ex(j)=gradf_norm3;
    vec_iter3_ex(j)=k3;
    vec_cg_iter3_ex(j)=sum(cgiterseq3)/k3; 
    vec_bt3_ex(j)=sum(btseq3)/k3; 
    vec_violations3_ex(j)=violations3;
    last_vals = conv_ord3_ex(max(end-11,1):end);
    mat_conv3_ex(:, j) = last_vals;
    
    for i=2:2:12
    h=10^(-i);
    
    % FINITE DIFFERENCES 1
    JF=@(x)JF_fd1(x,h);
    HF=@(x)HF_fd1(x,h);
    tic;
    [x3, f3, gradf_norm3, k3, xseq3, btseq3,cgiterseq3,conv_ord3_df1,flag3, converged3, violations3] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
    mat_times3_fd1(i/2,j)=toc;

    disp(['Finite differences (classical version) with h=1e-',num2str(i),' : ',flag3]) 
    mat_converged3_fd1(i/2,j)=converged3;
    mat_val3_fd1(i/2,j)=f3;
    mat_grad3_fd1(i/2,j)=gradf_norm3;
    mat_iter3_fd1(i/2,j)=k3;
    mat_cg_iter3_fd1(i/2,j)=sum(cgiterseq3)/k3; 
    mat_bt3_fd1(i/2,j)=sum(btseq3)/k3;
    mat_violations3_fd1(i/2,j)=violations3;
    last_vals = conv_ord3_df1(max(end-11,1):end);
    mat_conv3_fd1(i/2, j) = {last_vals};


    % FINITE DIFFERENCES 2
    JF=@(x) JF_fd2(x,h);
    HF=@(x) HF_fd2(x,h);
    tic;
    [x3, f3, gradf_norm3, k3, xseq3, btseq3,cgiterseq3,conv_ord3_df2,flag3, converged3, violations3] = truncated_newton(Mat_points(:,j), F, JF, HF, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
    mat_times3_fd2(i/2,j)=toc;

    disp(['Finite differences (new version) with h=1e-',num2str(i),' : ',flag3]) 
    mat_converged3_fd2(i/2,j)=converged3;
    mat_val3_fd2(i/2,j)=f3;
    mat_grad3_fd2(i/2,j)=gradf_norm3;
    mat_iter3_fd2(i/2,j)=k3;
    mat_cg_iter3_fd2(i/2,j)=sum(cgiterseq3)/k3; 
    mat_bt3_fd2(i/2,j)=sum(btseq3)/k3;
    mat_violations3_fd2(i/2,j)=violations3;
    last_vals = conv_ord3_df2(max(end-11,1):end);
    mat_conv3_fd2(i/2, j) = {last_vals};


    end
end


%%
num_initial_points = N + 1;

% Crea una figura
figure;
hold on;


% Plot per ciascuna condizione iniziale
for j = 1:num_initial_points
    % Estrai l'ordine di convergenza per la j-esima condizione iniziale
    conv_ord_ex = mat_conv3_ex(:,j); % Derivate esatte

    plot(1:12,conv_ord_ex, 'Color', 'b', 'LineWidth', 1.5);
    hold on;
    for i =1:6 
        conv_ord_fd1 = mat_conv3_fd1{i, j}; % Differenze finite classiche

       
        conv_ord_fd2 = mat_conv3_fd2{i, j}; % Differenze finite adattative
        plot(1:12,conv_ord_fd1, '-', 'Color', 'r', 'LineWidth', 1.5);

        hold on;
        plot(1:12,conv_ord_fd2, '-o', 'Color', 'g', 'LineWidth', 1.5);
        hold on;
    end

   
end

% Aggiungi titolo e legenda
title('F16, n=10^5, superlinear');
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

%% VIOLATION

% Metodo Exact (derivate esatte) - media unica
vec_times_ex_clean = vec_violations3_ex; % copia dei tempi
vec_times_ex_clean(vec_converged3_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_v3 = mean(vec_times_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_times_fd1_clean = mat_violations3_fd1;
mat_times_fd1_clean(mat_converged3_fd1 == 0) = NaN;
avg_fd1 = mean(mat_times_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_times_fd2_clean = mat_violations3_fd2;
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
data = [ fd1_vals, avg_exact_v3; fd2_vals, avg_exact_v3;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T18 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation violation  table (only for successful runs): F16, n=10^5, suplin');
disp(T18);


%% BT-SEQ
% Metodo Exact (derivate esatte) - media unica
vec_bt_ex_clean = vec_bt3_ex; % copia dei tempi
vec_bt_ex_clean(vec_converged3_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_bt3 = mean(vec_bt_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_bt_fd1_clean = mat_bt3_fd1;
mat_bt_fd1_clean(mat_converged3_fd1 == 0) = NaN;
avg_fd1 = mean(mat_bt_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_bt_fd2_clean = mat_bt3_fd2;
mat_bt_fd2_clean(mat_converged3_fd2 == 0) = NaN;
avg_fd2 = mean(mat_bt_fd2_clean, 2, 'omitnan'); % media per ogni h (6x1)

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
data = [ fd1_vals, avg_exact_bt3; fd2_vals, avg_exact_bt3;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T19 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation bt iteration table (only for successful runs): F16, n=10^5, suplin');
disp(T19);



%% CG-SEQ
% Metodo Exact (derivate esatte) - media unica
vec_bt_ex_clean = vec_cg_iter3_ex; % copia dei tempi
vec_bt_ex_clean(vec_converged3_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_cg3 = mean(vec_bt_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_bt_fd1_clean = mat_cg_iter3_fd1;
mat_bt_fd1_clean(mat_converged3_fd1 == 0) = NaN;
avg_fd1 = mean(mat_bt_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_bt_fd2_clean = mat_cg_iter3_fd2;
mat_bt_fd2_clean(mat_converged3_fd2 == 0) = NaN;
avg_fd2 = mean(mat_bt_fd2_clean, 2, 'omitnan'); % media per ogni h (6x1)

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
data = [ fd1_vals, avg_exact_cg3; fd2_vals, avg_exact_cg3;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T20 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation cg iteration table (only for successful runs): F16, n=10^5, suplin');
disp(T20);



%% Calcolo quanti a convergenza


% Creazione delle etichette per i valori di h
h_exponents = [2, 4, 6, 8, 10, 12];  % Solo valori di h (senza h=0)
h_labels = arrayfun(@(e) sprintf('h=1e-%d', e), h_exponents, 'UniformOutput', false);

% Preparazione dei dati per la tabella
% FD1 e FD2 hanno le medie per ogni h, mentre Exact è ripetuto in tutte le colonne
%exact_vals = [avg_exact, avg_exact]; % Esatto in tutte le colonne
fd1_vals = sum(mat_converged3_fd1,2)';  % Trasposto per allineare con le colonne
fd2_vals = sum(mat_converged3_fd2,2)';  % Trasposto per allineare con le colonne

% Costruzione della tabella
rowNames = {'FD1', 'FD2'};
columnNames = [ h_labels,'Exact']; % Prima colonna "Exact", poi gli h
data = [ fd1_vals, sum(vec_converged3_ex); fd2_vals, sum(vec_converged3_ex);];

% Creiamo la tabella con i nomi delle colonne e delle righe
T21 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Number of converged : F16, n=10^5, suplin');
disp(T21);

writetable(T7, 'results_f16_suplin.xlsx', 'Sheet', 'time_5','WriteRowNames', true);
writetable(T8, 'results_f16_suplin.xlsx', 'Sheet', 'niter_5','WriteRowNames', true);
writetable(T9, 'results_f16_suplin.xlsx', 'Sheet', 'f_val_5','WriteRowNames', true);
writetable(T18, 'results_f16_suplin.xlsx', 'Sheet', 'viol_5','WriteRowNames', true);
writetable(T19, 'results_f16_suplin.xlsx', 'Sheet', 'bt_5','WriteRowNames', true);
writetable(T20, 'results_f16_suplin.xlsx', 'Sheet', 'cg_5','WriteRowNames', true);
writetable(T21, 'results_f16_suplin.xlsx', 'Sheet', 'n_conv5','WriteRowNames', true);



%%
% Creazione della tabella
data = [avg_exact_t1, avg_exact_t2, avg_exact_t3;
        avg_exact_i1, avg_exact_i2, avg_exact_i3;
        avg_exact_f1, avg_exact_f2, avg_exact_f3;
        avg_exact_v1, avg_exact_v2, avg_exact_v3;
        avg_exact_bt1, avg_exact_bt2, avg_exact_bt3;
        avg_exact_cg1, avg_exact_cg2, avg_exact_cg3;
        sum(vec_converged1_ex),sum(vec_converged2_ex),sum(vec_converged3_ex)];

% Definizione delle intestazioni
rowNames = {'Average Time', 'Average Iter', 'Average fval','Violation','Average iter Bt','Average iter cg', 'N converged'};
columnNames = {'n=10^3', 'n=10^4', 'n=10^5'};

% Creazione tabella MATLAB
T_compare = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);
disp(T_compare)

% Salvataggio su Excel
writetable(T_compare, 'results_f16_suplin.xlsx', 'Sheet', 'ExactComparison', 'WriteRowNames', true);

