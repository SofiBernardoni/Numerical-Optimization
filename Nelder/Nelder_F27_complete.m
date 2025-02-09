%% FUNZIONE F27 n=10
format long
rng(345989);
n = 10;
tol = 1e-14;
max_iter = 1e06;
rho = 1.1;
mu = 2.5;
gamma = 0.8;
sigma = 0.9;
delta = 1; % delta del simplesso iniziale

F = @(x) F27(x);

N = 10;
x0 = ones(n, 1);
Mat_points=repmat(x0,1,N+1);
rand_mat=2*(rand([n, N+1]) - 0.5);
Mat_points=Mat_points + rand_mat;

times_10 = zeros(1,N+1);
vec_10 = zeros(1,N+1);
vec_iter_10 = zeros(1,N+1);

for j = 1:N+1
    tic;
    [xk_27_10, fk_27_10, n_iter] = Nelder_mead(Mat_points(:,j), F, rho, mu, gamma, sigma, tol, max_iter, delta);
    times_10(j) = toc;
    vec_10(j) = fk_27_10(end);
    vec_iter_10(j) = n_iter;
    %disp(['fatta iterazione ', num2str(j)]);
end

results_n10 = table(["x0"; "x1"; "x2"; "x3"; "x4"; "x5"; "x6"; "x7"; "x8"; "x9"; "x10"], ...
                     times_10', vec_10', vec_iter_10', ...
                     'VariableNames', {'Problem', 'Time', 'FinalValue', 'Iterations'});

% Calcola la media delle colonne numeriche
mean_time = mean(results_n10.Time);
mean_final_value = mean(results_n10.FinalValue);
mean_iterations = mean(results_n10.Iterations);

% Crea una nuova riga con le medie
mean_row = table("Mean", mean_time, mean_final_value, mean_iterations, ...
                 'VariableNames', results_n10.Properties.VariableNames);

% Aggiungi la riga alla tabella originale
results_n10 = [results_n10; mean_row];

% Visualizza la tabella aggiornata
disp(results_n10);
% Creazione tabella per Excel

writetable(results_n10, 'Risultati_F27_Nelder.xlsx', 'Sheet', 'n_10');

%% FUNZIONE F27 n=25
n = 25;
tol = 1e-14;
max_iter = 1e08;
rho = 1.1;
mu = 1.8;
gamma = 0.8;
sigma = 0.9;
delta = 0.1;

F = @(x) F27(x);

x0 = ones(n, 1);
Mat_points=repmat(x0,1,N+1);
rand_mat=2*(rand([n, N+1]) - 0.5);
Mat_points=Mat_points + rand_mat;

times_25 = zeros(1,N+1);
vec_25 = zeros(1,N+1);
vec_iter_25 = zeros(1,N+1);

for j = 1:N+1
    tic;
    [xk_27_25, fk_27_25, n_iter] = Nelder_mead(Mat_points(:,j), F, rho, mu, gamma, sigma, tol, max_iter, delta);
    times_25(j) = toc;
    vec_25(j) = fk_27_25(end);
    vec_iter_25(j) = n_iter;
    %disp(['fatta iterazione ', num2str(j)]);
end

results_n25 = table(["x0"; "x1"; "x2"; "x3"; "x4"; "x5"; "x6"; "x7"; "x8"; "x9"; "x10"], ...
                     times_25', vec_25', vec_iter_25', ...
                     'VariableNames', {'Problem', 'Time', 'FinalValue', 'Iterations'});

% Calcola la media delle colonne numeriche
mean_time = mean(results_n25.Time);
mean_final_value = mean(results_n25.FinalValue);
mean_iterations = mean(results_n25.Iterations);

% Crea una nuova riga con le medie
mean_row = table("Mean", mean_time, mean_final_value, mean_iterations, ...
                 'VariableNames', results_n25.Properties.VariableNames);

% Aggiungi la riga alla tabella originale
results_n25 = [results_n25; mean_row];

% Visualizza la tabella aggiornata
disp(results_n25);
% Creazione tabella per Excel


writetable(results_n25, 'Risultati_F27_Nelder.xlsx', 'Sheet', 'n_25');

%% FUNZIONE F27 n=50
n = 50;
tol = 1e-10;
max_iter = 1e06;
rho = 1.1;
mu = 1.8;
gamma = 0.8;
sigma = 0.9;
delta = 0.1;

F = @(x) F27(x);

x0 = ones(n, 1);
Mat_points=repmat(x0,1,N+1);
rand_mat=2*(rand([n, N+1]) - 0.5);
Mat_points=Mat_points + rand_mat;

times_50 = zeros(1,N+1);
vec_50 = zeros(1,N+1);
vec_iter_50 = zeros(1,N+1);

for j = 1:N+1
    tic;
    [xk_27_50, fk_27_50, n_iter] = Nelder_mead(Mat_points(:,j), F, rho, mu, gamma, sigma, tol, max_iter, delta);
    times_50(j) = toc;
    vec_50(j) = fk_27_50(end);
    vec_iter_50(j) = n_iter;
    %disp(['fatta iterazione ', num2str(j)]);
end

results_n50 = table(["x0"; "x1"; "x2"; "x3"; "x4"; "x5"; "x6"; "x7"; "x8"; "x9"; "x10"], ...
                     times_50', vec_50', vec_iter_50', ...
                     'VariableNames', {'Problem', 'Time', 'FinalValue', 'Iterations'});

% Calcola la media delle colonne numeriche
mean_time = mean(results_n50.Time);
mean_final_value = mean(results_n50.FinalValue);
mean_iterations = mean(results_n50.Iterations);

% Crea una nuova riga con le medie
mean_row = table("Mean", mean_time, mean_final_value, mean_iterations, ...
                 'VariableNames', results_n50.Properties.VariableNames);

% Aggiungi la riga alla tabella originale
results_n50 = [results_n50; mean_row];

% Visualizza la tabella aggiornata
disp(results_n50);
% Creazione tabella per Excel

writetable(results_n50, 'Risultati_F27_Nelder.xlsx', 'Sheet', 'n_50');

disp('Tutti i risultati sono stati salvati in Risultati_F27.xlsx.');