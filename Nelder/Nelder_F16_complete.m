%% FUNZIONE F16 n=10
% Impostazione parametri iniziali
format long
rng(345989);
n = 10;
tol = 1e-14;
max_iter = 1e08;
rho = 1.1;
mu = 2.5;
gamma = 0.8;
sigma = 0.9;
delta = 1; % delta del simplesso iniziale

% Definizione della funzione F16
F = @(x) F16(x);



N=10; %numero di punti iniziali da generare
x0 = ones(n, 1);  % Punto iniziale
%voglio creare i  punti iniziali:
Mat_points=repmat(x0,1,N+1);
rand_mat=2*(rand([n, N+1]) - 0.5);
Mat_points=Mat_points + rand_mat;
%vettore per salvare i tempi
times_10=zeros(1,N+1);
%vettore per salvare i minimi
vec_10=zeros(1,N+1);
vec_iter_10=zeros(1,N+1);



for j = 1:N+1
    tic;
    [xk_16_10, fk_16_10, n_iter] = Nelder_mead(Mat_points(:,j), F, rho, mu, gamma, sigma, tol, max_iter, delta);
    times_10(j) = toc;
    vec_10(j) = fk_16_10(end);
    vec_iter_10(j) = n_iter;
end

% Creazione tabella per Excel
results_n10 = table(["x0"; "x1"; "x2"; "x3"; "x4"; "x5"; "x6"; "x7"; "x8"; "x9"; "x10"], ...
                     times_10', vec_10', vec_iter_10', ...
                     'VariableNames', {'Initial condition', 'Time', 'FinalValue', 'Iterations'});

% Scrittura su Excel
writetable(results_n10, 'Risultati_F16_Nelder.xlsx', 'Sheet', 'n_10');


%% FUNZIONE F16 n=25
% Impostazione parametri
n = 25;
tol = 1e-14;
max_iter = 1e08;
rho = 1.1;
mu = 1.8;
gamma = 0.8;
sigma = 0.9;
delta = 0.1; % delta del simplesso iniziale

F = @(x) F16(x);

x0 = ones(n, 1);
Mat_points = repmat(x0,1,N+1) + 2*(rand([n, N+1]) - 0.5);

times_25 = zeros(1,N+1);
vec_25 = zeros(1,N+1);
vec_iter_25 = zeros(1,N+1);

for j = 1:N+1
    tic;
    [xk_16_25, fk_16_25, n_iter] = Nelder_mead(Mat_points(:,j), F, rho, mu, gamma, sigma, tol, max_iter, delta);
    times_25(j) = toc;
    vec_25(j) = fk_16_25(end);
    vec_iter_25(j) = n_iter;
end

% Creazione tabella per Excel
results_n25 = table(["x0"; "x1"; "x2"; "x3"; "x4"; "x5"; "x6"; "x7"; "x8"; "x9"; "x10"], ...
                     times_25', vec_25', vec_iter_25', ...
                     'VariableNames', {'Problem', 'Time', 'FinalValue', 'Iterations'});

% Scrittura su Excel
writetable(results_n25, 'Risultati_F16_Nelder.xlsx', 'Sheet', 'n_25');


%% FUNZIONE F16 n=50
% Impostazione parametri
n = 50;
tol = 1e-10;
max_iter = 1e06;
rho = 1.1;
mu = 1.8;
gamma = 0.8;
sigma = 0.9;
delta = 0.1; % delta del simplesso iniziale

F = @(x) F16(x);

x0 = ones(n, 1);
Mat_points = repmat(x0,1,N+1) + 2*(rand([n, N+1]) - 0.5);

times_50 = zeros(1,N+1);
vec_50 = zeros(1,N+1);
vec_iter_50 = zeros(1,N+1);

for j = 1:N+1
    tic;
    [xk_16_50, fk_16_50, n_iter] = Nelder_mead(Mat_points(:,j), F, rho, mu, gamma, sigma, tol, max_iter, delta);
    times_50(j) = toc;
    vec_50(j) = fk_16_50(end);
    vec_iter_50(j) = n_iter;
end

% Creazione tabella per Excel
results_n50 = table(["x0"; "x1"; "x2"; "x3"; "x4"; "x5"; "x6"; "x7"; "x8"; "x9"; "x10"], ...
                     times_50', vec_50', vec_iter_50', ...
                     'VariableNames', {'Problem', 'Time', 'FinalValue', 'Iterations'});

% Scrittura su Excel
writetable(results_n50, 'Risultati_F16_Nelder.xlsx', 'Sheet', 'n_50');

disp('Tutti i risultati sono stati salvati in Risultati_F16.xlsx.');
