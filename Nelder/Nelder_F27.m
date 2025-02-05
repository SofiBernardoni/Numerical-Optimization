%% FUNZIONE F27
% partiamo da n=10
%inserisco tutti i parametri
format long
rng(345989);
n = 10;
tol = 1e-14;       % Tolleranza per la convergenza
max_iter = 1e08;  % Numero massimo di iterazioni
rho = 1.1;          % Parametro di espansione 
mu = 2.5;           % Parametro di contrazione
gamma = 0.8;      % Parametro di riflessione 
sigma = 0.9;      % Parametro di riduzione 
delta = 1; % delta del simplesso iniziale

% Definizione della funzione F27 come handle
F = @(x) F27(x);  % Passa x e n alla funzione F27semilogy

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

for j =1:N+1
    tic;
    [xk_27_10, fk_27_10, n_iter] = Nelder_mead(Mat_points(:,j), F, rho, mu, gamma, sigma, tol, max_iter,delta);
    times_10(j)=toc;
    vec_10(j)=fk_27_10(end);
    vec_iter_10(j)=n_iter;
end

% Raccolta dei dati per le tabelle
results_n10 = {
    'Problem','Times','Valore finale','Iterazioni';
    'x0', times_10(1), vec_10(1), vec_iter_10(1);
    'x1', times_10(2), vec_10(2), vec_iter_10(2);
    'x2', times_10(3), vec_10(3), vec_iter_10(3);
    'x3', times_10(4), vec_10(4), vec_iter_10(4);
    'x4', times_10(5), vec_10(5), vec_iter_10(5);
    'x5', times_10(6), vec_10(6), vec_iter_10(6);
    'x6', times_10(7), vec_10(7), vec_iter_10(7);
    'x7', times_10(8), vec_10(8), vec_iter_10(8);
    'x8', times_10(9), vec_10(9), vec_iter_10(9);
    'x9', times_10(10), vec_10(10), vec_iter_10(10);
    'x10', times_10(11), vec_10(11), vec_iter_10(11);

};

% Visualizzazione delle tabelle
disp('Risultati per n = 10:');
disp(results_n10);

subplot(3, 3, 1);
plot(1:length(times_10), times_10, '-.', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'Tempo F16 e n=10');
subplot(3,3,2);
plot(1:length(vec_10), vec_10, '-.', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'Valore minimo F16 e n=10');
subplot(3,3,3);
plot(1:length(vec_iter_10), vec_iter_10, '-.', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'Numero iterazioni F16 e n=10');

%% FUNZIONE F27 n=25
% partiamo da n=25
%inserisco tutti i parametri
format long
rng(345989);
n = 25;
tol = 1e-14;       % Tolleranza per la convergenza
max_iter = 1e08;  % Numero massimo di iterazioni
rho = 1.1;          % Parametro di espansione 
mu = 1.8;           % Parametro di contrazione
gamma = 0.8;      % Parametro di riflessione 
sigma = 0.9;      % Parametro di riduzione 
delta = 0.1; % delta del simplesso iniziale

% Definizione della funzione F27 come handle
F = @(x) F27(x);  % Passa x e n alla funzione F27semilogy




N=10; %numero di punti iniziali da generare
x0 = ones(n, 1);  % Punto iniziale
%voglio creare i  punti iniziali:
Mat_points=repmat(x0,1,N+1);
rand_mat=2*(rand([n, N+1]) - 0.5);
Mat_points=Mat_points + rand_mat;
%vettore per salvare i tempi
times_25=zeros(1,N+1);
%vettore per salvare i minimi
vec_25=zeros(1,N+1);
vec_iter_25=zeros(1,N+1);

for j =1:N+1
    tic;
    [xk_27_25, fk_27_25, n_iter] = Nelder_mead(Mat_points(:,j), F, rho, mu, gamma, sigma, tol, max_iter,delta);
    times_25(j)=toc;
    vec_25(j)=fk_27_25(end);
    vec_iter_25(j)=n_iter;
end

% Raccolta dei dati per le tabelle
results_n25 = {
    'Problem','Times','Valore finale','Iterazioni';
    'x0', times_25(1), vec_25(1), vec_iter_25(1);
    'x1', times_25(2), vec_25(2), vec_iter_25(2);
    'x2', times_25(3), vec_25(3), vec_iter_25(3);
    'x3', times_25(4), vec_25(4), vec_iter_25(4);
    'x4', times_25(5), vec_25(5), vec_iter_25(5);
    'x5', times_25(6), vec_25(6), vec_iter_25(6);
    'x6', times_25(7), vec_25(7), vec_iter_25(7);
    'x7', times_25(8), vec_25(8), vec_iter_25(8);
    'x8', times_25(9), vec_25(9), vec_iter_25(9);
    'x9', times_25(10), vec_25(10), vec_iter_25(10);
    'x10', times_25(11), vec_25(11), vec_iter_25(11);

};

% Visualizzazione delle tabelle
disp('Risultati per n = 25:');
disp(results_n25);

subplot(3, 3, 4);
plot(1:length(times_25), times_25, '-.', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'Tempo F27 e n=25');
subplot(3,3,5);
plot(1:length(vec_25), vec_25, '-.', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'Valore minimo F27 e n=25');
subplot(3,3,6);
plot(1:length(vec_iter_25), vec_iter_25, '-.', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'Numero iterazioni F27 e n=25');


%% FUNZIONE F27 n=50
% partiamo da n=50
%inserisco tutti i parametri
format long
rng(345989);
n = 50;
tol = 1e-10;       % Tolleranza per la convergenza
max_iter = 1e06;  % Numero massimo di iterazioni
rho = 1.1;          % Parametro di espansione 
mu = 1.8;           % Parametro di contrazione
gamma = 0.8;      % Parametro di riflessione 
sigma = 0.9;      % Parametro di riduzione 
delta = 0.1; % delta del simplesso iniziale

% Definizione della funzione F27 come handle
F = @(x) F27(x);  % Passa x e n alla funzione F27semilogy




N=10; %numero di punti iniziali da generare
x0 = ones(n, 1);  % Punto iniziale
%voglio creare i  punti iniziali:
Mat_points=repmat(x0,1,N+1);
rand_mat=2*(rand([n, N+1]) - 0.5);
Mat_points=Mat_points + rand_mat;
%vettore per salvare i tempi
times_50=zeros(1,N+1);
%vettore per salvare i minimi
vec_50=zeros(1,N+1);
vec_iter_50=zeros(1,N+1);

for j =1:N+1
    tic;
    [xk_27_50, fk_27_50, n_iter] = Nelder_mead(Mat_points(:,j), F, rho, mu, gamma, sigma, tol, max_iter,delta);
    times_50(j)=toc;
    vec_50(j)=fk_27_50(end);
    vec_iter_50(j)=n_iter;
end

% Raccolta dei dati per le tabelle
results_n50 = {
    'Problem','Times','Valore finale','Iterazioni';
    'x0', times_50(1), vec_50(1), vec_iter_50(1);
    'x1', times_50(2), vec_50(2), vec_iter_50(2);
    'x2', times_50(3), vec_50(3), vec_iter_50(3);
    'x3', times_50(4), vec_50(4), vec_iter_50(4);
    'x4', times_50(5), vec_50(5), vec_iter_50(5);
    'x5', times_50(6), vec_50(6), vec_iter_50(6);
    'x6', times_50(7), vec_50(7), vec_iter_50(7);
    'x7', times_50(8), vec_50(8), vec_iter_50(8);
    'x8', times_50(9), vec_50(9), vec_iter_50(9);
    'x9', times_50(10), vec_50(10), vec_iter_50(10);
    'x10', times_50(11), vec_50(11), vec_iter_50(11);

};

% Visualizzazione delle tabelle
disp('Risultati per n = 50:');
disp(results_n50);

subplot(3, 3, 7);
plot(1:length(times_50), times_50, '-.', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'Tempo F27 e n=50');
subplot(3,3,8);
plot(1:length(vec_50), vec_50, '-.', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'Valore minimo F27 e n=50');
subplot(3,3,9);
plot(1:length(vec_iter_50), vec_iter_50, '-.', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'Numero iterazioni F27 e n=50');

