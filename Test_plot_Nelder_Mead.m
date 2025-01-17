%% FILE CON N=10, arrivano tutti a convergenza

%% PROBLEM 16 1.1 2.5 0.8 0.9 sono 664 iterazioni
format long
%rng(345989);
%questa potrebbe avere minimi locali e il metodo non è in grado di trovare
%quelli globali
n = 10;
x0 = ones(n, 1);  % Punto iniziale
tol = 1e-14;       % Tolleranza per la convergenza
max_iter = 1e08;  % Numero massimo di iterazioni
rho = 1.1;          % Parametro di espansione 
mu = 2.5;           % Parametro di contrazione
gamma = 0.8;      % Parametro di riflessione 
sigma = 0.9;      % Parametro di riduzione 
delta=1; % delta del simplesso iniziale

% Definizione della funzione F16 come handle
F = @(x) F16(x);  % Passa x e n alla funzione F16semilogy

% Chiamata del metodo Nelder_mead

options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'TolFun', tol);
[x_min, F_min] = fminunc(F, x0, options);
[xk_16_10, fk_16_10, n_iter] = Nelder_mead(x0, F, rho, mu, gamma, sigma, tol, max_iter,delta);

%errore relativo
%err_rel_16_10=abs(F_min -fk_16_10)/abs(F_min);
err_rel_16_10=abs(F_min -fk_16_10);
%semilogy(1:n_iter+1, err_rel_16_10)
%hold on


%plot(1:n_iter+1, fk)
format long
disp(['min reale vale', num2str(F_min), 'min trovato da me ', num2str(fk_16_10(end)), 'con ', num2str(n_iter),'iterazioni' ])
%disp(xk_16_10(:,end))
%disp(x_min);
% problema: non arriva dove deve perchè cade in un minimo locale il
% simplesso rimae li per come è costriuto l'algoritmo
disp(fk_16_10(612))
disp(fk_16_10(end))


colors = lines(3); % Tre colori distinti per i tre problemi in ogni grafico

% Grafico per n = 10
figure;
%hold on;
semilogy(1:length(err_rel_16_10), err_rel_16_10, '-', 'Color', colors(1, :), 'LineWidth', 1.5, 'DisplayName', 'Problem 16');

%% PROBLEM 27 n=10 e 0.1 a tutti --> 2020 iter e 10^-5
rng(345989);
n=10;
x0=(1:n)';
F=@(x) F27(x);

tol = 1e-10;       % Tolleranza per la convergenza
max_iter = 1e08;  % Numero massimo di iterazioni
rho = 1.1;          % Parametro di espansione
mu = 2.1;           % Parametro di contrazione
gamma = 0.6;      % Parametro di riflessione
sigma = 0.6;      % Parametro di riduzione
delta=1; %delta del simplesso iniziale

%[xk_27_10, fk, n_iter] = Nelder_mead(x0, F, rho, mu, gamma, sigma, tol, max_iter);
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'TolFun', tol);
[x_min, F_min] = fminunc(F, x0, options);
[xk_27_10, fk_27_10, n_iter] = Nelder_mead(x0, F, rho, mu, gamma, sigma, tol, max_iter, delta);

%errore relativo
err_rel_27_10=abs(F_min -fk_27_10);
%err_rel_27_10=abs(F_min -fk_27_10)/abs(F_min);
%semilogy(1:n_iter+1, err_rel_27_10)
%hold on
disp(['min reale vale ', num2str(F_min), ' min trovato da me ', num2str(fk_27_10(end)), ' con', num2str(n_iter), 'iterazioni' ])


%% PROBLEM 79 NON MODIFICARE PARAMETRI!!!! con n=10 rho=1.1, mu=2.7, gamma=0.8, sigma=0.6
%dubbiooo, trovo un valore minore con 836 iterazioni
rng(345989);
n=10;
x0=-1*ones(n,1);
F= @(x) F79(x);

tol = 1e-13;       % Tolleranza per la convergenza
max_iter = 1e08;  % Numero massimo di iterazioni
rho = 1.1;          % Parametro di espansione
mu = 2.7;           % Parametro di contrazione
gamma = 0.8;      % Parametro di riflessione
sigma = 0.6;      % Parametro di riduzione
delta=1; % delta del simplesso iniziale



options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'TolFun', tol);
[x_min, F_min] = fminunc(F, x0, options);
[xk_79_10, fk_79_10, n_iter] = Nelder_mead(x0, F, rho, mu, gamma, sigma, tol, max_iter, delta);

%errore relativo
%err_rel_79_10=abs(F_min -fk_79_10)/abs(F_min);
err_rel_79_10=abs(F_min -fk_79_10);
%semilogy(1:n_iter+1, err_rel_79_10)
%hold on

disp(['min reale vale ', num2str(F_min), ' min trovato da me ', num2str(fk_79_10(end)), ' con', num2str(n_iter), 'iterazioni' ])


%% PROBLEM 16 1.1 1.8 0.8 0.9 sono 2152 iterazioni
rng(345989);
%questa potrebbe avere minimi locali e il metodo non è in grado di trovare
%quelli globali
n = 25;
x0 = ones(n, 1);  % Punto iniziale
tol = 1e-10;       % Tolleranza per la convergenza
max_iter = 1e08;  % Numero massimo di iterazioni
rho = 1.1;          % Parametro di espansione  
mu = 1.8;           % Parametro di contrazione
gamma = 0.8;      % Parametro di riflessione 
sigma = 0.9;      % Parametro di riduzione
delta=1; %delta del simplesso iniziale

% Definizione della funzione F16 come handle
F = @(x) F16(x);  % Passa x e n alla funzione F16

% Chiamata del metodo Nelder_mead

options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'TolFun', tol);
[x_min, F_min] = fminunc(F, x0, options);
[xk_16_25, fk_16_25, n_iter] = Nelder_mead(x0, F, rho, mu, gamma, sigma, tol, max_iter, delta);

%errore relativo
%err_rel_16_25=abs(F_min -fk_16_25)/abs(F_min);
err_rel_16_25=abs(F_min -fk_16_25);
%semilogy(1:n_iter+1, err_rel_16_25)
%hold on

disp(['min reale vale ', num2str(F_min), ' min trovato da me ', num2str(fk_16_25(end)), ' con', num2str(n_iter), 'iterazioni' ])

% problema: non arriva dove deve perchè cade in un minimo locale il
% simplesso rimae li per come è costriuto l'algoritmo

%% PROBLEM 27 n=25 1, 2.1, 0.6, 0.6 quasi giusto ma non perfetto
rng(345989);
n=25;
x0=(1:n)';
F=@(x) F27(x);

tol = 1e-14;       % Tolleranza per la convergenza
max_iter = 1e08;  % Numero massimo di iterazioni
rho = 1;          % Parametro di espansione
mu = 2.1;           % Parametro di contrazione
gamma = 0.6;      % Parametro di riflessione
sigma = 0.6;      % Parametro di riduzione
delta=1; %delta del simplesso iniziale


options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'TolFun', tol);
[x_min, F_min] = fminunc(F, x0, options);

[xk_27_25, fk_27_25, n_iter] = Nelder_mead(x0, F, rho, mu, gamma, sigma, tol, max_iter,delta);

%errore relativo
%err_rel_27_25=abs(F_min -fk_27_25)/abs(F_min);
err_rel_27_25=abs(F_min -fk_27_25);
%semilogy(1:n_iter+1, err_rel_27_25)
%hold on
disp(['min reale vale ', num2str(F_min), ' min trovato da me ', num2str(fk_27_25(end)), ' con ', num2str(n_iter), 'iterazioni' ])


%% PROBLEM 79 confusa... 1.1, 2.7,0.8,0.6 trovo un valore minore che lalgoritmo con 14058
rng(345989);
n=25;
x0=-1*ones(n,1);
F= @(x) F79(x);

tol = 1e-13;       % Tolleranza per la convergenza
max_iter = 1e08;  % Numero massimo di iterazioni
rho = 1.1;          % Parametro di espansione
mu = 2.7;           % Parametro di contrazione
gamma = 0.8;      % Parametro di riflessione
sigma = 0.6;      % Parametro di riduzione
delta=1; %delta del simplesso iniziale



options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'TolFun', tol);
[x_min, F_min] = fminunc(F, x0, options);
[xk_79_25, fk_79_25, n_iter] = Nelder_mead(x0, F, rho, mu, gamma, sigma, tol, max_iter, delta);

%errore relativo
%err_rel_79_25=abs(F_min -fk_79_25)/abs(F_min);
err_rel_79_25=abs(F_min -fk_79_25);
%semilogy(1:n_iter+1, err_rel_79_25)
%hold on
disp(['min reale vale ', num2str(F_min), ' min trovato da me ', num2str(fk_79_25(end)), ' con', num2str(n_iter), 'iterazioni' ])


%% PROBLEM 16 1.1 1.8 0.8 0.9 sono 152 iterazioni
rng(345989);
%questa potrebbe avere minimi locali e il metodo non è in grado di trovare
%quelli globali
n = 50;
x0 = ones(n, 1);  % Punto iniziale
tol = 1e-10;       % Tolleranza per la convergenza
max_iter = 1e08;  % Numero massimo di iterazioni
rho = 1.1;          % Parametro di espansione  
mu = 1.8;           % Parametro di contrazione
gamma = 0.8;      % Parametro di riflessione 
sigma = 0.9;      % Parametro di riduzione 
delta=1; %delta del simplesso iniziale
% Definizione della funzione F16 come handle
F = @(x) F16(x);  % Passa x e n alla funzione F16

% Chiamata del metodo Nelder_mead

options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'TolFun', tol);
[x_min, F_min] = fminunc(F, x0, options);
[xk_16_50, fk_16_50, n_iter] = Nelder_mead(x0, F, rho, mu, gamma, sigma, tol, max_iter, delta);

%errore relativo
%err_rel_16_50=abs(F_min -fk_16_50)/abs(F_min);
err_rel_16_50=abs(F_min -fk_16_50);
%semilogy(1:n_iter+1, err_rel_16_50)
%hold on

disp(['min reale vale ', num2str(F_min), ' min trovato da me ', num2str(fk_16_50(end)), ' con', num2str(n_iter), 'iterazioni' ])

% problema: non arriva dove deve perchè cade in un minimo locale il
% simplesso rimae li per come è costriuto l'algoritmo

%% PROBLEM 27 n=25 1, 2.1, 0.6, 0.6 quasi giusto ma non perfetto
rng(345989);
n=50;
x0=(1:n)';
F=@(x) F27(x);

tol = 1e-14;       % Tolleranza per la convergenza
max_iter = 1e08;  % Numero massimo di iterazioni
rho = 1;          % Parametro di espansione
mu = 2.1;           % Parametro di contrazione
gamma = 0.6;      % Parametro di riflessione
sigma = 0.6;      % Parametro di riduzione
delta=1; %delta del simplesso iniziale

options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'TolFun', tol);
[x_min, F_min] = fminunc(F, x0, options);

[xk_27_50, fk_27_50, n_iter] = Nelder_mead(x0, F, rho, mu, gamma, sigma, tol, max_iter, delta);
%errore relativo
%err_rel_27_50=abs(F_min -fk_27_50)/abs(F_min);
err_rel_27_50=abs(F_min -fk_27_50);
%semilogy(1:n_iter+1, err_rel_27_50)
%hold  on
disp(['min reale vale ', num2str(F_min), ' min trovato da me ', num2str(fk_27_50(end)), ' con ', num2str(n_iter), 'iterazioni' ])


%% PROBLEM 79 confusa... 1.1, 2.7,0.8,0.6 trovo un valore minore che lalgoritmo con 14058
%grandi problemiiiii
rng(345989);
n=50;
x0=-1*ones(n,1);
F= @(x) F79(x);

tol = 1e-13;       % Tolleranza per la convergenza
max_iter = 1e08;  % Numero massimo di iterazioni
rho = 1 ;          % Parametro di espansione
mu = 2;           % Parametro di contrazione
gamma = 0.5;      % Parametro di riflessione
sigma = 0.5;      % Parametro di riduzione
delta=1; % delta del simplesso iniziale



options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'TolFun', tol);
[x_min, F_min] = fminunc(F, x0, options);
[xk_79_50, fk_79_50, n_iter] = Nelder_mead(x0, F, rho, mu, gamma, sigma, tol, max_iter, delta);

%errore relativo
%err_rel_79_50=abs(F_min -fk_79_50)/abs(F_min);
err_rel_79_50=abs(F_min -fk_79_50);
%semilogy(1:n_iter+1, err_rel_79_50)
%hold on
disp(['min reale vale ', num2str(F_min), ' min trovato da me ', num2str(fk_79_50(end)), ' con', num2str(n_iter), 'iterazioni' ])


%% plot 
% Prepara i colori per ogni serie
colors = lines(3); % Tre colori distinti per i tre problemi in ogni grafico

% Grafico per n = 10
figure;
%hold on;
semilogy(1:length(err_rel_16_10), err_rel_16_10, '-', 'Color', colors(1, :), 'LineWidth', 1.5, 'DisplayName', 'Problem 16');
figure;
semilogy(1:length(err_rel_27_10), err_rel_27_10, '--', 'Color', colors(2, :), 'LineWidth', 1.5, 'DisplayName', 'Problem 27');
figure;
semilogy(1:length(err_rel_79_10), err_rel_79_10, '-.', 'Color', colors(3, :), 'LineWidth', 1.5, 'DisplayName', 'Problem 79');
xlabel('Numero di Iterazioni', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Errore Relativo (Scala Log)', 'FontSize', 14, 'FontWeight', 'bold');
title('Convergenza Metodo Nelder-Mead (n = 10)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
set(gcf, 'Color', 'w');
box on;
saveas(gcf, 'convergence_n10.png');
hold off;

%%
% Grafico per n = 25
figure;
%hold on;
semilogy(1:length(err_rel_16_25), err_rel_16_25, '-', 'Color', colors(1, :), 'LineWidth', 1.5, 'DisplayName', 'Problem 16');
figure;
semilogy(1:length(err_rel_27_25), err_rel_27_25, '--', 'Color', colors(2, :), 'LineWidth', 1.5, 'DisplayName', 'Problem 27');
figure;
semilogy(1:length(err_rel_79_25), err_rel_79_25, '-.', 'Color', colors(3, :), 'LineWidth', 1.5, 'DisplayName', 'Problem 79');
xlabel('Numero di Iterazioni', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Errore Relativo (Scala Log)', 'FontSize', 14, 'FontWeight', 'bold');
title('Convergenza Metodo Nelder-Mead (n = 25)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
set(gcf, 'Color', 'w');
box on;
saveas(gcf, 'convergence_n25.png');
hold off;


%%
% Grafico per n = 50
figure;
%hold on;
semilogy(1:length(err_rel_16_50), err_rel_16_50, '-', 'Color', colors(1, :), 'LineWidth', 1.5, 'DisplayName', 'Problem 16');
figure;
semilogy(1:length(err_rel_27_50), err_rel_27_50, '--', 'Color', colors(2, :), 'LineWidth', 1.5, 'DisplayName', 'Problem 27');
figure;
semilogy(1:length(err_rel_79_50), err_rel_79_50, '-.', 'Color', colors(3, :), 'LineWidth', 1.5, 'DisplayName', 'Problem 79');
xlabel('Numero di Iterazioni', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Errore Relativo (Scala Log)', 'FontSize', 14, 'FontWeight', 'bold');
title('Convergenza Metodo Nelder-Mead (n = 50)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
set(gcf, 'Color', 'w');
box on;
%saveas(gcf, 'convergence_n50.png');
hold off;


%%
% Raccolta dei dati per le tabelle
results_n10 = {
    'Problem', 'Iterazioni', 'Errore Relativo Finale';
    '16', length(err_rel_16_10)-1, err_rel_16_10(end);
    '27', length(err_rel_27_10)-1, err_rel_27_10(end);
    '79', length(err_rel_79_10)-1, err_rel_79_10(end);
};

results_n25 = {
    'Problem', 'Iterazioni', 'Errore Relativo Finale';
    '16', length(err_rel_16_25)-1, err_rel_16_25(end);
    '27', length(err_rel_27_25)-1, err_rel_27_25(end);
    '79', length(err_rel_79_25)-1, err_rel_79_25(end);
};

results_n50 = {
    'Problem', 'Iterazioni', 'Errore Relativo Finale';
    '16', length(err_rel_16_50)-1, err_rel_16_50(end);
    '27', length(err_rel_27_50)-1, err_rel_27_50(end);
    '79', length(err_rel_79_50)-1, err_rel_79_50(end);
};

% Visualizzazione delle tabelle
disp('Risultati per n = 10:');
disp(results_n10);

disp('Risultati per n = 25:');
disp(results_n25);

disp('Risultati per n = 50:');
disp(results_n50);

% Salvataggio dei risultati come file CSV per analisi esterne
%writecell(results_n10, 'results_n10.csv');
%writecell(results_n25, 'results_n25.csv');
%writecell(results_n50, 'results_n50.csv');
