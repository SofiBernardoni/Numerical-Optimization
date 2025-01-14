%%


x0 = linspace(-1,1,1000)';
plot(x0, F16(x0))



%% PROBLEM 16

%questa potrebbe avere minimi locali e il metodo non è in grado di trovare
%quelli globali
n = 50;
x0 = ones(n, 1);  % Punto iniziale
tol = 1e-10;       % Tolleranza per la convergenza
max_iter = 1e08;  % Numero massimo di iterazioni
rho = 1.2;          % Parametro di espansione
mu = 2.5;           % Parametro di contrazione
gamma = 0.7;      % Parametro di riflessione
sigma = 0.8;      % Parametro di riduzione

% Definizione della funzione F16 come handle
F = @(x) F16(x);  % Passa x e n alla funzione F16

% Chiamata del metodo Nelder_mead

options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'TolFun', tol);
[x_min, F_min] = fminunc(F, x0, options);
[xk, fk, n_iter] = Nelder_mead(x0, F, rho, mu, gamma, sigma, tol, max_iter, x_min);

disp(['min reale vale', num2str(F_min), 'min trovato da me ', num2str(fk), 'con ', num2str(n_iter),'iterazioni' ])


% problema: non arriva dove deve perchè cade in un minimo locale il
% simplesso rimae li per come è costriuto l'algoritmo

%% PROBLEM 27 n=10 e 0.1 a tutti --> 2020 iter e 10^-5
n=50;
x0=(1:n)';
F=@(x) F27(x);

tol = 1e-10;       % Tolleranza per la convergenza
max_iter = 1e08;  % Numero massimo di iterazioni
rho = 1.1;          % Parametro di espansione
mu = 2.1;           % Parametro di contrazione
gamma = 0.6;      % Parametro di riflessione
sigma = 0.6;      % Parametro di riduzione

[xk, fk, n_iter] = Nelder_mead(x0, F, rho, mu, gamma, sigma, tol, max_iter);
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'TolFun', tol);
[x_min, F_min] = fminunc(F, x0, options);

disp(['min reale vale', num2str(F_min), 'min trovato da me ', num2str(fk) ])


%% PROBLEM 79 NON MODIFICARE PARAMETRI!!!! con n=10 rho=1.1, mu=2.7, gamma=0.8, sigma=0.6
rng(347900);
n=10^5;
x0=-1*ones(n,1);
F= @(x) F79(x);

tol = 1e-13;       % Tolleranza per la convergenza
max_iter = 1e08;  % Numero massimo di iterazioni
rho = 1.1;          % Parametro di espansione
mu = 2.7;           % Parametro di contrazione
gamma = 0.8;      % Parametro di riflessione
sigma = 0.6;      % Parametro di riduzione



options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'TolFun', tol);
[x_min, F_min] = fminunc(F, x0, options);
[xk, fk, n_iter] = Nelder_mead(x0, F, rho, mu, gamma, sigma, tol, max_iter, x_min);
disp(['min reale vale', num2str(F_min), 'min trovato da me ', num2str(fk) ])

