rng(345989);

f=@(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

% Punto iniziale
x0 = [1.2; 1.2];  

% Parametri di Nelder-Mead
rho = 1;    % Coefficiente di riflessione
mu = 0.5;   % Coefficiente di contrazione
gamma = 2;  % Coefficiente di espansione
sigma = 0.5; % Coefficiente di riduzione
tol = 1e-16;  % Tolleranza per la convergenza
max_iter = 1000000; % Numero massimo di iterazioni

% Chiamata alla funzione Nelder-Mead
[x_opt, f_opt, n_iter] = Nelder_mead(x0, f, rho, mu, gamma, sigma, tol, max_iter);

% Visualizzare il risultato
disp('Punto ottimale:');
disp(x_opt);
disp('Valore della funzione ottimizzata:');
disp(f_opt);
disp('Numero di iterazioni:');
disp(n_iter);


x1 = [-1.2; 1]; 


% Chiamata alla funzione Nelder-Mead
[x_opt_1, f_opt_1, n_iter_1] = Nelder_mead(x1, f, rho, mu, gamma, sigma, tol, max_iter);

% Visualizzare il risultato
disp('Punto ottimale:');
disp(x_opt_1);
disp('Valore della funzione ottimizzata:');
disp(f_opt_1);
disp('Numero di iterazioni:');
disp(n_iter_1);


