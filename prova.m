rng(345989);

f=@(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

% Punto iniziale
x0 = [1.2; 1.2];  

% Parametri di Nelder-Mead
rho = 1;    % Coefficiente di riflessione
mu = 2;   % Coefficiente di contrazione
gamma = 0.5;  % Coefficiente di espansione
sigma = 0.5; % Coefficiente di riduzione
tol = 1e-16;  % Tolleranza per la convergenza
max_iter = 1000; % Numero massimo di iterazioni

% Chiamata alla funzione Nelder-Mead
[x_0, f_0, n_iter_0] = Nelder_mead(x0, f, rho, mu, gamma, sigma, tol, max_iter);

% Visualizzare il risultato

disp('-------------NELDER-MEAD CON VALORI STANDARD---------------')
disp('Punto ottimale da x0 = [1.2, 1.2]:');
disp(x_0);
disp('Valore della funzione ottimizzata:');
disp(f_0);
disp('Numero di iterazioni:');
disp(n_iter_0);


x1 = [-1.2; 1]; 


% Chiamata alla funzione Nelder-Mead
[x_1, f_1, n_iter_1] = Nelder_mead(x1, f, rho, mu, gamma, sigma, tol, max_iter);

% Visualizzare il risultato
disp('Punto ottimale da x1=[-1.2,1]:');
disp(x_1);
disp('Valore della funzione ottimizzata:');
disp(f_1);
disp('Numero di iterazioni:');
disp(n_iter_1);



%% CAMBIO PARANTRI
disp('change dei parametri:');
disp('rho=1.1, gamma=0.6, sigma=0.4,mu=3')

% Parametri di Nelder-Mead
rho = 1.1;    % Coefficiente di riflessione
mu = 3;   % Coefficiente di contrazione
gamma = 0.6;  % Coefficiente di espansione
sigma = 0.4; % Coefficiente di riduzione

% Chiamata alla funzione Nelder-Mead
[x_1c, f_1c, n_iter_1c] = Nelder_mead(x1, f, rho, mu, gamma, sigma, tol, max_iter);
[x_0c, f_0c, n_iter_0c] = Nelder_mead(x0, f, rho, mu, gamma, sigma, tol, max_iter);

% Visualizzare il risultato
disp('----------------- NEALDER-MEAD CON NUOVI VALORI')
% Visualizzare il risultato
disp('Punto ottimale da x1=[-1.2,1]:');
disp(x_1c);
disp('Valore della funzione ottimizzata:');
disp(f_1c);
disp('Numero di iterazioni:');
disp(n_iter_1c);

% Visualizzare il risultato
disp('Punto ottimale da x0=[1.2,1.2]:');
disp(x_0c);
disp('Valore della funzione ottimizzata:');
disp(f_0c);
disp('Numero di iterazioni:');
disp(n_iter_0c);

%% CAMBIO PARAMETRI
disp('change dei parametri:');
disp('rho=1.2, gamma=0.7, sigma=0.3,mu=4')


% Parametri di Nelder-Mead
rho = 1.2;    % Coefficiente di riflessione
mu = 4;   % Coefficiente di contrazione
gamma = 0.7;  % Coefficiente di espansione
sigma = 0.3; % Coefficiente di riduzione

[x_1_1, f_1_1, n_iter_1_1] = Nelder_mead(x1, f, rho, mu, gamma, sigma, tol, max_iter);
[x_0_1, f_0_1, n_iter_0_1] = Nelder_mead(x0, f, rho, mu, gamma, sigma, tol, max_iter);


disp('----------------- NEALDER-MEAD CON NUOVI VALORI')
% Visualizzare il risultato
disp('Punto ottimale da x1=[-1.2,1]:');
disp(x_1_1);
disp('Valore della funzione ottimizzata:');
disp(f_1_1);
disp('Numero di iterazioni:');
disp(n_iter_1_1);

% Visualizzare il risultato
disp('Punto ottimale da x0=[1.2,1.2]:');
disp(x_0_1);
disp('Valore della funzione ottimizzata:');
disp(f_0_1);
disp('Numero di iterazioni:');
disp(n_iter_0_1);

%% NUOVO CAMBIO
disp('change dei parametri:');
disp('rho=1.2, gamma=1.7, sigma=0.3,mu=0.3')


% Parametri di Nelder-Mead
rho = 1.4;    % Coefficiente di riflessione
mu = 5;   % Coefficiente di contrazione
gamma = 0.8;  % Coefficiente di espansione
sigma = 0.2; % Coefficiente di riduzione

[x_1_1, f_1_1, n_iter_1_1] = Nelder_mead(x1, f, rho, mu, gamma, sigma, tol, max_iter);
[x_0_1, f_0_1, n_iter_0_1] = Nelder_mead(x0, f, rho, mu, gamma, sigma, tol, max_iter);


disp('----------------- NEALDER-MEAD CON NUOVI VALORI')
% Visualizzare il risultato
disp('Punto ottimale da x1=[-1.2,1]:');
disp(x_1_1);
disp('Valore della funzione ottimizzata:');
disp(f_1_1);
disp('Numero di iterazioni:');
disp(n_iter_1_1);

% Visualizzare il risultato
disp('Punto ottimale da x0=[1.2,1.2]:');
disp(x_0_1);
disp('Valore della funzione ottimizzata:');
disp(f_0_1);
disp('Numero di iterazioni:');
disp(n_iter_0_1);

%% NUOVO CAMBIO
disp('change dei parametri:');
disp('rho=1.5, gamma=3, sigma=0.2,mu=0.2')


% Parametri di Nelder-Mead
rho = 1.5;    % Coefficiente di riflessione
mu = 4.5;   % Coefficiente di contrazione
gamma = 0.8;  % Coefficiente di espansione
sigma = 0.1; % Coefficiente di riduzione

[x_1_1, f_1_1, n_iter_1_1] = Nelder_mead(x1, f, rho, mu, gamma, sigma, tol, max_iter);
[x_0_1, f_0_1, n_iter_0_1] = Nelder_mead(x0, f, rho, mu, gamma, sigma, tol, max_iter);


disp('----------------- NEALDER-MEAD CON NUOVI VALORI')
% Visualizzare il risultato
disp('Punto ottimale da x1=[-1.2,1]:');
disp(x_1_1);
disp('Valore della funzione ottimizzata:');
disp(f_1_1);
disp('Numero di iterazioni:');
disp(n_iter_1_1);

% Visualizzare il risultato
disp('Punto ottimale da x0=[1.2,1.2]:');
disp(x_0_1);
disp('Valore della funzione ottimizzata:');
disp(f_0_1);
disp('Numero di iterazioni:');
disp(n_iter_0_1);


%% NUOVO CAMBIO
disp('change dei parametri:');
disp('rho=1.6, gamma=0.9, sigma=0.1,mu=4.5')


% Parametri di Nelder-Mead
rho = 1.6;    % Coefficiente di riflessione
mu = 4.5;   % Coefficiente di contrazione
gamma = 0.9;  % Coefficiente di espansione
sigma = 0.1; % Coefficiente di riduzione

[x_1_1, f_1_1, n_iter_1_1] = Nelder_mead(x1, f, rho, mu, gamma, sigma, tol, max_iter);
[x_0_1, f_0_1, n_iter_0_1] = Nelder_mead(x0, f, rho, mu, gamma, sigma, tol, max_iter);


disp('----------------- NEALDER-MEAD CON NUOVI VALORI')
% Visualizzare il risultato
disp('Punto ottimale da x1=[-1.2,1]:');
disp(x_1_1);
disp('Valore della funzione ottimizzata:');
disp(f_1_1);
disp('Numero di iterazioni:');
disp(n_iter_1_1);

% Visualizzare il risultato
disp('Punto ottimale da x0=[1.2,1.2]:');
disp(x_0_1);
disp('Valore della funzione ottimizzata:');
disp(f_0_1);
disp('Numero di iterazioni:');
disp(n_iter_0_1);

%% NUOVO CAMBIO
disp('change dei parametri:');
disp('rho=1.6, gamma=0.9, sigma=0.1,mu=4.5')


% Parametri di Nelder-Mead
rho = 1.65;    % Coefficiente di riflessione
mu = 4.55;   % Coefficiente di contrazione
gamma = 0.95;  % Coefficiente di espansione
sigma = 0.15; % Coefficiente di riduzione

[x_1_1, f_1_1, n_iter_1_1] = Nelder_mead(x1, f, rho, mu, gamma, sigma, tol, max_iter);
[x_0_1, f_0_1, n_iter_0_1] = Nelder_mead(x0, f, rho, mu, gamma, sigma, tol, max_iter);


disp('----------------- NEALDER-MEAD CON NUOVI VALORI')
% Visualizzare il risultato
disp('Punto ottimale da x1=[-1.2,1]:');
disp(x_1_1);
disp('Valore della funzione ottimizzata:');
disp(f_1_1);
disp('Numero di iterazioni:');
disp(n_iter_1_1);

% Visualizzare il risultato
disp('Punto ottimale da x0=[1.2,1.2]:');
disp(x_0_1);
disp('Valore della funzione ottimizzata:');
disp(f_0_1);
disp('Numero di iterazioni:');
disp(n_iter_0_1);