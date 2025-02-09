rng(345989);
format short;

% Rosenbrock function
f = @(x) 100*(x(2,:) - x(1,:).^2).^2 + (1 - x(1,:)).^2;

% initial points
x0 = [1.2; 1.2];  
x1 = [-1.2; 1];
tol = 1e-14;  %tolerance
max_iter = 1e5; % Number of max iteration 

% parameters with some tuning
parametri=[
    1,   2,   0.5, 0.5, 0.1;
    1.2, 4,   0.7, 0.3, 0.1;
    1.2, 4,   0.7, 0.3, 0.5;
    1.4, 5,   0.8, 0.2, 0.1;
    1.65, 4.55, 0.95, 0.15, 0.1;
    2,   4,   0.7, 0.5, 0.5;
    ];
% Cration of a table
results = table;

%test the function with several parameteres in the two different inital
%points
for i = 1:size(parametri, 1)
    rho = parametri(i, 1);
    chi = parametri(i, 2);
    gamma = parametri(i, 3);
    sigma = parametri(i, 4);
    Delta = parametri(i, 5);  
    
    tic;
    [x_0, f_0, n_iter_0] = Nelder_mead(x0, f, rho, chi, gamma, sigma, tol, max_iter, Delta);
    tempo_x0=toc;
    tic;
    [x_1, f_1, n_iter_1] = Nelder_mead(x1, f, rho, chi, gamma, sigma, tol, max_iter, Delta);
    tempo_x1=toc;

    % Add the result in the table
    results = [results; table(rho, chi, gamma, sigma, Delta, ...
                x_0(1,end), x_0(2,end), f_0(end), n_iter_0, tempo_x0, ...
                x_1(1,end), x_1(2,end), f_1(end), n_iter_1, tempo_x1)];
end

% columns' name
results.Properties.VariableNames = {'Rho', 'Chi', 'Gamma', 'Sigma', 'Delta', ...
    'X_0(1)', 'X_0(2)', 'F_0', 'Iter_0','Time_0' ...
    'X_1(1)', 'X_1(2)', 'F_1', 'Iter_1', 'Time_1'};


% disp table
disp('Table with the results of Nelder-Mead method with Rosenbrock function ');
disp(results);

% save the results ona cvs file
writetable(results, 'Nelder_Rosenbrock_with_Delta.xlsx', 'WriteRowNames', true);
