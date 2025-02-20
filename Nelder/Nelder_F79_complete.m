%% FUNZIONE F79
% setting parameters 
format long
rng(345989);
n = 10;
tol = 1e-13;      %tollerance
max_iter = 1e06;  %max iteration
rho = 1.1;        %reflection parameter  
mu = 2.5;         %expansion parameter  
gamma = 0.8;      %contraction parameter
sigma = 0.9;      %shirinking parameter
delta = 1;        %Initialization of the simplex

% function
F = @(x) F79(x);  

N=10; %number of starting points
x0 = ones(n, 1);  % starting point
Mat_points=repmat(x0,1,N+1);
rand_mat=2*(rand([n, N+1]) - 0.5); % random matrix between [-1,1]
Mat_points=Mat_points + rand_mat; % starting points 
%vector for saving times
times_10=zeros(1,N+1);
%vector for saving minimum point
vec_10=zeros(1,N+1);
%vector for saving iterations
vec_iter_10=zeros(1,N+1);

for j =1:N+1
    %applying the function F79 to the 11 starting points
    tic;
    [xk_79_10, fk_79_10, n_iter] = Nelder_mead(Mat_points(:,j), F, rho, mu, gamma, sigma, tol, max_iter,delta);
    %saving results
    times_10(j)=toc;
    vec_10(j)=fk_79_10(end);
    vec_iter_10(j)=n_iter;
end

% Creation of a table with the results
results_n10 = table(["x0"; "x1"; "x2"; "x3"; "x4"; "x5"; "x6"; "x7"; "x8"; "x9"; "x10"], ...
                     times_10', vec_10', vec_iter_10', ...
                     'VariableNames', {'Initial condition', 'Time', 'FinalValue', 'Iterations'});
% Computation of mean of the three values saved
mean_time = mean(results_n10.Time);
mean_final_value = mean(results_n10.FinalValue);
mean_iterations = mean(results_n10.Iterations);

% Insert the mean in the tables
mean_row = table("Mean", mean_time, mean_final_value, mean_iterations, ...
                 'VariableNames', results_n10.Properties.VariableNames);
results_n10 = [results_n10; mean_row];

% Display the table
disp(results_n10);
% Creation an excel table
writetable(results_n10, 'Risultati_F79_Nelder.xlsx', 'Sheet', 'n_10');




%% FUNZIONE F79 n=25
%the same structure of n=10
format long
rng(345989);
n = 25;
tol = 1e-13;       
max_iter = 1e06;  
rho = 1.1;          
mu = 1.9;           
gamma = 0.8;       
sigma = 0.9;       
delta = 0.1; 

F = @(x) F79(x);  

N=10; 
x0 = ones(n, 1);
Mat_points=repmat(x0,1,N+1);
rand_mat=2*(rand([n, N+1]) - 0.5);
Mat_points=Mat_points + rand_mat;
times_25=zeros(1,N+1);
vec_25=zeros(1,N+1);
vec_iter_25=zeros(1,N+1);

for j =1:N+1
    tic;
    [xk_79_25, fk_79_25, n_iter] = Nelder_mead(Mat_points(:,j), F, rho, mu, gamma, sigma, tol, max_iter,delta);
    times_25(j)=toc;
    vec_25(j)=fk_79_25(end);
    vec_iter_25(j)=n_iter;
end

results_n25 = table(["x0"; "x1"; "x2"; "x3"; "x4"; "x5"; "x6"; "x7"; "x8"; "x9"; "x10"], ...
                     times_25', vec_25', vec_iter_25', ...
                     'VariableNames', {'Problem', 'Time', 'FinalValue', 'Iterations'});

mean_time = mean(results_n25.Time);
mean_final_value = mean(results_n25.FinalValue);
mean_iterations = mean(results_n25.Iterations);

mean_row = table("Mean", mean_time, mean_final_value, mean_iterations, ...
                 'VariableNames', results_n25.Properties.VariableNames);

results_n25 = [results_n25; mean_row];

disp(results_n25);

writetable(results_n25, 'Risultati_F79_Nelder.xlsx', 'Sheet', 'n_25');


%% FUNZIONE F79 n=50
% The same structure of n=10
format long
rng(345989);
n = 50;
tol = 1e-13;       
max_iter = 1e06;  
rho = 1.1;          
mu = 1.8;           
gamma = 0.8;      
sigma = 0.9;      
delta = 0.1; 

F = @(x) F79(x);  

N=10; 
x0 = ones(n, 1); 
Mat_points=repmat(x0,1,N+1);
rand_mat=2*(rand([n, N+1]) - 0.5);
Mat_points=Mat_points + rand_mat;
times_50=zeros(1,N+1);
vec_50=zeros(1,N+1);
vec_iter_50=zeros(1,N+1);

for j =1:N+1
    tic;
    [xk_79_50, fk_79_50, n_iter] = Nelder_mead(Mat_points(:,j), F, rho, mu, gamma, sigma, tol, max_iter,delta);
    times_50(j)=toc;
    vec_50(j)=fk_79_50(end);
    vec_iter_50(j)=n_iter;
end

results_n50 = table(["x0"; "x1"; "x2"; "x3"; "x4"; "x5"; "x6"; "x7"; "x8"; "x9"; "x10"], ...
                     times_50', vec_50', vec_iter_50', ...
                     'VariableNames', {'Problem', 'Time', 'FinalValue', 'Iterations'});

mean_time = mean(results_n50.Time);
mean_final_value = mean(results_n50.FinalValue);
mean_iterations = mean(results_n50.Iterations);

mean_row = table("Mean", mean_time, mean_final_value, mean_iterations, ...
                 'VariableNames', results_n50.Properties.VariableNames);

results_n50 = [results_n50; mean_row];

disp(results_n50);

writetable(results_n50, 'Risultati_F79_Nelder.xlsx', 'Sheet', 'n_50');

disp('All the results have been saved in Risultati_F79_Nelder.xlsx.');
