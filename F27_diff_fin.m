%% FUNCTION 27 - PENALTY FUNCTION 1 (with specific truncated newton function useful for big dimensions) (with different initial points)
%differenze finite normali
rng(345989);

exact=false;
fin_dif_2=false;
h=[2,4,6,8,10,12];
h_values =10.^(-h);

F = @(x) F27(x);  % Definizione della funzione F27 come handle

% The hessian cannot be defined as sparse. It is directly used in the newton method (with a specific version suitable to function 27)

load forcing_terms.mat % termini per tolleranza adattiva
%% n=10^3 (1e3)

rng(345989);

n=1e3; 

kmax=1e3;
tolgrad=5e-7;
cg_maxit=50;

z0=zeros(n,1);
c1=1e-4;
rho=0.5;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

x0=(1:n)';  % Initial point
N=10; % number of initial points to be generated

% Initial points:
Mat_points=repmat(x0,1,N+1); 
rand_mat=2*rand(n, N)-1;
Mat_points(:,2:end)=Mat_points(:,2:end) + rand_mat; % matrix with columns=initial points

results_n3 = struct();  % Struttura per salvare i risultati

for h_idx = 1:length(h_values)
    h = h_values(h_idx);  % Assegna il valore corrente di h
    fprintf('Test con h = %.4e\n', h);

    % Vettori per salvare i risultati per ogni h
    vec_times1 = zeros(1, N+1);
    vec_val1 = zeros(1, N+1);
    vec_grad1 = zeros(1, N+1);
    vec_iter1 = zeros(1, N+1);
    %vec_converged1 = zeros(1, N+1);
    vec_violations1 = zeros(1, N+1);
    conv_order_cell1 = cell(1, N+1);

    JF = @(x) JF27(x,exact,fin_dif_2,h); % Definizione della funzione JF27 come handle (derivata esatta)

    for j = 1:N+1
        tic;
        [x1, f1, gradf_norm1, k1, xseq1, btseq1, cgiterseq1, conv_ord1, flag1, converged1, violations1] = ...
            truncated_newton_27(Mat_points(:,j), F, JF, exact, fin_dif_2, h, kmax, tolgrad, fterms_suplin, cg_maxit, z0, c1, rho, btmax);
        vec_times1(j) = toc;

        disp(['Tentativo n. ', num2str(j), ': ', flag1])
        %vec_converged1(j) = conv_ord1;
        conv_order_cell1{j} = conv_ord1;
        vec_val1(j) = f1;
        vec_grad1(j) = gradf_norm1;
        vec_iter1(j) = k1;
        vec_violations1(j) = violations1;
    end

    % Salva i risultati per questo valore di h
    results_n3(h_idx).h = h;
    results_n3(h_idx).times = vec_times1;
    results_n3(h_idx).values = vec_val1;
    results_n3(h_idx).gradients = vec_grad1;
    results_n3(h_idx).iterations = vec_iter1;
    %results_n3(h_idx).converged = vec_converged1;
    results_n3(h_idx).violations = vec_violations1;
end

for j = 1:length(conv_order_cell1)
    fprintf('Esecuzione %d: ordine di convergenza medio = %f\n', j, mean(conv_order_cell1{j}, 'omitnan'));
    %fprintf('Contenuto della cella %d:\n', j);
    %disp(conv_order_cell1{j});
end


%% n=10^4 (1e4)

rng(345989);

n=1e4; 

kmax=1.5e3;
tolgrad=5e-7;
cg_maxit=100; %50 uguale

z0=zeros(n,1);
c1=1e-4;
rho=0.5;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

x0=(1:n)';  % Initial point
N=10; % number of initial points to be generated

% Initial points:
Mat_points=repmat(x0,1,N+1); 
rand_mat=2*rand(n, N)-1;
Mat_points(:,2:end)=Mat_points(:,2:end) + rand_mat; % matrix with columns=initial points

% Definire i valori di h da testare
results_n4 = struct();  % Struttura per salvare i risultati

for h_idx = 1:length(h_values)
    h = h_values(h_idx);
    fprintf('Test con h = %.4e\n', h);
    
    % Vettori per salvare i risultati per ogni h
    vec_times2 = zeros(1, N+1);
    vec_val2 = zeros(1, N+1);
    vec_grad2 = zeros(1, N+1);
    vec_iter2 = zeros(1, N+1);
    %vec_converged2 = zeros(1, N+1);
    vec_violations2 = zeros(1, N+1);
    conv_order_cell2 = cell(1, N+1);

    for j = 1:N+1
        tic;
        [x2, f2, gradf_norm2, k2, xseq2, btseq2, cgiterseq2, conv_ord2, flag2, converged2, violations2] = ...
            truncated_newton_27(Mat_points(:,j), F, JF, exact, fin_dif_2, h, kmax, tolgrad, fterms_suplin, cg_maxit, z0, c1, rho, btmax);
        vec_times2(j) = toc;

        disp(['Tentativo n. ', num2str(j), ': ', flag2])
        %vec_converged2(j) = converged2;
        conv_order_cell2{j} = conv_ord2;
        vec_val2(j) = f2;
        vec_grad2(j) = gradf_norm2;
        vec_iter2(j) = k2;
        vec_violations2(j) = violations2;
    end

    % Salva i risultati per questo valore di h
    results_n4(h_idx).h = h;
    results_n4(h_idx).times = vec_times2;
    results_n4(h_idx).values = vec_val2;
    results_n4(h_idx).gradients = vec_grad2;
    results_n4(h_idx).iterations = vec_iter2;
    %results_n4(h_idx).converged = vec_converged2;
    results_n4(h_idx).violations = vec_violations2;
end
for j = 1:length(conv_order_cell2)
    fprintf('Esecuzione %d: ordine di convergenza medio = %f\n', j, mean(conv_order_cell2{j}, 'omitnan'));
end

%% n=10^5 (1e5)
rng(345989);

n = 1e5; 

kmax = 1.5e3; 
tolgrad = 5e-7; 
cg_maxit = 100; 

z0 = zeros(n,1);
c1 = 1e-4;
rho = 0.5;
btmax = 50; 

x0 = (1:n)';  % Initial point
N = 10; % Number of initial points to be generated

% Initial points:
Mat_points = repmat(x0,1,N+1); 
rand_mat = 2*rand(n, N)-1;
Mat_points(:,2:end) = Mat_points(:,2:end) + rand_mat; % matrix with columns = initial points

results_n5 = struct();  % Struttura per salvare i risultati

for h_idx = 1:length(h_values)
    h = h_values(h_idx);
    fprintf('Test con h = %.4e\n', h);
    
    % Vettori per salvare i risultati per ogni h
    vec_times3 = zeros(1, N+1);
    vec_val3 = zeros(1, N+1);
    vec_grad3 = zeros(1, N+1);
    vec_iter3 = zeros(1, N+1);
    %vec_converged3 = zeros(1, N+1);
    conv_order_cell3 = cell(1, N+1);
    vec_violations3 = zeros(1, N+1);

    for j = 1:N+1
        tic;
        [x3, f3, gradf_norm3, k3, xseq3, btseq3, cgiterseq3, conv_ord3, flag3, converged3, violations3] = ...
            truncated_newton_27(Mat_points(:,j), F, JF, exact, fin_dif_2, h, kmax, tolgrad, fterms_suplin, cg_maxit, z0, c1, rho, btmax);
        vec_times3(j) = toc;

        disp(['Tentativo n. ', num2str(j), ': ', flag3])
        %vec_converged3(j) = converged3;
        conv_order_cell3{j} = conv_ord3;
        vec_val3(j) = f3;
        vec_grad3(j) = gradf_norm3;
        vec_iter3(j) = k3;
        vec_violations3(j) = violations3;
    end

    % Salva i risultati per questo valore di h
    results_n5(h_idx).h = h;
    results_n5(h_idx).times = vec_times3;
    results_n5(h_idx).values = vec_val3;
    results_n5(h_idx).gradients = vec_grad3;
    results_n5(h_idx).iterations = vec_iter3;
    %results_n5(h_idx).converged = vec_converged3;
    results_n5(h_idx).violations = vec_violations3;
end
for j = 1:length(conv_order_cell3)
    fprintf('Esecuzione %d: ordine di convergenza medio = %f\n', j, mean(conv_order_cell3{j}, 'omitnan'));
end

%%
% Assicurarsi che tutte le variabili siano colonne
h_col = h_values(:); % Converte h_values in colonna
mean_time_col = cellfun(@mean, {results_n3.times})'; % Tempo medio
mean_value_col = cellfun(@mean, {results_n3.values})'; % Valore medio della funzione
mean_iter_col = cellfun(@mean, {results_n3.iterations})'; % Iterazioni medie

% Creazione della tabella
T = table(h_col, mean_time_col, mean_value_col, mean_iter_col, ...
    'VariableNames', {'h', 'Tempo_Medio', 'Valore_Medio_Funzione', 'Iterazioni_Medie'});

% Stampa la tabella in console
disp('===== Risultati Medi per h =====');
disp(T);




%% FUNCTION 27 - PENALTY FUNCTION 1 (with specific truncated newton function useful for big dimensions) (with different initial points)
%differenze finite di secondo tipo
rng(345989);

exact=false;
fin_dif_2=true;
h=[2,4,6,8,10,12];
h_values =10.^(-h);

F = @(x) F27(x);  % Definizione della funzione F27 come handle
JF = @(x) JF27(x,exact,0); % Definizione della funzione JF27 come handle (derivata esatta)
% The hessian cannot be defined as sparse. It is directly used in the newton method (with a specific version suitable to function 27)

load forcing_terms.mat % termini per tolleranza adattiva


%% n=10^3 (1e3)

rng(345989);

n=1e3; 

kmax=1e3;
tolgrad=5e-7;
cg_maxit=50;

z0=zeros(n,1);
c1=1e-4;
rho=0.5;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

x0=(1:n)';  % Initial point
N=10; % number of initial points to be generated

% Initial points:
Mat_points=repmat(x0,1,N+1); 
rand_mat=2*rand(n, N)-1;
Mat_points(:,2:end)=Mat_points(:,2:end) + rand_mat; % matrix with columns=initial points

results_n3 = struct();  % Struttura per salvare i risultati

for h_idx = 1:length(h_values)
    h = h_values(h_idx);  % Assegna il valore corrente di h
    fprintf('Test con h = %.4e\n', h);

    % Vettori per salvare i risultati per ogni h
    vec_times1 = zeros(1, N+1);
    vec_val1 = zeros(1, N+1);
    vec_grad1 = zeros(1, N+1);
    vec_iter1 = zeros(1, N+1);
    %vec_converged1 = zeros(1, N+1);
    conv_order_cell1=cell(1,N+1);
    vec_violations1 = zeros(1, N+1);

    for j = 1:N+1
        tic;
        [x1, f1, gradf_norm1, k1, xseq1, btseq1, cgiterseq1, conv_ord1, flag1, converged1, violations1] = ...
            truncated_newton_27(Mat_points(:,j), F, JF, exact, fin_dif_2, h, kmax, tolgrad, fterms_suplin, cg_maxit, z0, c1, rho, btmax);
        vec_times1(j) = toc;

        disp(['Tentativo n. ', num2str(j), ': ', flag1])
        %vec_converged1(j) = converged1;
        conv_order_cell1{j}=conv_ord1;
        vec_val1(j) = f1;
        vec_grad1(j) = gradf_norm1;
        vec_iter1(j) = k1;
        vec_violations1(j) = violations1;
    end

    % Salva i risultati per questo valore di h
    results_n3(h_idx).h = h;
    results_n3(h_idx).times = vec_times1;
    results_n3(h_idx).values = vec_val1;
    results_n3(h_idx).gradients = vec_grad1;
    results_n3(h_idx).iterations = vec_iter1;
    %results_n3(h_idx).converged = vec_converged1;
    results_n3(h_idx).violations = vec_violations1;
end

for j = 1:length(conv_order_cell1)
    fprintf('Esecuzione %d: ordine di convergenza medio = %f\n', j, mean(conv_order_cell1{j}, 'omitnan'));
   
end

%% n=10^4 (1e4)

rng(345989);

n=1e4; 

kmax=1.5e3;
tolgrad=5e-7;
cg_maxit=100; %50 uguale

z0=zeros(n,1);
c1=1e-4;
rho=0.5;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

x0=(1:n)';  % Initial point
N=10; % number of initial points to be generated

% Initial points:
Mat_points=repmat(x0,1,N+1); 
rand_mat=2*rand(n, N)-1;
Mat_points(:,2:end)=Mat_points(:,2:end) + rand_mat; % matrix with columns=initial points

% Definire i valori di h da testare
results_n4 = struct();  % Struttura per salvare i risultati

for h_idx = 1:length(h_values)
    h = h_values(h_idx);
    fprintf('Test con h = %.4e\n', h);
    
    % Vettori per salvare i risultati per ogni h
    vec_times2 = zeros(1, N+1);
    vec_val2 = zeros(1, N+1);
    vec_grad2 = zeros(1, N+1);
    vec_iter2 = zeros(1, N+1);
    vec_converged2 = zeros(1, N+1);
    conv_order_cell2=cell(1, N+1);
    vec_violations2 = zeros(1, N+1);

    for j = 1:N+1
        tic;
        [x2, f2, gradf_norm2, k2, xseq2, btseq2, cgiterseq2, conv_ord2, flag2, converged2, violations2] = ...
            truncated_newton_27(Mat_points(:,j), F, JF, exact, fin_dif_2, h, kmax, tolgrad, fterms_suplin, cg_maxit, z0, c1, rho, btmax);
        vec_times2(j) = toc;

        disp(['Tentativo n. ', num2str(j), ': ', flag2])
        vec_converged2(j) = converged2;
        vec_val2(j) = f2;
        vec_grad2(j) = gradf_norm2;
        vec_iter2(j) = k2;
        vec_violations2(j) = violations2;
    end

    % Salva i risultati per questo valore di h
    results_n4(h_idx).h = h;
    results_n4(h_idx).times = vec_times2;
    results_n4(h_idx).values = vec_val2;
    results_n4(h_idx).gradients = vec_grad2;
    results_n4(h_idx).iterations = vec_iter2;
    %results_n4(h_idx).converged = vec_converged2;
    results_n4(h_idx).violations = vec_violations2;
end

for j = 1:length(conv_order_cell2)
    fprintf('Esecuzione %d: ordine di convergenza medio = %f\n', j, mean(conv_order_cell2{j}, 'omitnan'));
end


%% n=10^5 (1e5)
rng(345989);

n = 1e5; 
h=[4,6,8,10,12];
h_values =10.^(-h);


kmax = 1.5e3; 
tolgrad = 5e-7; 
cg_maxit = 100; 

z0 = zeros(n,1);
c1 = 1e-4;
rho = 0.5;
btmax = 50; 

x0 = (1:n)';  % Initial point
N = 10; % Number of initial points to be generated

% Initial points:
Mat_points = repmat(x0,1,N+1); 
rand_mat = 2*rand(n, N)-1;
Mat_points(:,2:end) = Mat_points(:,2:end) + rand_mat; % matrix with columns = initial points

results_n5 = struct();  % Struttura per salvare i risultati

for h_idx = 1:length(h_values)
    h = h_values(h_idx);
    fprintf('Test con h = %.4e\n', h);
    
    % Vettori per salvare i risultati per ogni h
    vec_times3 = zeros(1, N+1);
    vec_val3 = zeros(1, N+1);
    vec_grad3 = zeros(1, N+1);
    vec_iter3 = zeros(1, N+1);
    vec_converged3 = zeros(1, N+1);
    vec_violations3 = zeros(1, N+1);

    for j = 1:N+1
        tic;
        [x3, f3, gradf_norm3, k3, xseq3, btseq3, cgiterseq3, conv_ord3, flag3, converged3, violations3] = ...
            truncated_newton_27(Mat_points(:,j), F, JF, exact, fin_dif_2, h, kmax, tolgrad, fterms_suplin, cg_maxit, z0, c1, rho, btmax);
        vec_times3(j) = toc;

        disp(['Tentativo n. ', num2str(j), ': ', flag3])
        vec_converged3(j) = converged3;
        vec_val3(j) = f3;
        vec_grad3(j) = gradf_norm3;
        vec_iter3(j) = k3;
        vec_violations3(j) = violations3;
    end

    % Salva i risultati per questo valore di h
    results_n5(h_idx).h = h;
    results_n5(h_idx).times = vec_times3;
    results_n5(h_idx).values = vec_val3;
    results_n5(h_idx).gradients = vec_grad3;
    results_n5(h_idx).iterations = vec_iter3;
    results_n5(h_idx).converged = vec_converged3;
    results_n5(h_idx).violations = vec_violations3;
end

for j = 1:length(conv_order_cell)
    fprintf('Esecuzione %d: ordine di convergenza medio = %f\n', j, mean(conv_order_cell{j}, 'omitnan'));
end

%%
% Assicurarsi che tutte le variabili siano colonne
h_col = h_values(:); % Converte h_values in colonna
mean_time_col = cellfun(@mean, {results_n4.times})'; % Tempo medio
mean_value_col = cellfun(@mean, {results_n4.values})'; % Valore medio della funzione
mean_iter_col = cellfun(@mean, {results_n4.iterations})'; % Iterazioni medie
mean_conv_ord_col = cellfun(@mean, {results_n4.conv_ord})'; % Ordine medio di convergenza

% Creazione della tabella
T = table(h_col, mean_time_col, mean_value_col, mean_iter_col, mean_conv_ord_col, ...
    'VariableNames', {'h', 'Tempo_Medio', 'Valore_Medio_Funzione', 'Iterazioni_Medie', 'Ordine_Convergenza'});

% Stampa la tabella in console
disp('===== Risultati Medi per h =====');
disp(T);
