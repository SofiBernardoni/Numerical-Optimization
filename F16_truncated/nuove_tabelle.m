%% TABELLE NUOVE
%% VIOLATION

% Metodo Exact (derivate esatte) - media unica
vec_times_ex_clean = vec_violations1_ex; % copia dei tempi
vec_times_ex_clean(vec_converged1_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_f1 = mean(vec_times_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_times_fd1_clean = mat_violations1_fd1;
mat_times_fd1_clean(mat_converged1_fd1 == 0) = NaN;
avg_fd1 = mean(mat_times_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_times_fd2_clean = mat_violations1_fd2;
mat_times_fd2_clean(mat_converged1_fd2 == 0) = NaN;
avg_fd2 = mean(mat_times_fd2_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Creazione delle etichette per i valori di h
h_exponents = [2, 4, 6, 8, 10, 12];  % Solo valori di h (senza h=0)
h_labels = arrayfun(@(e) sprintf('h=1e-%d', e), h_exponents, 'UniformOutput', false);

% Preparazione dei dati per la tabella
% FD1 e FD2 hanno le medie per ogni h, mentre Exact è ripetuto in tutte le colonne
%exact_vals = [avg_exact, avg_exact]; % Esatto in tutte le colonne
fd1_vals = avg_fd1';  % Trasposto per allineare con le colonne
fd2_vals = avg_fd2';  % Trasposto per allineare con le colonne

% Costruzione della tabella
rowNames = {'FD1', 'FD2'};
columnNames = [ h_labels,'Exact']; % Prima colonna "Exact", poi gli h
data = [ fd1_vals, avg_exact_f1; fd2_vals, avg_exact_f1;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T10 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation violation  table (only for successful runs): F16, n=10^3, suplin');
disp(T10);


%% BT-SEQ
% Metodo Exact (derivate esatte) - media unica
vec_bt_ex_clean = vec_bt1_ex; % copia dei tempi
vec_bt_ex_clean(vec_converged1_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_i1 = mean(vec_bt_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_bt_fd1_clean = mat_bt1_fd1;
mat_bt_fd1_clean(mat_converged1_fd1 == 0) = NaN;
avg_fd1 = mean(mat_bt_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_bt_fd2_clean = mat_bt1_fd2;
mat_bt_fd2_clean(mat_converged1_fd2 == 0) = NaN;
avg_fd2 = mean(mat_bt_fd2_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Creazione delle etichette per i valori di h
h_exponents = [2, 4, 6, 8, 10, 12];  % Solo valori di h (senza h=0)
h_labels = arrayfun(@(e) sprintf('h=1e-%d', e), h_exponents, 'UniformOutput', false);

% Preparazione dei dati per la tabella
% FD1 e FD2 hanno le medie per ogni h, mentre Exact è ripetuto in tutte le colonne
%exact_vals = [avg_exact, avg_exact]; % Esatto in tutte le colonne
fd1_vals = avg_fd1';  % Trasposto per allineare con le colonne
fd2_vals = avg_fd2';  % Trasposto per allineare con le colonne

% Costruzione della tabella
rowNames = {'FD1', 'FD2'};
columnNames = [ h_labels,'Exact']; % Prima colonna "Exact", poi gli h
data = [ fd1_vals, avg_exact_i1; fd2_vals, avg_exact_i1;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T11 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation bt iteration table (only for successful runs): F16, n=10^3, suplin');
disp(T11);



%% CG-SEQ
% Metodo Exact (derivate esatte) - media unica
vec_bt_ex_clean = vec_cg_iter1_ex; % copia dei tempi
vec_bt_ex_clean(vec_converged1_ex == 0) = NaN; % sostituisce con NaN i non convergenti
avg_exact_i1 = mean(vec_bt_ex_clean, 'omitnan');  % calcola la media (scalare)

% Metodo FD1 (Finite Differences classiche)
mat_bt_fd1_clean = mat_cg_iter1_fd1;
mat_bt_fd1_clean(mat_converged1_fd1 == 0) = NaN;
avg_fd1 = mean(mat_bt_fd1_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Metodo FD2 (Finite Differences nuove)
mat_bt_fd2_clean = mat_cg_iter1_fd2;
mat_bt_fd2_clean(mat_converged1_fd2 == 0) = NaN;
avg_fd2 = mean(mat_bt_fd2_clean, 2, 'omitnan'); % media per ogni h (6x1)

% Creazione delle etichette per i valori di h
h_exponents = [2, 4, 6, 8, 10, 12];  % Solo valori di h (senza h=0)
h_labels = arrayfun(@(e) sprintf('h=1e-%d', e), h_exponents, 'UniformOutput', false);

% Preparazione dei dati per la tabella
% FD1 e FD2 hanno le medie per ogni h, mentre Exact è ripetuto in tutte le colonne
%exact_vals = [avg_exact, avg_exact]; % Esatto in tutte le colonne
fd1_vals = avg_fd1';  % Trasposto per allineare con le colonne
fd2_vals = avg_fd2';  % Trasposto per allineare con le colonne

% Costruzione della tabella
rowNames = {'FD1', 'FD2'};
columnNames = [ h_labels,'Exact']; % Prima colonna "Exact", poi gli h
data = [ fd1_vals, avg_exact_i1; fd2_vals, avg_exact_i1;];

% Creiamo la tabella con i nomi delle colonne e delle righe
T12 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Average computation cg iteration table (only for successful runs): F16, n=10^3, suplin');
disp(T12);



%% Calcolo quanti a convergenza


% Creazione delle etichette per i valori di h
h_exponents = [2, 4, 6, 8, 10, 12];  % Solo valori di h (senza h=0)
h_labels = arrayfun(@(e) sprintf('h=1e-%d', e), h_exponents, 'UniformOutput', false);

% Preparazione dei dati per la tabella
% FD1 e FD2 hanno le medie per ogni h, mentre Exact è ripetuto in tutte le colonne
%exact_vals = [avg_exact, avg_exact]; % Esatto in tutte le colonne
fd1_vals = sum(mat_converged1_fd1,2)';  % Trasposto per allineare con le colonne
fd2_vals = sum(mat_converged1_fd2,2)';  % Trasposto per allineare con le colonne

% Costruzione della tabella
rowNames = {'FD1', 'FD2'};
columnNames = [ h_labels,'Exact']; % Prima colonna "Exact", poi gli h
data = [ fd1_vals, sum(vec_converged1_ex); fd2_vals, sum(vec_converged1_ex);];

% Creiamo la tabella con i nomi delle colonne e delle righe
T13 = array2table(data, 'VariableNames', columnNames, 'RowNames', rowNames);

% Visualizza la tabella
disp('Number of converged : F16, n=10^3, suplin');
disp(T13);

writetable(T1, 'results_f16_suplin.xlsx', 'Sheet', 'time_3','WriteRowNames', true);
writetable(T2, 'results_f16_suplin.xlsx', 'Sheet', 'niter_3','WriteRowNames', true);
writetable(T3, 'results_f16_suplin.xlsx', 'Sheet', 'f_val_3','WriteRowNames', true);
writetable(T11, 'results_f16_suplin.xlsx', 'Sheet', 'bt_3','WriteRowNames', true);
writetable(T12, 'results_f16_suplin.xlsx', 'Sheet', 'cg_3','WriteRowNames', true);
writetable(T13, 'results_f16_suplin.xlsx', 'Sheet', 'n_conv3','WriteRowNames', true);
