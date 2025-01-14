function [xk, fk, n_iter]=Nelder_mead(x0,f,rho,mu, gamma, sigma, tol, max_iter)
n=length(x0);
%NOTAZIONE
%n=dimesione vettori
% S=simplesso = matrice n*n+1= n+1 vettori colonna di lunghezza n
% f_val_S= vettore riga di lunghezza n+1 con i valori assunti dalla
% funzione nei n+1 del simplesso


%create the initial simplex:
%S=repmat(x0,1,n+1);
%rand_mat= 2 * (rand([n, n+1]) - 0.5);
%S=S+rand_mat;

Delta = 1;  % Perturbazione: variazione di ogni componente
S = [x0, x0 + Delta * eye(n)];


%f_values in the vertices of the simplex
f_val_S=zeros(1,n+1);
for i = 1:n+1
    f_val_S(i) = f(S(:,i));
end
idx=1:n+1;
iter=0;


while iter < max_iter
    %------------ordiniamo i punti in ordine crescente:
    [f_val_S, sort_idx] = sort(f_val_S); %vettore con gli indici ordinati per valore della funzione 
    idx=idx(sort_idx);
    %risistemiamo il simplesso:
    %S=S(:, idx);

    %---------REFLECTION PHASE ----------------------
    %calcolo del baricentro
    x_bar=mean(S(:, idx(1:n)),2); % 

    %calcolo riflesso
    x_r= x_bar + rho*(x_bar - S(:,idx(n+1)));
    f_r=f(x_r);

    if f_r < f_val_S(1)
        %-----------EXPANSION-----------
        x_e=x_bar + mu*(x_r-x_bar);
        f_e=f(x_e);
        if f_e < f_r
            S(:,idx(n+1))=x_e;
            f_val_S(n+1)=f_e;
        else
            % rimani con reflecion
            S(:, idx(n+1)) = x_r;
            f_val_S(n+1) = f_r;
        end
    elseif f_r < f_val_S(n)
        %se sta in mezzo accetti il punto x_r
        S(:, idx(n+1)) = x_r;
        f_val_S(n+1) = f_r;
    else
        %-------------------CONTRACTION------------
        if f_r < f_val_S(n+1)
            x_c= x_bar + gamma *(x_bar - x_r);
        else 
            x_c= x_bar + gamma *(x_bar - S(:,idx(n+1)));
        end
        f_c=f(x_c);
        if f_c < f_val_S(n+1)
            S(:, idx(n+1)) = x_c;
            f_val_S(n+1) = f_c;
        else
            %----------------------SHRINKAGE------------
            for i = 2:n+1
                    S(:, idx(i)) = S(:, idx(1)) + sigma * (S(:, idx(i)) - S(:,idx(1) ));
                    f_val_S(i) = f(S(:, idx(i)));
            end
        end
    end

    % ---------------- CRITERI DI ARRESTO ----------------
    % 1. Differenza nei valori della funzione obiettivo
    if max(abs(f_val_S - f_val_S(1))) < tol
        disp('esco per f_val')
        break;
    end

    % 2. Diametro del simplesso: troppo costoso
    %max_dist = 0;
    %for i = 1:n+1
        %for j = i+1:n+1
           % max_dist = max(max_dist, norm(S(:, idx(i)) - S(:, idx(j))));
        %end
    %end
    %if max_dist < tol
        %break;
    %end

    %distanza dal punto ottimale:
    %if max(abs(x_min - S(:,idx(1)))) < tol
        %break;
    %end

    iter = iter + 1;
end





xk=S(:,idx(1));
fk=f_val_S(1);
n_iter=iter;

end


