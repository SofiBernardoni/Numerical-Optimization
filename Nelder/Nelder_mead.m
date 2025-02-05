function [x, f_k, n_iter]=Nelder_mead(x0,f,rho,mu, gamma, sigma, tol, max_iter, Delta)
n=length(x0);
%NOTAZIONE
%n=dimesione vettori
% S=simplesso = matrice n*n+1= n+1 vettori colonna di lunghezza n
% f_val_S= vettore riga di lunghezza n+1 con i valori assunti dalla
% funzione nei n+1 del simplesso

x=x0;
f_k=f(x0);

S = [x0, x0 + Delta * eye(n)];


%f_values in the vertices of the simplex
f_val_S=zeros(1,n+1);
for i = 1:n+1
   f_val_S(i) = f(S(:,i));
end


idx=1:n+1;
iter=0;

[f_val_S, sort_idx] = sort(f_val_S); %vettore con gli indici ordinati per valore della funzione 
idx=idx(sort_idx); % sono le posizioni dei punti nel simplesso

while iter < max_iter
    %------------ordiniamo i punti in ordine crescente:

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
            %disp(['numero iterazione ',num2str(iter), 'espansione migliore 1'])
        else
            % rimani con reflecion
            S(:, idx(n+1)) = x_r;
            f_val_S(n+1) = f_r;
            %disp(['numero iterazione ',num2str(iter), 'riflesso migliore 1'])
            
        end

        [f_val_S, sort_idx] = sort(f_val_S); %vettore con gli indici ordinati per valore della funzione 
        idx=idx(sort_idx); % sono le posizioni dei punti nel simplesso

    elseif f_r < f_val_S(n)
        %se sta in mezzo accetti il punto x_r
        S(:, idx(n+1)) = x_r;
        f_val_S(n+1) = f_r;
        %disp(['numero iterazione ',num2str(iter), 'riflesso  tra 1 e n'])

        [f_val_S, sort_idx] = sort(f_val_S); %vettore con gli indici ordinati per valore della funzione 
        idx=idx(sort_idx); % sono le posizioni dei punti nel simplesso

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
            %disp(['numero iterazione ',num2str(iter), 'contraction'])
        else
            %----------------------SHRINKAGE------------
            %disp(['numero iterazione ',num2str(iter), 'shrink'])
            for i = 2:n+1
                    S(:, idx(i)) = S(:, idx(1)) + sigma * (S(:, idx(i)) - S(:,idx(1) ));
                    f_val_S(i) = f(S(:, idx(i)));
            end
            % controllo su grandezza simplesso

        end
        [f_val_S, sort_idx] = sort(f_val_S); %vettore con gli indici ordinati per valore della funzione 
        idx=idx(sort_idx); % sono le posizioni dei punti nel simplesso

    end

    % ---------------- CRITERI DI ARRESTO ----------------
    % 1. Differenza nei valori della funzione obiettivo
    if max(abs(f_val_S(n+1) - f_val_S(1))) < tol
        %disp('esco per f_val')
        break;
    end

    
    iter = iter + 1;
    %%%%% errore qua
    x=[x, S(:,idx(1))];
    f_k=[f_k;f_val_S(1)];
end

n_iter=iter;

end


