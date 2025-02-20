function [x, f_k, n_iter]=Nelder_mead(x0,f,rho,mu, gamma, sigma, tol, max_iter, Delta)
n=length(x0);
% NOTATION
%n = dimension of vectors
% S = Simplex = matrix n x n+1 --> n+1 vectors of lenght n
% f_val_S = row vector, lenght=n+1 containing the values of the function in
% the n+1 point of the simplex

% Function that performs the Nelder-Mead optimization method, for a
% given function f

% INPUTS:
% x0 = n-dimensional column vector. Initial point;
% f = function handle that describes a function R^n->R;
% rho = reflection parameter
% mu= expansion parameter
% gamma = contraction parameter
% sigma = shrinking parameter
% tol= tolerance for stopping criteria
% max_iter= maximum number of iteration permitted;

% OUTPUTS:
% xk = the sequence of the best xk ad every iteration;
% fk = the sequence of f(xk) ad every iteration value;
% n_iter = number of iteration 

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

[f_val_S, sort_idx] = sort(f_val_S); %vector with ordered index eith respect to the values of the function
idx=idx(sort_idx); 

while iter < max_iter && abs(f_val_S(n+1) - f_val_S(1)) > tol

    %---------REFLECTION PHASE ----------------------
    %computation of barycenter point
    x_bar=mean(S(:, idx(1:n)),2); % 

    %computation and evaluation of reflection point
    x_r= x_bar + rho*(x_bar - S(:,idx(n+1)));
    f_r=f(x_r);

    if f_r < f_val_S(1)
        %-----------EXPANSION-----------
        %computation and evaluation of expansion point
        x_e=x_bar + mu*(x_r-x_bar);
        f_e=f(x_e);

        if f_e < f_r
            %hold expansion point
            S(:,idx(n+1))=x_e;
            f_val_S(n+1)=f_e;
        else
            % hold reflection point
            S(:, idx(n+1)) = x_r;
            f_val_S(n+1) = f_r;
        end

    elseif f_r < f_val_S(n)

        %hold reflexion point x_r
        S(:, idx(n+1)) = x_r;
        f_val_S(n+1) = f_r;
        
    else
        %-------------------CONTRACTION------------
        %calculate x_c with the best point between x_r, x_n+1
        if f_r < f_val_S(n+1)
            x_c= x_bar + gamma *(x_bar - x_r);
        else 
            x_c= x_bar + gamma *(x_bar - S(:,idx(n+1)));
        end
        f_c=f(x_c);
        if f_c < f_val_S(n+1)
            %hold x_c 
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
    [f_val_S, sort_idx] = sort(f_val_S);  
    idx=idx(sort_idx); 
    
    iter = iter + 1;
    
    x=[x, S(:,idx(1))];
    f_k=[f_k;f_val_S(1)];
end

n_iter=iter;

end


