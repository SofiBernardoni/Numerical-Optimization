%% PROVA PRECONDITIONING:
%%% Compute pk solving Hessf(xk)pk=-gradk with Coniugate Gradient method. %%%
    % Hessf(xk)=A, pk=z, -gradk=b
    A=Hessf(xk); % computing Hessian (if A sparse products with dense vectors will be dense)
    % Initialization of zj and j
    zj = z0; 
    j= 0; 
    % Inizializzazione del residuo relativo e della direzione di discesa
    res = A*zj-gradk ; % initialize relative residual res=b-Ax 

    % M da definire:
    % ci sono 2 opzioni:
    % M=ichol(A); % se è simmetrica e definita positiva
    % M=(D+L); % D diag L triangolarew iferiore

    D = diag(diag(A));  % Matrice diagonale (D)
    L = tril(A, -1);    % Matrice inferiore (L)
    M=D+L;
    M_inv=inv(M);
    y=pcg(M,res);
    p = -y; % initialize descent direction


    norm_b = gradfk_norm; % norm(b) where b=-gradk
    norm_r = norm(res); % norm of the residual
    
    %neg_curv= false; % boolean checking negative curvature condition
    
    while (j<cg_maxit && norm_r>ftol(j,norm_b)*norm_b ) %adaptive tolerance based on the norm of the gradient
        z = A*p; %product of A and descent direction 
       a = (res'*y)/(p'*z); % update exact step along the direction
       zj = zj+ a*p; % update solution 
       res1 = res - a*z; %update residual

       %risolvere il sistema Myk+1=rk+1
       %y1=M\res1

       beta = (res1'*y)/(res'*y);
       p = -y1 + beta*p; % update descent direction
       
       % se vuoi sposta qui calcolo z per usarlo in condizione %%%%%%%%%%%%%%%%%%%%%%%%%
       sign_curve=sign(res'*M_inv*res1);
       if sign_curve ~= 1 % negative curvature condition  p'*A*p <= 0
           violations =violations+1; %%%%%%%% togli
           break;
       end

       res=res1;
       y=y1;

       norm_r = norm(res);
       j = j+1;
    
    end
    