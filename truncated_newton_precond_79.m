function [xk, fk, gradfk_norm, k, xseq, btseq,cgiterseq,convergence_order,flag, converged, violations] = truncated_newton_precond_79(x0, f, gradf, Hessf, kmax, tolgrad, ftol, cg_maxit,z0, c1, rho, btmax)
% Function that performs the truncated Newton optimization method, for a
% given function f, with backtracking. This version can use both the exact derivatives and the approximated version.
%
% INPUTS:
% x0 = n-dimensional column vector. Initial point;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f;
% Hessf= function handle that describes the Hessian of f;
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the gradient
% ftol= function handle of relative tolerance depending on the norm of the gradient (for coniugate gradient method)
% cg_maxit = maximum number of iterations of coniugate gradient method
% c1= factor for the Armijo condition in (0,1);
% rho= fixed (for simplicity) factor less than 1 used to reduce alpha in
% backtracking;
% btmax= maximum number of backtracks permitted;
%%%%%%%%%%%%%%%%% AGGIUNGERE: z0 (valore iniziale per risoluzione sistema)
%
% OUTPUTS:
% xk = the last x computed by the function;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf(xk)
% k = index of the last iteration performed
% xseq = n-by-k matrix where the columns are the elements xk of the sequence
% btseq = row vector with the number of backtracks done at every iteration
% cgiterseq=
% convergence_order = estimated order of convergence
% flag= string that says how the method has ended
% converged= bool. True if the method has converged
% violations=number of violations of positive curvature condition


% Initializations
xseq = zeros(length(x0), kmax);
cgiterseq = zeros(1, kmax); 
btseq = zeros(1,kmax);
convergence_order=zeros(1,kmax);

xk = x0; % assigning the initial point
k = 0;
gradk= gradf(xk); % assigning the initial gradient of f(xk)
gradfk_norm = norm(gradk); % assigning the initial gradient norm
flag=nan;

violations=0; %%%%%%%%%% togli

%%%%%%%%%% armijo function
while k < kmax && gradfk_norm > tolgrad
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
    %M=ichol(A); % se è simmetrica e definita positiva
    % M=(D+L); % D diag L triangolarew iferiore

    D = diag(diag(A));  % Matrice diagonale (D)
    L = tril(A, -1);    % Matrice inferiore (L)
    M=D+L;
    %M_inv=inv(M);
    y=M\res;
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
       y1=M\res1;
       %disp('sitema risolto');

       beta = (res1'*y)/(res'*y);
       p = -y1 + beta*p; % update descent direction
       
       % se vuoi sposta qui calcolo z per usarlo in condizione %%%%%%%%%%%%%%%%%%%%%%%%%
       sign_curve=sign(p'*A*p);
       if sign_curve ~= 1 % negative curvature condition  p'*A*p <= 0
           violations =violations+1; %%%%%%%% togli
           break;
       end

       res=res1;
       y=y1;

       norm_r = norm(res);
       j = j+1;
    
    end
    
    pk=zj; % descent direction computed (considering the negative curvature condition)

    
    % Backtracking to compute the steplength
    bt=0;
    alpha=1; % initial steplenght=1
    xnew = xk + alpha * pk; % Compute the new value for x with alpha
    while bt<btmax && (f(xnew)>(f(xk)+c1*alpha*(gradk'*pk))) % Armijo condition
        alpha=rho*alpha;
        xnew = xk + alpha * pk; % Compute the new value for x with alpha
        bt = bt +1;
    end
     
    if bt==btmax && f(xnew)>(f(xk)+c1*alpha*(gradk'*pk))  % Break if armijo not satisfied
        flag='Procedure stopped because the Armijo condition was NOT satisfied';
        converged=false;
        break;
    else 
        xk=xnew;
    end

    gradk= gradf(xk); % assigning the initial gradient of f(xk)
    gradfk_norm = norm(gradk); % assigning the initial gradient norm
    k = k + 1; % Increase the step by one
    
    xseq(:, k) = xk; % Store current xk in xseq
    btseq(k)=bt; % Store number of backtracking iterations
    cgiterseq(k)=j; % Store coniugate gradient iterations in pcgiterseq
    if k>3
        convergence_order(k)=log(norm(xseq(:, k)-xseq(:, k-1))/norm(xseq(:, k-1)-xseq(:, k-2)))/log(norm(xseq(:, k-1)-xseq(:, k-2))/norm(xseq(:, k-2)-xseq(:, k-3)));
    end
end

if isnan(flag)
    if k==kmax && gradfk_norm > tolgrad 
        flag='Procedure stopped because the maximum number of iterations was reached';
        converged=false;
    else
        flag=['Procedure stopped in ',num2str(k), ' steps, with gradient norm ', num2str(gradfk_norm)];
        converged=true;
    end
end
fk = f(xk); % Compute f(xk)

xseq = xseq(:, 1:k); % "Cut" xseq to the correct size
xseq = [x0, xseq]; % "Add" x0 at the beginning of xseq (otherwise the first el. is x1)
btseq = btseq(1:k); % "Cut" btseq to the correct size
cgiterseq = cgiterseq(1:k); % "Cut" cgiterseq to the correct size
convergence_order=convergence_order(1:k); % "Cut" convergence order


end