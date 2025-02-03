function [xk, fk, gradfk_norm, k, xseq, btseq,cgiterseq,convergence_order,flag, converged, violations] = truncated_newton_27(x0, f, gradf,exact, fin_dif_2, h, kmax, tolgrad, ftol, cg_maxit,z0, c1, rho, btmax)
% Function that performs the truncated Newton optimization method, for for function F27, with backtracking. 
% INPUTS:
% x0 = n-dimensional column vector. Initial point;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f;
% exact = bool. True if exact version of the hessian. False= approximated version with finite differences
% h= increment for finite differences. if exact=true put h=0.
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
    % Hessf27(x)(i,j)= 4*x_i*x_j 
    % Hessf27(x)(i,i)= (2/100000 -1+ 4*(sum(x.^2)) + 8*x_i^2)/2
    % the matrix is NOT sparse BUT with large n cannot be stored. So we
    % compute directly the matrix vector products.
    % EXACT VERSION: Hessf27(x)*z= 4*s*v1 -4*v2 +v3   with s=sum(x), v1=x.*z, v2= (x.^2).*z, v3=diag(Hessf27(x)).*z
    % diag(Hessf27(x))= (2/100000 -1+ 4*s + 8*x.^2)/2;
    % APPROXIMATED VERSION:  Hessf27_approx(x,h)*z= Hessf27(x)*z + 2*n*(h^2)*z


    %with finite difference 1: Hessf27(x)(i,j)= 4*x_i*x_j + h^2 +2hx_j+2hx_i
    %APPROXIMATED VERSION: Hessf27(x)(i,i)= (2/100000 -1+ 4*(sum(x.^2)) + 8*x_i^2 + 2*h^2)/2
    %APPROXIMATED VERSION:  Hessf27_approx(x,h)*z= Hessf27(x)*z +(h^2)*sum_z*ones(n,1) - (h^2)*z + 2*h*sum_z*x + 2*h*(x'*z)-4*h*(x.*z)

    %with finite difference 2: Hessf27(x)(i,j)= 4*x_i*x_j*mod(x_i)*mod(x_j)
    %+ h^2*mod(x_i)^2*mod(x_j)^2
    %+2hx_j*mod(x_i)^2*mod(x_j)+2hx_i*mod(x_j)^2*mod(x_i)

    diagA=(2/100000 -1+ 4*sum(xk.^2) + 8*xk.^2)/2;
    
    % Initialization of zj and j
    zj = z0; 
    j= 0; 
    
    % Initialization of relative residual and of descent direction
    Azj= 4*(xk'*zj)*xk-4*(xk.^2).*zj+diagA.*zj; % A*zj
    sum_z=sum(zj);
    if ~exact %approximation with finite difference (not exact)
        if fin_dif_2
            Azj= (diagA.*zj +(h^2.*abs(zj).^2)) +( h^2*abs(zj).^2*((abs(xk).^2)'*zj) - h^2*abs(zj).^4.*zj )+(2*h*xk.*abs(xk)*((abs(xk).^2)'*zj) - 2*h*xk.*(abs(xk).^2).*zj) +(2*h*abs(xk).^2*sum(xk.*abs(xk).*zj)-2*h*abs(xk).^3.*xk.*zj ) + (4*xk.*abs(xk)*sum(xk.*abs(xk).*zj)-4*abs(xk).^2.*xk.^2.*zj);
        else
            %Azj= Azj+2*(h^2)*sum(zj)*ones(n,1);
            Azj=Azj+ (h^2)*sum_z*ones(length(x0),1) - (h^2)*zj + 2*h*sum_z*xk + 2*h*(xk'*zj)-4*h*(xk.*zj);
        end
    end



    res = -gradk - Azj; % initialize relative residual res=b-Ax 
    p = res; % initialize descent direction
    norm_b = gradfk_norm; % norm(b) where b=-gradk
    norm_r = norm(res); % norm of the residual
    
    %neg_curv= false; % boolean checking negative curvature condition
    
    while (j<cg_maxit && norm_r>ftol(j,norm_b)*norm_b ) %adaptive tolerance based on the norm of the gradient
       z= 4*(xk'*p)*xk-4*(xk.^2).*p+diagA.*p; % A*p : product of A and descent direction
       sum_p=sum(p);
       if ~exact %approximation with finite difference (not exact)
           %z= z+2*(h^2)*sum(p)*ones(n,1);
           if fin_dif_2
                 z= (diagA.*p +(h^2.*abs(p).^2)) +( h^2*abs(p).^2*((abs(xk).^2)'*p) - h^2*abs(p).^4.*p )+(2*h*xk.*abs(xk)*((abs(xk).^2)'*p) - 2*h*xk.*(abs(xk).^2).*p) +(2*h*abs(xk).^2*sum(xk.*abs(xk).*p)-2*h*abs(xk).^3.*xk.*p ) + (4*xk.*abs(xk)*sum(xk.*abs(xk).*p)-4*abs(xk).^2.*xk.^2.*p);
           else
               %%%%%% Da controllare se ci va z o zj
                z=z+ (h^2)*sum_p*ones(length(x0),1) - (h^2)*p + 2*h*sum_p*xk + 2*h*(xk'*p)-4*h*(xk.*p);
           end
       end
       a = (res'*p)/(p'*z); % update exact step along the direction
       zj = zj+ a*p; % update solution 
       res = res - a*z; %update residual
       beta = -(res'*z)/(p'*z);
       p = res + beta*p; % update descent direction
       
       % se vuoi sposta qui calcolo z per usarlo in condizione %%%%%%%%%%%%%%%%%%%%%%%%%
       z_new=4*(xk'*p)*xk'-4*((xk.^2).*p)'+(diagA.*p)'; % p'*A (as A*p because A symmetric but as a row vector)  --> needed for curvature condition
       sum_p=sum(p);
       if ~exact %approximation with finite difference (not exact)
            if fin_dif_2  
                z_new= ((diagA.*p +(h^2.*abs(p).^2)) +( h^2*abs(p).^2*((abs(xk).^2)'*p) - h^2*abs(p).^4.*p )+(2*h*xk.*abs(xk)*((abs(xk).^2)'*p) - 2*h*xk.*(abs(xk).^2).*p) +(2*h*abs(xk).^2*sum(xk.*abs(xk).*p)-2*h*abs(xk).^3.*xk.*p ) + (4*xk.*abs(xk)*sum(xk.*abs(xk).*p)-4*abs(xk).^2.*xk.^2.*p))';
            else
                z_new=z_new+ ((h^2)*sum_p*ones(length(x0),1))' - ((h^2)*p+ 2*h*sum_p*xk)' + (2*h*(xk'*p)-4*h*(xk.*p))';
            end
       end
       sign_curve=sign(z_new*p); %%%%%%%%%%%%%% CONTROLLA (forse rimettere A*p)
       if sign_curve ~= 1 % negative curvature condition  p'*A*p <= 0
           violations =violations+1; %%%%%%%% togli
           break;
       end

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