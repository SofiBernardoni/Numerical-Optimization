function HF=HF79(x,sparse,exact,h)
% Function that computes the Hessian of function 79
% sparse= bool. True= computes the sparse version
% exact= bool. True= computes the exact version, False= computes the approximated version with finite differences
% h= increment for the approximated version (if exact=true put h=0)

x_next=[x(2:end);0];
x_prev=[0;x(1:end-1)];
n= length(x);

f=(3-x/10).*x+1-x_prev-2*x_next;
diag0=5-f/5+(3-x/5).^2;
diag_1= -2*(3-x(1:end-1)/5)-(3-x_next(1:end-1)/5);
if sparse % sparse
    diag_up=[0;diag_1] ;
    diag_down=[diag_1 ; 0];
    HF=spdiags([diag_down, diag0,diag_up],[-1,0,1],n,n); %exact version
    if ~exact  %approximation with finite difference (not exact)
        HF= HF + spdiags((h^2)*ones(n,1)/100,0,n,n);
        new_coef=3*h/10;
        HF= HF + spdiags(new_coef*ones(n-1,1),1,n,n) +spdiags(new_coef*ones(n-1,1),-1,n,n);
    end
    
else % NOT sparse
    HF=diag(diag0)+diag(diag_1,1)+diag(diag_1,-1); %exact version
    if ~exact  %approximation with finite difference (not exact)
        HF= HF + (h^2)*eye(n)/100 ;
        new_coef=3*h/10;
        HF= HF + diag(new_coef*ones(n-1,1),1) +diag(new_coef*ones(n-1,1),-1);
    end
end

end