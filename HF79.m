function HF=HF79(x,sparse)
% Function that computes the Hessian of function 79
% sparse= bool. True= computes the sparse version

x_next=[x(2:end);0];
x_prev=[0;x(1:end-1)];
n= length(x);

f=(3-x/10).*x+1-x_prev-2*x_next;
if sparse
    diag0=5-f/5+(3-x/5).^2;
    diag_1= -2*(3-x(1:end-1)/5)-(3-x_next(1:end-1)/5);
    diag_up=[0;diag_1] ;
    diag_down=[diag_1 ; 0];
    HF=spdiags([diag_down, diag0,diag_up],[-1,0,1],n,n); % ci saranno problemi per dimensioni delle due sovradiagonali?? (in teorai ora no)
else
    diag0=5-f/5+(3-x/5).^2;
    diag_up= -2*(3-x(1:end-1)/5)-(3-x_next(1:end-1)/5);
    HF=diag(diag0)+diag(diag_up,1)+diag(diag_up,-1);
end

end