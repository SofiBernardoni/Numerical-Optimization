function H=HF79(x)
% Function that computes the Hessian of function 79

x_next=[x(2:end);0];
x_prev=[0;x(1:end-1)];

f=(3-x/10).*x+1-x_prev-2*x_next;
diag0=5-f/5+(3-x/5).^2;
diag_up= -2*(3-x/5)-(3-x_next/5);
H=diag(diag0)+diag(diag_up,1);

end