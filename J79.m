function J=J79(x)
% Function that computes the gradient of function 79
x_next=[x(2:end);0];
x_prev=[0;x(1:end-1)];

f=(3-x/10).*x+1-x_prev-2*x_next;
f_next=[f(2:end);0];
f_prev=[0;f(1:end-1)];

J=-2*f_prev+f.*(3-x/5)-f_next;

end