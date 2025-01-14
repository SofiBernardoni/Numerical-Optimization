function F=F79(x)
term1=(3-x/10).*x;
x_next=[x(2:end);0];
x_prev=[0;x(1:end-1)];
term2=1-x_prev-2*x_next;
F=sum((term1+term2).^2);
F=F/2;
end