function F1=F16(x,n)
term1=1-cos(x);
sin_prec=[0; x(1:end-1)];
sin_next=[x(2:end); 0];
term2=sin(sin_prec)-sin(sin_next);
indici=(1:n)';
F1=sum(indici.*(term1+term2));
end