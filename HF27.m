function HF=HF27(x,sparse)
% Function that computes the Hessian of function 27
% sparse= bool. True= computes the sparse version

n=length(x);
s=sum(x.^2);
HF=4*x*x';
HF(1:n+1:end)=(2/100000 -1+ 4*s + 8*x.^2)/2;
if sparse
    HF=sparse(HF); %  non risolve il problema CAMBIA
end

end
