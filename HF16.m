function HF=HF16(x,sparse)
% Function that computes the Hessian of function 79
% sparse= bool. True= computes the sparse version

n=length(x);
indices=(1:n)';
if sparse
    D=indices.*cos(x)-2*sin(x);
    HF=spdiags(D,0,n,n);
    HF(n,n)=(n-1)*sin(x(n))+n*cos(x(n));
else
    HF=diag(indices.*cos(x)-2*sin(x));
    HF(n,n)=(n-1)*sin(x(n))+n*cos(x(n));
end

end