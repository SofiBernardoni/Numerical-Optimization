function HF=HF16(x)
n=length(x);
indices=(1:n)';
HF=diag(indices.*cos(x)-2*sin(x));
HF(n,n)=(n-1)*sin(x(n))+n*cos(x(n));