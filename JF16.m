function JF=JF16(x)
n=length(x);
JF=zeros(length(x),1);
indices=(1:n)';
JF=indices.*sin(x)-2*cos(x);
JF(n)=(n-1)*cos(x(n))+n*sin(x(n));