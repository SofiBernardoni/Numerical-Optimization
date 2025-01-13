function JF=JF27(x)
n=length(x);
s=sum(x.^2);
JF=((2*x-2)/100000 + 4*x*s-x)/2;