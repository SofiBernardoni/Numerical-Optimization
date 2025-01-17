function HF=HF27(x)
n=length(x);
s=sum(x.^2);
HF=4*x*x';
HF(1:n+1:end)=(2/100000 -1+ 4*s + 8*x.^2)/2;

end


