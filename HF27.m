function HF=HF27(x)
n=length(x);
s=sum(x.^2);
HF=diag((2/100000 -1+ 4*s + 8*x.^2)/2);
for i=1:n
    for j=i:n
        HF(i,j)=4*x(i)*x(j);
        HF(j,i)=4*x(i)*x(j);
    end
end
