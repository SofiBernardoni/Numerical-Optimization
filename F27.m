function F=F27(x,n)
f_k=@(x,k) (x(k)-1)/sqrt(100000) ;
F=0;
for j=1:n
    F=F+f_k(x,j)^2;
end
F=F+ (sum(x.^2)-0.25)^2;