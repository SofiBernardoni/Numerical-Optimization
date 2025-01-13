function F=F27(x)
fk=(x-1)/sqrt(100000);
F=sum(fk.^2);
F=F+(sum(x.^2)-0.25)^2;
F=F/2;