function F=F16(n,x)
F=0;
x=[0;x; 0];
for i = 2:n
    F = F + (i * (  (1-cos(x(i)))  +sin(x(i-1))  - sin(x(i+1))    ) );
end 
end