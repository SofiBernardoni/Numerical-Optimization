function JF=JF79(x,exact,h)
% Function that computes the gradient of function 79
% exact= bool. True= computes the exact version, False= computes the approximated version with finite differences
% h= increment for the approximated version (if exact=true put h=0)

x_next=[x(2:end);0];
x_prev=[0;x(1:end-1)];

f=(3-x/10).*x+1-x_prev-2*x_next;
f_next=[f(2:end);0];
f_prev=[0;f(1:end-1)];

JF=-2*f_prev+f.*(3-x/5)-f_next;
if ~exact %approximation with finite difference (not exact)
   JF=JF-(3-x/5)*(h^2)/10;  
end

end