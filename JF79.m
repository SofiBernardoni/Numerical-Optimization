function JF=JF79(x,exact,fin_dif_2,h)
% Function that computes the gradient of function 79
% exact= bool. True= computes the exact version, False= computes the approximated version with finite differences
% fin_dif_2= bool. True if exact=false and finite differences are computed using as increment h*abs(x_j) for the derivative with respect to j
% h= increment for the approximated version (if exact=true put h=0)

x_next=[x(2:end);0];
x_prev=[0;x(1:end-1)];

f=(3-x/10).*x+1-x_prev-2*x_next;
f_next=[f(2:end);0];
f_prev=[0;f(1:end-1)];

JF=-2*f_prev+f.*(3-x/5)-f_next;
if ~exact %approximation with finite difference (not exact)
   if fin_dif_2 % version of finite differences with abs(xj)
       JF=JF-(3-x/5).*(x.^2)*(h^2)/10;
   else % classic version of finite differences
       JF=JF-(3-x/5)*(h^2)/10;
   end
     
end

end