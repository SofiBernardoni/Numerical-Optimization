function JF=JF16(x,exact,fin_dif_2,h)
% Function that computes the gradient of function 27
% exact= bool. True= computes the exact version, False= computes the approximated version with finite differences
% fin_dif_2= bool. True if exact=false and finite differences are computed using as increment h*abs(x_j) for the derivative with respect to j
% h= increment for the approximated version (if exact=true put h=0)

n=length(x);

indices=(1:n)';
JF=indices.*sin(x)+2*cos(x);
JF(n)=(1-n)*cos(x(n))+n*sin(x(n));
if ~exact %approximation with finite difference (not exact)
   if fin_dif_2 % version of finite differences with abs(xj)
       JF=(sin(h*abs(x))./(h*abs(x))).*JF;
   else % classic version of finite differences
       JF=sin(h)/h*JF;
   end
end

end