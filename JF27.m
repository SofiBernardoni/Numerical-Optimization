function JF=JF27(x,exact,h)
% Function that computes the gradient of function 27
% exact= bool. True= computes the exact version, False= computes the approximated version with finite differences
% h= increment for the approximated version (if exact=true put h=0)

n=length(x);
s=sum(x.^2);
JF=((2*x-2)/100000 + 4*s*x-x)/2;
if ~exact %approximation with finite difference (not exact)
   JF=JF+2*(h^2)*x;  
end

end