function JF=JF16(x,exact,h)
% Function that computes the gradient of function 27
% exact= bool. True= computes the exact version, False= computes the approximated version with finite differences
% h= increment for the approximated version (if exact=true put h=0)

n=length(x);
%JF=zeros(length(x),1);
indices=(1:n)';
JF=indices.*sin(x)+2*cos(x);
JF(n)=(1-n)*cos(x(n))+n*sin(x(n));
if ~exact %approximation with finite difference (not exact)
   JF=sin(h)/h*JF;  
end

end