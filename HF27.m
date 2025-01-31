function HF=HF27(x,sparse,exact,h)
% Function that computes the Hessian of function 27
% sparse= bool. True= computes the sparse version
% exact= bool. True= computes the exact version, False= computes the approximated version with finite differences
% h= increment for the approximated version (if exact=true put h=0)

%%%%% PROBLEMA: matrice non sparsa

n=length(x);
s=sum(x.^2);
HF=4*(x*x');
HF(1:n+1:end)=(2/100000 -1+ 4*s + 8*x.^2)/2;
if ~exact %approximation with finite difference (not exact)
    HF= HF+2*(h^2);
end
if sparse
    HF=sparse(HF); %  non risolve il problema CAMBIA!!!!!!!!!!!!! 
end

end
